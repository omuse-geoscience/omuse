from __future__ import print_function
from amuse.support.exceptions import AmuseException
from amuse.units import units
from amuse.units.quantities import is_quantity

try:
    import scipy.ndimage as ndimage
    import scipy.interpolate as interpolate
    import scipy.spatial as spatial
except:
    raise Exception("EddyTracker requires scipy to be installed")

try:
    from matplotlib import pyplot
except:
    raise Exception("EddyTracker requires matplotlib to be installed")

try:
    from mpl_toolkits.basemap import Basemap, maskoceans
except:
    raise Exception("EddyTracker requires basemap to be installed")

import numpy

from .py_eddy_tracker_classes import plt, np, dt, Dataset, time, \
                                    datestr2datetime, gaussian_resolution, \
                                    get_cax, collection_loop, track_eddies, \
                                    anim_figure, pcol_2dxy
from .py_eddy_tracker_property_classes import SwirlSpeed
from .make_eddy_track_AVISO import *



class EddyTracker(object):

    def __init__(self, grid=None, lon=None, lat=None, domain='Regional',
                 lonmin=0., lonmax=50., latmin=-45., latmax=-20., days_between=7):

        #initialize using either lon lat or StructuredGrid passed
        if grid is not None:
            if type(grid) is not StructuredGrid:
                raise Exception("Error: grid passed to EddyTracker should be of type StructuredGrid")
            lon = grid.lon
            lat = grid.lat

        if is_quantity(lon):
            lon = numpy.array(lon.value_in(units.deg)).T
        if is_quantity(lat):
            lat = numpy.array(lat.value_in(units.deg)).T

        #force lon range into [-360.0, 0] as required by eddy tracker
        if lon[:,0].max() > 0.0:
            lon = lon - 360.0

        #perform the initialization phase and create eddy tracking list objects
        import os
        THIS_DIR = os.path.dirname(os.path.realpath(__file__))

        self.RW_PATH = THIS_DIR + '/rossrad.dat'
        self.SAVE_DIR = THIS_DIR + '/'
        self.DATA_DIR = THIS_DIR + '/'
        self.DIAGNOSTIC_TYPE = 'SLA'
        self.A_SAVEFILE = self.SAVE_DIR + 'eddy_tracks_SLA_OMUSE_anticyclonic.nc'
        self.C_SAVEFILE = self.SAVE_DIR + 'eddy_tracks_SLA_OMUSE_cyclonic.nc'

        config = {}
        self.config = config 
        config['DATA_DIR'] = self.DATA_DIR
        config['THE_DOMAIN'] = domain
        config['LONMIN'] = lonmin
        config['LONMAX'] = lonmax
        config['LATMIN'] = latmin
        config['LATMAX'] = latmax

        if len(lon.shape) != 2:
            raise Exception("Expecting longitude and latitude arrays to be 2 dimensional")
        if lon.shape != lat.shape:
            raise Exception("Expecting longitude and latitude arrays to be of same shape")

        print('lon.shape', lon.shape)

        self.grd = grd = GenericGrid(lon.shape[1], lon.shape[0], lons=lon, lats=lat, sla=None, **config)

        self._mean_ssh = None
        self.first_record = True

        self.days_between = days_between
        search_ellipse = eddy_tracker.SearchEllipse(grd.THE_DOMAIN,
                                            grd, days_between,
                                            self.RW_PATH)

        config['TRACK_DURATION_MIN'] = 0
        config['DAYS_BTWN_RECORDS'] = days_between
        config['MAX_LOCAL_EXTREMA'] = 1 # Mason et al use 1, Chelton has unlimited
        config['SAVE_FIGURES'] = False
        config['VERBOSE'] = False
        config['AMPMAX'] = 150.0 #150 is default
        config['AMPMIN'] = 1.0 #1 by default, but Mason seems to be using 0.02 for AVISO and ROMS

        self.CONTOUR_PARAMETER = np.arange(-100., 101, 1)  #used to be -100, 101
        config['CONTOUR_PARAMETER'] = self.CONTOUR_PARAMETER
        config['SHAPE_ERROR'] = np.full(self.CONTOUR_PARAMETER.size, 55.)   #Ben: Mason et al use 55 

        # Initialise two eddy objects to hold data
        self.A_eddy = A_eddy = eddy_tracker.TrackList('Anticyclonic', self.A_SAVEFILE,
                                grd, search_ellipse, **config)
        self.C_eddy = C_eddy = eddy_tracker.TrackList('Cyclonic', self.C_SAVEFILE,
                                grd, search_ellipse, **config)

        A_eddy.search_ellipse = search_ellipse
        C_eddy.search_ellipse = search_ellipse

        #constants that define the minimal and maximal eddy radius in degrees
        config['RADMIN'] = 0.35 # 0.35 is default
        config['RADMAX'] = 4.461

        #config['RADMAX'] = 8.0   #Ben testing with bigger radius
        #to see if big anticyclonic eddies in the Aghulas leakage can be found
        #result: increasing maximum radius does not help

        # See Chelton section B2 (0.4 degree radius)
        # These should give 8 and 1000 for 0.25 deg resolution
        PIXMIN = np.round((np.pi * config['RADMIN'] ** 2) /
                   grd.get_resolution() ** 2)
        PIXMAX = np.round((np.pi * config['RADMAX'] ** 2) /
                   grd.get_resolution() ** 2)
        print('--- Pixel range = %s-%s' % (np.int(PIXMIN),
                                   np.int(PIXMAX)))

        print("resolution", grd.get_resolution())

        A_eddy.PIXEL_THRESHOLD = [PIXMIN, PIXMAX]
        C_eddy.PIXEL_THRESHOLD = [PIXMIN, PIXMAX]

        # Create nc files for saving of eddy tracks
        A_eddy.create_netcdf(self.DATA_DIR, self.A_SAVEFILE)
        C_eddy.create_netcdf(self.DATA_DIR, self.C_SAVEFILE)

        # Get parameters for smoothing SLA field
        ZWL = numpy.atleast_1d(20.) # degrees, zonal wavelength (see Chelton etal 2011)
        MWL = numpy.atleast_1d(10.) # degrees, meridional wavelength
        self.ZRES, self.MRES = gaussian_resolution(grd.get_resolution(),
                                 ZWL, MWL)

        #create figure and axes for contour computation
        contfig = plt.figure(99)
        self.ax = contfig.add_subplot(111)

        #create figure, axes, and colorbar for output plot
        animfig = plt.figure(999)
        self.animax = animfig.add_subplot(111)
        # Colorbar axis
        self.animax_cbar = get_cax(self.animax, dx=0.03,
                       width=.05, position='b')

    def find_eddies(self, ssh=None, sla=None, rtime=0.0):
        #if ssh is passed the mean ssh is subtracted from the ssh
        #if the mean_ssh is currently unknown it is read in and
        #interpolated to passed grid positions
        #sla is used as passed
        grd = self.grd
        A_eddy = self.A_eddy
        C_eddy = self.C_eddy

        #rtime is a float that denotes the number of days since 01-01-0001 plus one
        #as per matplotlib.dates.date2num See:
        #http://matplotlib.org/api/dates_api.html#matplotlib.dates.date2num
        if is_quantity(rtime):
            rtime = rtime.value_in(units.day)

        #if both sla and ssh are supplied use sla
        if (ssh is not None) and (sla is not None):
            ssh=None

        #convert to cm
        if sla is not None:
            if is_quantity(sla):
                sla = numpy.array(sla.value_in(units.cm)).T
        if ssh is not None:
            if is_quantity(ssh):
                ssh = numpy.array(ssh.value_in(units.cm)).T

        #convert ssh to sla
        if ssh is not None:
            if self._mean_ssh is None:
                self._mean_ssh = grd.get_interpolated_mean_ssh()
            sla = ssh - self._mean_ssh


        #reduce sla to the area of interest (includes padding)
        sla = sla[grd.jp0:grd.jp1, grd.ip0:grd.ip1]

        # so far the preprocessing of parameters, eddy tracking starts here

        # Holding variables
        A_eddy.reset_holding_variables()
        C_eddy.reset_holding_variables()

        # Apply Gaussian smoothing
        sla -= ndimage.gaussian_filter(sla, [self.MRES, self.ZRES])

        #print "gaussian parameters", [self.MRES, self.ZRES]
        #sla = ndimage.gaussian_filter(sla, 2) #Ben: attempt to smooth the field a bit, does not help at all

        # Apply the landmask
        sla = np.ma.masked_where(grd.mask == 0, sla)

        #print "sla max min", sla.max(), sla.min()

        # Multiply by 0.01 for m
        grd.set_geostrophic_velocity(sla * 0.01)
        # Calculate EKE
        grd.getEKE()
        # Get scalar speed
        uspd = np.sqrt(grd.u ** 2 + grd.v ** 2)
        uspd = np.ma.masked_where(
            grd.mask[grd.jup0:grd.jup1,grd.iup0:grd.iup1] == 0,
            uspd)

        # Remove padded boundary
        sla = sla[grd.jup0:grd.jup1, grd.iup0:grd.iup1]

        try:
            A_eddy.sla[:] = sla.copy()
            C_eddy.sla[:] = sla.copy()
            A_eddy.slacopy[:] = sla.copy()
            C_eddy.slacopy[:] = sla.copy()
            A_eddy.uspd[:] = uspd.copy()
            C_eddy.uspd[:] = uspd.copy()
        except Exception:
            A_eddy.sla = sla.copy()
            C_eddy.sla = sla.copy()
            A_eddy.slacopy = sla.copy()
            C_eddy.slacopy = sla.copy()
            A_eddy.uspd = uspd.copy()
            C_eddy.uspd = uspd.copy()

        # Set interpolation coefficients
        grd.set_interp_coeffs(sla, uspd)
        A_eddy.sla_coeffs = grd.sla_coeffs
        A_eddy.uspd_coeffs = grd.uspd_coeffs
        C_eddy.sla_coeffs = grd.sla_coeffs
        C_eddy.uspd_coeffs = grd.uspd_coeffs

        # Get contours of Q/sla parameter
        if self.first_record:
            print('------ processing SLA contours for eddies')

        A_CS = self.ax.contour(grd.lon(),
                  grd.lat(),
                  A_eddy.sla, A_eddy.CONTOUR_PARAMETER)
        # Note that C_CS is in reverse order
        C_CS = self.ax.contour(grd.lon(),
                  grd.lat(),
                  C_eddy.sla, C_eddy.CONTOUR_PARAMETER)

        # clear the current axis
        self.ax.cla()

        # Set contour coordinates and indices for calculation of
        # speed-based radius
        A_eddy.swirl = SwirlSpeed(A_CS)
        C_eddy.swirl = SwirlSpeed(C_CS)

        # Now we loop over the CS collection
        A_eddy = collection_loop(A_CS, grd, rtime,
                         A_list_obj=A_eddy, C_list_obj=None,
                         sign_type=A_eddy.SIGN_TYPE,
                         VERBOSE=A_eddy.VERBOSE)
        # Note that C_CS is reverse order
        C_eddy = collection_loop(C_CS, grd, rtime,
                         A_list_obj=None, C_list_obj=C_eddy,
                         sign_type=C_eddy.SIGN_TYPE,
                         VERBOSE=C_eddy.VERBOSE)

        if self.first_record:
            # Set old variables equal to new variables
            A_eddy.set_old_variables()
            C_eddy.set_old_variables()

        # Track the eddies
        A_eddy = track_eddies(A_eddy, self.first_record)
        C_eddy = track_eddies(C_eddy, self.first_record)

        # Save inactive eddies to nc file
        if not self.first_record:
            if A_eddy.VERBOSE:
                print(('--- saving to nc', A_eddy.SAVE_DIR))
                print(('--- saving to nc', C_eddy.SAVE_DIR))
                print('+++')

            A_eddy.write2netcdf(rtime)
            C_eddy.write2netcdf(rtime)

        # mark the end of the first record
        self.first_record = False



    def rtime_to_ymdstr(self, rtime):
        # Get timing
        try:
            thedate = dt.num2date(rtime)[0]
        except:
            thedate = dt.num2date(rtime)
        yr = thedate.year
        mo = thedate.month
        da = thedate.day
        ymd_str = ''.join((str(yr), str(mo).zfill(2), str(da).zfill(2)))

        return ymd_str


    def plot_eddies(self, rtime=0.0):
        if is_quantity(rtime):
            rtime = rtime.value_in(units.day)

        # Get timing
        ymd_str = self.rtime_to_ymdstr(rtime)

        # Set coordinates for figures
        grd = self.grd
        Mx, My = grd.M(grd.lon(), grd.lat())

        MMx, MMy = Mx, My

        anim_figure(self.A_eddy, self.C_eddy, Mx, My, MMx, MMy, plt.cm.RdBu_r, rtime,
                self.DIAGNOSTIC_TYPE, self.SAVE_DIR, 'SLA ' + ymd_str,
                self.animax, self.animax_cbar)


    def plot_result(self, rtime=0.0):
        if is_quantity(rtime):
            rtime = rtime.value_in(units.day)

        # Get timing
        ymd_str = self.rtime_to_ymdstr(rtime)

        # Set coordinates for figures
        grd = self.grd
        Mx, My = grd.M(grd.lon(), grd.lat())

        MMx, MMy = Mx, My

        anim_figure(self.A_eddy, self.C_eddy, Mx, My, MMx, MMy, plt.cm.RdBu_r, rtime,
                self.DIAGNOSTIC_TYPE, self.SAVE_DIR, 'ALL ' + ymd_str,
                self.animax, self.animax_cbar, track_length=7/self.days_between, plot_all=True) #plot all tracks of at least 28 days




    def stop(self, rtime=0.0):
        if is_quantity(rtime):
            rtime = rtime.value_in(units.day)

        self.A_eddy.kill_all_tracks()
        self.C_eddy.kill_all_tracks()

        self.A_eddy.write2netcdf(rtime, stopper=1)
        self.C_eddy.write2netcdf(rtime, stopper=1)

        print('Outputs saved to', self.SAVE_DIR)



class GenericGrid (PyEddyTracker):
    """
    Class to satisfy the need of the eddy tracker
    to have a grid class
    """
    def __init__(self, xsize, ysize, lons=None, lats=None, sla=None, THE_DOMAIN='Global',
                 LONMIN=0.0, LONMAX=50., LATMIN=-45., LATMAX=-20., with_pad=True, DATA_DIR=''):
        """
        Initialise the grid object
        """
        super(GenericGrid, self).__init__()

        self.PRODUCT = 'Generic'
        self.DATA_DIR = DATA_DIR
        self.THE_DOMAIN = THE_DOMAIN
        self.LONMIN = LONMIN
        self.LONMAX = LONMAX
        self.LATMIN = LATMIN
        self.LATMAX = LATMAX
        self.FILLVAL = 0.0

        dims = [xsize, ysize]

        #see if lon and lat are known otherwise make something up
        if lons == None or lats == None:
            ind = numpy.indices(dims)
            X = (ind[1]/float(dims[1]) -0.5)*2*numpy.pi
            Y = (ind[0]/float(dims[0]) -0.5) * numpy.pi
            X = numpy.degrees(X)
            Y = 0.001 + numpy.degrees(Y)

        if lons == None:
            self._lon = X
        else:
            self._lon = lons

        if lats == None:
            self._lat = Y
        else:
            self._lat = lats

        if sla != None:
            self._sla = sla

        self.sla_coeffs = None
        self.uspd_coeffs = None
        self.set_initial_indices()
        self.set_index_padding()
        self.get_AVISO_f_pm_pn()
        self.set_basemap(with_pad=with_pad)
        self.set_u_v_eke()
        self.shape = self.lon().shape

        #use basemap.maskoceans to generate a mask
        mask = numpy.ones_like(self._lon)
        mask = maskoceans(self._lon, self._lat, mask, inlands=False,
                                          resolution='f', grid=1.25)
        mask = numpy.ones_like(self._lon) - mask
        self.mask = numpy.array(mask[self.jp0:self.jp1, self.ip0:self.ip1],
                                                         dtype=int)
        #generate mask at U and V points
        self.uvmask()

    #methods copied from AvisoGrid
    def lon(self):
        if self.ZERO_CROSSING:
            # TO DO: These concatenations are possibly expensive, they
            # shouldn't need to happen with every call to self.lon()
            lon0 = self._lon[self.j0:self.j1, self.i1:]
            lon1 = self._lon[self.j0:self.j1, :self.i0]
            return np.concatenate((lon0 - 360., lon1), axis=1)
        else:
            return self._lon[self.j0:self.j1, self.i0:self.i1]

    def lat(self):
        if self.ZERO_CROSSING:
            lat0 = self._lat[self.j0:self.j1, self.i1:]
            lat1 = self._lat[self.j0:self.j1, :self.i0]
            return np.concatenate((lat0, lat1), axis=1)
        else:
            return self._lat[self.j0:self.j1, self.i0:self.i1]

    def lonpad(self):
        if self.ZERO_CROSSING:
            lon0 = self._lon[self.jp0:self.jp1, self.ip1:]
            lon1 = self._lon[self.jp0:self.jp1, :self.ip0]
            return np.concatenate((lon0 - 360., lon1), axis=1)
        else:
            return self._lon[self.jp0:self.jp1, self.ip0:self.ip1]

    def latpad(self):
        if self.ZERO_CROSSING:
            lat0 = self._lat[self.jp0:self.jp1, self.ip1:]
            lat1 = self._lat[self.jp0:self.jp1, :self.ip0]
            return np.concatenate((lat0, lat1), axis=1)
        else:
            return self._lat[self.jp0:self.jp1, self.ip0:self.ip1]

    def get_resolution(self):
        #there may be a way to do this more efficiently, but after the first
        #access the answer is cached anyway
        #problem with the original code was that it could not handle gaps
        #in the domain, for example due to land blocks that have been masked

        if not hasattr(self, "_resolution"):

            lon = self.lon()[:]
            print(lon.shape)
            count = 0
            sum = 0.0
            for i in range(len(lon)):
                for j in range(len(lon[0])-1):
                    lon_p0 = lon[i,j]
                    lon_p1 = lon[i,j+1]
                    if lon_p0 != 0.0 and lon_p1 != 0.0:
                        sum += lon_p1 - lon_p0
                        count += 1
            diff_lon = sum / float(count)

            lat = self.lat()[:]
            count = 0
            sum = 0.0
            for i in range(len(lat)-1):
                for j in range(len(lat[0])):
                    lat_p0 = lat[i,j]
                    lat_p1 = lat[i+1,j]
                    if lat_p0 != 0.0 and lat_p1 != 0.0:
                        sum += lat_p1 - lat_p0
                        count += 1
            diff_lat = sum / float(count)

            self._resolution = (diff_lon + diff_lat) / 2.0

        #original code computes a geometric mean of diff_lon and diff_lat
        #but I have no clue why that would be necessary
        #return np.sqrt(np.absolute(diff_lon) * np.absolute(diff_lat)).mean()

        return self._resolution


    def umask(self):  # Mask at U points
        return self._umask

    def vmask(self):  # Mask at V points
        return self._vmask

    def f(self):  # Coriolis
        return self._f

    def gof(self):  # Gravity / Coriolis
        return self._gof

    def dx(self):  # Grid spacing along X direction
        return self._dx

    def dy(self):  # Grid spacing along Y direction
        return self._dy

    def pm(self):  # Reciprocal of dx
        return self._pm

    def pn(self):  # Reciprocal of dy
        return self._pn

    def sla(self): # SLA
        return self._sla[self.jp0:self.jp1, self.ip0:self.ip1]




    def get_interpolated_mean_ssh(self):
        given_lon = self._lon
        given_lat = self._lat

        #mean_ssh = numpy.fromfile(DATA_DIR + "mean_ssh.dat")
        #mean_lon = numpy.fromfile(DATA_DIR + "mean_ssh_lon.dat")
        #mean_lat = numpy.fromfile(DATA_DIR + "mean_ssh_lat.dat")
        #mean_ssh = mean_ssh.reshape(mean_lon.size, mean_lat.size, order='F')
        #mean_ssh *= 100.0 #m to cm

        from netCDF4 import Dataset
        filename = self.DATA_DIR + "SSH_t.t0.1_42l_nccs01.avg0315-0325.nc"
        f = Dataset(filename, "r")
        mean_lon = numpy.array(f.variables['t_lon'][:])
        mean_lat = numpy.array(f.variables['t_lat'][:])
        mean_ssh = numpy.array(f.variables['SSH'][:]).T

        fillvalue = f.variables['SSH']._FillValue
        mean_ssh[mean_ssh == fillvalue] = 0.0

        #print 'mean_ssh.shape', mean_ssh.shape

        #print 'mean_lon.shape', mean_lon.shape
        #print 'mean_lon', mean_lon
        #print 'given_lon.shape', given_lon.shape
        #print 'given_lon', given_lon

        #print 'mean_lat.shape', mean_lat.shape
        #print 'mean_lat', mean_lat
        #print 'given_lat.shape', given_lat.shape
        #print 'given_lat', given_lat

        interp_func = interpolate.RectBivariateSpline(mean_lon,
                                                    mean_lat, mean_ssh)

        interpolated_mean_ssh = numpy.zeros(given_lon.shape)

        #doing this one at a time because of memory error
        for i in range(len(given_lon)):
            for j in range(len(given_lon[0])):
                lon = given_lon[i,j]
                if lon < 0.0:
                    lon += 360.0
                if lon > 360.0:
                    lon -= 360.0
                lat = given_lat[i,j]

                if lon == 0.0 and lat == 0.0:
                    interpolated_mean_ssh[i,j] = 0.0
                else:
                    interpolated_mean_ssh[i,j] = interp_func(lon, lat)

        return interpolated_mean_ssh




#go interactive
#import readline
#import rlcompleter
#readline.parse_and_bind("tab: complete")
#import code
#code.interact(local=dict(globals(), **locals()) )


if __name__ == "__main__":

    dims = [1440, 720]

    # generate some fake SLA field
    sla = numpy.random.random(dims) * 100.0 #for m to cm
    #smooth for slightly more realistic field
    sla = ndimage.gaussian_filter(sla, 2)

    grd = GenericGrid(dims[0], dims[1], sla=sla, **config)

    find_eddies(grd)
