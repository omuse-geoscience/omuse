from __future__ import print_function
#!/usr/bin/env python

import scipy.ndimage as ndimage
import scipy.interpolate as interpolate
import scipy.spatial as spatial

from .py_eddy_tracker_classes import plt, np, dt, Dataset, time, \
                                    datestr2datetime, gaussian_resolution, \
                                    get_cax, collection_loop, track_eddies, \
                                    anim_figure, pcol_2dxy
from .py_eddy_tracker_property_classes import SwirlSpeed

from mpl_toolkits.basemap import Basemap, maskoceans

import numpy
from matplotlib import pyplot

from .make_eddy_track_AVISO import *


#config = {}
#config['THE_DOMAIN'] = 'Regional'
#config['LONMIN'] = 0.
#config['LONMAX'] = 50.
#config['LATMIN'] = -45.
#config['LATMAX'] = -25.

import os
THIS_DIR = os.path.dirname(os.path.realpath(__file__))

RW_PATH = THIS_DIR + '/rossrad.dat'
SAVE_DIR = THIS_DIR + '/'
DATA_DIR = THIS_DIR + '/'
DIAGNOSTIC_TYPE = 'SLA'
A_SAVEFILE = "".join([SAVE_DIR, 'eddy_tracks_SLA_OMUSE_anticyclonic.nc'])
C_SAVEFILE = "".join([SAVE_DIR, 'eddy_tracks_SLA_OMUSE_cyclonic.nc'])

class GenericGrid (PyEddyTracker):
    """
    Class to satisfy the need of the eddy tracker
    to have a grid class
    """
    def __init__(self, xsize, ysize, lons=None, lats=None, sla=None, THE_DOMAIN='Global',
                 LONMIN=0.0, LONMAX=360, LATMIN=-85, LATMAX=85, with_pad=True):
        """
        Initialise the grid object
        """
        super(GenericGrid, self).__init__()

        self.PRODUCT = 'Generic'
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
        #pyplot.pcolormesh(self._lon, self._lat, mask, cmap=pyplot.cm.bone)
        #pyplot.show()
        #raw_input()

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
        #there may be was to do this more efficiently, but after the first
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

            print("resolution=", self._resolution)

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



def find_eddies(grd):
    X = grd.lon()
    Y = grd.lat()
    sla = grd.sla()

    search_ellipse = eddy_tracker.SearchEllipse(grd.THE_DOMAIN,
                                            grd, 7,
                                            RW_PATH)

    # Initialise two eddy objects to hold data
    A_eddy = eddy_tracker.TrackList('Anticyclonic', A_SAVEFILE,
                                grd, search_ellipse, **config)
    C_eddy = eddy_tracker.TrackList('Cyclonic', C_SAVEFILE,
                                grd, search_ellipse, **config)

    A_eddy.search_ellipse = search_ellipse
    C_eddy.search_ellipse = search_ellipse

    config['RADMIN'] = 0.35
    config['RADMAX'] = 4.461

    # See Chelton section B2 (0.4 degree radius)
    # These should give 8 and 1000 for 0.25 deg resolution
    PIXMIN = np.round((np.pi * config['RADMIN'] ** 2) /
                   grd.get_resolution() ** 2)
    PIXMAX = np.round((np.pi * config['RADMAX'] ** 2) /
                   grd.get_resolution() ** 2)
    print('--- Pixel range = %s-%s' % (np.int(PIXMIN),
                                   np.int(PIXMAX)))

    A_eddy.PIXEL_THRESHOLD = [PIXMIN, PIXMAX]
    C_eddy.PIXEL_THRESHOLD = [PIXMIN, PIXMAX]

    # Create nc files for saving of eddy tracks
    A_eddy.create_netcdf(DATA_DIR, A_SAVEFILE)
    C_eddy.create_netcdf(DATA_DIR, C_SAVEFILE)

    # Holding variables
    A_eddy.reset_holding_variables()
    C_eddy.reset_holding_variables()



    ZWL = numpy.atleast_1d(20.) # degrees, zonal wavelength (see Chelton etal 2011)
    MWL = numpy.atleast_1d(10.) # degrees, meridional wavelength
    ZRES, MRES = gaussian_resolution(grd.get_resolution(),
                                 ZWL, MWL)

    #end of initialization

    # Apply Gaussian smoothing
    sla -= ndimage.gaussian_filter(sla, [MRES, ZRES])

    # Apply the landmask
    sla = np.ma.masked_where(grd.mask == 0, sla)

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
    if 'first_record' not in locals():

        print('------ processing SLA contours for eddies')
        contfig = plt.figure(99)
        ax = contfig.add_subplot(111)

        animfig = plt.figure(999)
        animax = animfig.add_subplot(111)
        # Colorbar axis
        animax_cbar = get_cax(animax, dx=0.03,
                          width=.05, position='b')

    A_CS = ax.contour(grd.lon(),
                  grd.lat(),
                  A_eddy.sla, A_eddy.CONTOUR_PARAMETER)
    # Note that C_CS is in reverse order
    C_CS = ax.contour(grd.lon(),
                  grd.lat(),
                  C_eddy.sla, C_eddy.CONTOUR_PARAMETER)

    # clear the current axis
    ax.cla()

    # Set contour coordinates and indices for calculation of
    # speed-based radius
    A_eddy.swirl = SwirlSpeed(A_CS)
    C_eddy.swirl = SwirlSpeed(C_CS)

    rtime = 0

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


    #this block should only be executed the first iteration
    first_record = True
    # Set old variables equal to new variables
    A_eddy.set_old_variables()
    C_eddy.set_old_variables()

    # Track the eddies
    A_eddy = track_eddies(A_eddy, first_record)
    C_eddy = track_eddies(C_eddy, first_record)


    # Set coordinates for figures
    Mx, My = grd.M(grd.lon(), grd.lat())

    MMx, MMy = Mx, My 
    ymd_str = '20151119'
    anim_figure(A_eddy, C_eddy, Mx, My, MMx, MMy, plt.cm.RdBu_r, rtime,
            DIAGNOSTIC_TYPE, SAVE_DIR, 'SLA ' + ymd_str,
            animax, animax_cbar)






def get_interpolated_mean_ssh(given_lon, given_lat):

    #mean_ssh = numpy.fromfile(DATA_DIR + "mean_ssh.dat")
    #mean_lon = numpy.fromfile(DATA_DIR + "mean_ssh_lon.dat")
    #mean_lat = numpy.fromfile(DATA_DIR + "mean_ssh_lat.dat")
    #mean_ssh = mean_ssh.reshape(mean_lon.size, mean_lat.size, order='F')
    #mean_ssh *= 100.0 #m to cm

    from netCDF4 import Dataset
    filename = DATA_DIR + "SSH_t.t0.1_42l_nccs01.avg0315-0325.nc"
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
