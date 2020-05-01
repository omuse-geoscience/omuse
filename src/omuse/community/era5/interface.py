import os
import hashlib
import numpy
import datetime

from omuse.units import units
from omuse.units import quantities
from amuse.datamodel import new_cartesian_grid

from . import era5

from netCDF4 import Dataset

_era5_units_to_omuse={ 
                       "~" : units.none, 
                       "dimensionless" : units.none, 
                       "(0-1)" : units.none, 
                       "unknown" : units.none, 
                       "none" : units.none, 
                       "Proportion" : units.none, 
                       "code table (4.201)": units.none,
                       "s" : units.s, 
                       "K" : units.K, 
                       "m" : units.m, 
                       "m s**-1" : units.m/units.s,
                       "W m**-2" : units.W/units.m**2,
                       "W m**-1" : units.W/units.m**1,
                       "J m**-2" : units.J/units.m**2,
                       "J m**-1" : units.J/units.m**1,
                       "m of water equivalent": units.m, # make named unit
                       "Pa" : units.Pa,
                       "kg m**-2 s**-1" : units.kg/units.m**2/units.s,
                       "kg m**-2" : units.kg/units.m**2,
                       "kg m**-3" : units.kg/units.m**3,
                       "K kg m**-2" : units.K*units.kg/units.m**2,
                       "kg m**-1 s**-1" : units.kg/units.m/units.s,
                       "m**2 m**-2" : units.m**2/units.m**2,
                       "m**3 m**-3" : units.m**3/units.m**3,
                       "degrees" : units.deg,
                       "radians" : units.rad,
                       "N m**-2" : units.N/units.m**2,
                       "N m**-2 s" : units.N*units.s/units.m**2,
                       "m**-1" : units.m**-1,
                       "J kg**-1" : units.J/units.kg,
                       }

class ERA5(object):

    def __init__(self, maxcache=None, cachedir="./__era5_cache", 
                 start_datetime=datetime.datetime(1979,1,2), variables=[], 
                 grid_resolution=None, nwse_boundingbox=None):
        self.maxcache=maxcache
        if not os.path.exists(cachedir):
            os.mkdir(cachedir)
        self.cachedir=cachedir
        self.download_timespan="day"
        self.nwse_boundingbox=nwse_boundingbox

        self.variables=variables
        self.start_datetime=start_datetime
        self.tnow=0. | units.day

        if grid_resolution is None:
          self.grid_resolution=0.25 | units.deg

        self.grid_resolution=grid_resolution
        
        self.grid=self.generate_initial_grid()

        self.update_grid()

    def generate_initial_grid(self):
        lsm=self.get_dataset("land_sea_mask", 
                             nwse_boundingbox=self.nwse_boundingbox,
                             grid_resolution=self.grid_resolution)
        
        lat=lsm["latitude"][:]
        lon=lsm["longitude"][:]
        
        dx=lat[1]-lat[0]
        dy=lat[1]-lat[0]
        
        assert dx==dy
        
        shortname=era5.SHORTNAME["land_sea_mask"]
        grid=new_cartesian_grid( lsm[shortname][0,:,:].shape, dx | units.deg, 
                                       offset=[lat[0],lon[0]]|units.deg, axes_names=["lat","lon"])
        grid.land_sea_mask=lsm[shortname][0,:,:]

        return grid


    def evolve_model(self, endtime):
        self.tnow=endtime
        self.update_grid()

    def update_grid(self):
        for v in self.variables:
            self.update_variable(v)
            
    def update_variable(self, var):
      
        time=self.start_datetime+datetime.timedelta(days=self.tnow.value_in(units.day))
      
        dataset=self.get_dataset(var, time=time, download_timespan=self.download_timespan,
                                  nwse_boundingbox=self.nwse_boundingbox,
                                  grid_resolution=self.grid_resolution)

        # the netcdf hortname do not match the CDS website documentation
        # little code to extract variable shortname 

        shortname=era5.SHORTNAME[var]
        notvars=['longitude', 'latitude', 'time']
        netcdf_keys=[x for x in dataset.variables.keys() if (x not in notvars)]
        if shortname not in netcdf_keys:
            for x in dataset.variables.keys():
              if var==dataset.variables[x].long_name.lower().strip().replace(" ","_"):
                shortname=x              
        if shortname not in netcdf_keys and len(netcdf_keys)==1:
            shortname=netcdf_keys[0]
            print('new shortname for "{0}" : "{1}"'.format(var, shortname))
            era5.SHORTNAME[var]=shortname
      
        if self.download_timespan=="day":
            value=dataset[shortname][time.hour,...]
        else:
            value=dataset[shortname][0,...]

        value=value | _era5_units_to_omuse[era5.UNITS[var]]

        setattr(self.grid, "_"+var, value)


    @property
    def model_time(self):
        return self.tnow

    @staticmethod
    def generate_outputfile(name, request, directory="./"):
        filename="_era5_cache_"+hashlib.sha1((name+repr(request)).encode()).hexdigest() + ".nc"
        return os.path.join(directory, filename)
  
    def maintain_cache(self, datafile):
        pass
  
    def get_dataset(self, variable="land_sea_mask", time=datetime.datetime(1979,1,1), 
                    download_timespan=None, grid_resolution=None, nwse_boundingbox=None):        
        if download_timespan=="day":
            hour=["%2.2i:00"%i for i in range(24)]
        else:
            hour="%2.2i:00"%(time.hour+1)
        year="%4.4i"%time.year
        month="%2.2i"%time.month
        day="%2.2i"%time.day
        
        if nwse_boundingbox is not None:
            area=quantities.value_in(nwse_boundingbox, units.deg)
        else:
            area=None
                
        if grid_resolution is None:
            g=0.25
        else:
            g=grid_resolution.value_in(units.deg)
        if g not in [0.25,0.5,1.0]:
            raise Exception("unsupported grid resolution %s"%str(g))
            
        name,request=era5.build_request(variable=variable, year=year, month=month, 
                            day=day, hour=hour, area=area, grid=[g, g])
        datafile=self.generate_outputfile(name,request, self.cachedir)
        if not os.path.isfile(datafile):
            era5.fetch(name,request, datafile)
            self.maintain_cache(datafile)
            self._append_history(name, request, datafile)
        result=Dataset(datafile)
        
        if download_timespan=="day":
            assert len(result["time"])==24 
        else:
            assert len(result["time"])==1
            
        return result
        

    def _append_history(self, name, request, filename):
        """Append download information. """
        dtime = datetime.datetime.utcnow().strftime(
              "%Y-%m-%d %H:%M:%S %Z")
        appendtxt = "ERA5 OMUSE time: {}, name: {} request: {}".format(dtime, name, request)
        _, extension = os.path.splitext(filename)
        d = Dataset(filename, 'r+')
        try:
            d.history = "{} ;; {}".format(appendtxt, d.history)
        except AttributeError:
            d.history = appendtxt
        d.close()

era5cached=ERA5
