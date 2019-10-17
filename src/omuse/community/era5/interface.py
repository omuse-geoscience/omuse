import os
import sha
import numpy
import datetime

from omuse.units import units
from omuse.units import quantities
from amuse.datamodel import new_cartesian_grid

import era5

from netCDF4 import Dataset

_era5_units_to_omuse={ "none" : units.none, "K" : units.K, "m" : units.m}

class era5cached(object):

    def __init__(self, maxcache=None, cachedir="./", 
                 start_datetime=datetime.datetime(1979,1,2), variables=[], nwse_boundingbox=None):
        self.maxcache=maxcache
        self.cachedir=cachedir
        self.download_timespan="day"
        self.nwse_boundingbox=nwse_boundingbox

        self.variables=variables
        self.start_datetime=start_datetime
        self.tnow=0. | units.day

        self.grid=self.generate_initial_grid()

        self.update_grid()

    def generate_initial_grid(self):
        lsm=self.get_dataset("land_sea_mask", nwse_boundingbox=self.nwse_boundingbox)
        
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
                                  nwse_boundingbox=self.nwse_boundingbox)
      
        if self.download_timespan=="day":
            value=dataset[era5.SHORTNAME[var]][time.hour,...]
        else:
            value=dataset[era5.SHORTNAME[var]][0,...]

        value=value | _era5_units_to_omuse[era5.UNITS[var]]

        setattr(self.grid, "_"+var, value)


    @property
    def model_time(self):
        return self.tnow

    @staticmethod
    def generate_outputfile(name, request, directory="./"):
        filename="_era5_cache_"+sha.sha(name+repr(request)).hexdigest() + ".nc"
        return os.path.join(directory, filename)
  
    def maintain_cache(self, datafile):
        pass
  
    def get_dataset(self, variable="land_sea_mask", time=datetime.datetime(1979,1,1), 
                    download_timespan=None, nwse_boundingbox=None):        
        if download_timespan=="day":
            hour=["%2.2i:00"%i for i in range(24)]
        else:
            hour="%2.2i:00"%(time.hour+1)
        year="%4.4i"%time.year
        month="%2.2i"%time.month
        day="%2.2i"%time.day
        
        if nwse_boundingbox is not None:
            area=quantities.value_in(nwse_boundingbox, units.deg)
        
        name,request=era5.build_request(variable=variable, year=year, month=month, 
                            day=day, hour=hour, area=area)
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

if __name__=="__main__":
  
    e=era5cached(variables=["2m_temperature", "total_precipitation"], nwse_boundingbox=[40, 0, 0, 60]| units.deg)
    print(e.grid) # note grid has prepended the names with _ (because long names are not valid python var names)

    t1=e.grid._2m_temperature.copy()
    e.evolve_model(24. | units.hour)
    t2=e.grid._2m_temperature
        
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib import pyplot
    #~ pyplot.imshow(e.grid.land_sea_mask)
    pyplot.imshow((t2-t1).value_in(units.K))
    pyplot.show()
