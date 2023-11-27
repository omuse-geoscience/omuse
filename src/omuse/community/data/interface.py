import os
import hashlib
import numpy
import datetime
import calendar
import xarray as xr

from omuse.units import units
from omuse.units import quantities

from amuse.datamodel import new_cartesian_grid

from amuse.support.literature import LiteratureReferencesMixIn

from netCDF4 import Dataset

class Data(LiteratureReferencesMixIn):
    """Data interface
    
    Arguments
    ---------
    source_file: str
        path to the data source file
    source_dir: str
        path to the data source directory
    start_datetime: datetime object, optional
        starting date from which to load the data. Must be datatime object, 
        default: datetime.datetime(1979,1,2)
    variables: list of str, optional
        variables to be loaded. Must be list of valid strings. default: []
    units: list of omuse units or None, optional
        units to asociate with each variable. default: []
    invariate_variables : list of str, optional
        invariate variables to load. Will be included on the grid (but not updated). default: ["land_sea_mask"]
    center_meridian: bool, optional
        whether or not to shift the dataset by 180 degrees to center the meridian. default: True
    """

    def __init__(self, source_dir, source_file, start_datetime=datetime.datetime(1979, 1, 2), 
        variables=[], var_units = [], invariate_variables = [], center_meridian = True):

        self.source_path = os.path.join(source_dir, source_file)

        self.variables = variables
        self.units = {variables[i]: var_units[i] for i in range(len(variables))}
        self.start_datetime = start_datetime
        self.tnow=0. | units.day

        self.invariate_variables = invariate_variables
        self.center_meridian = center_meridian
        self.dataset = self.get_dataset()

        self.grid=self.generate_initial_grid()

        self.update_grid()
        
        super(Data, self).__init__()

    def generate_initial_grid(self):
        grid=None
        data = self.dataset
        for variable in self.variables:

            shortname = variable


            lat=data["latitude"][:]
            lon=data["longitude"][:]
            
            dx=float(lon[1]-lon[0])
            dy=float(lat[1]-lat[0])
            
            assert lon[1]>lon[0]                
            assert lat[1]>lat[0]                
            assert dx==dy
            self.shape=data[shortname][0,:,:].shape
   
            if grid is None:
                grid=new_cartesian_grid(self.shape, 
                                        dx|units.deg,
                                        offset=([lat[0]-dx/2,lon[0]-dx/2]|units.deg),
                                        axes_names=["lat","lon"])
 
            value=numpy.zeros(self.shape)
            
            setattr(grid, variable, value | self.units[variable])

        return grid

    def evolve_model(self, endtime):
        self.tnow=endtime
        self.update_grid()

    def update_grid(self):
        for v in self.variables:
            self.update_variable(v)

    def update_variable(self, var):
        time_now=self.start_datetime+datetime.timedelta(days=self.tnow.value_in(units.day))
        shortname = var
        dataset=self.dataset
        _value=dataset[shortname].sel(time=time_now).values
        assert len(_value.shape)==len(self.shape)
            
        setattr(self.grid, var, _value | self.units[var])

    @property
    def model_time(self):
        return self.tnow

    def shift_coordinates(self, ds):
        # Shift longitude values to enforce (-180, 180) versus eg. (0, 360) 
        # which corresponds to default model representations
        ds["_lon_temp"] = xr.where(
            ds["longitude"] > 180,
            ds["longitude"] - 360,
            ds["longitude"])
        ds["_lon_temp"] = xr.where(
            ds["longitude"] < -180,
            ds["longitude"] + 360,
            ds["longitude"])

        ds = (
            ds
            .swap_dims({"longitude": "_lon_temp"})
            .sel(**{"_lon_temp": sorted(ds._lon_temp)})
            .drop("longitude"))

        ds = ds.rename({"_lon_temp": "longitude"}) 
        return ds  

    def flip_latitude(self, ds):
        ds = ds.sortby("latitude")
        return ds

    def get_dataset(self):
        if not os.path.isfile(self.source_path):
            raise Exception("this file does not exist")

        dataset = xr.open_dataset(self.source_path)
        if self.center_meridian:
            dataset = self.shift_coordinates(dataset)
        dataset = self.flip_latitude(dataset)
        return dataset
        

datacached=Data
