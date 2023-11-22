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
    invariate_variables : list of str, optional
        invariate variables to load. Will be included on the grid (but not updated). default: ["land_sea_mask"]
    """

    def __init__(self, source_dir, source_file, start_datetime=datetime.datetime(1979, 1, 2), 
        variables=[], invariate_variables = []):

        self.source_path = os.path.join(source_dir, source_file)

        self.variables=variables
        self.start_datetime=start_datetime
        self.tnow=0. | units.day

        self.invariate_variables=invariate_variables

        self.dataset = self.get_dataset()

        self.grid=self.generate_initial_grid()

        self.update_grid()
        
        super(Data, self).__init__()


    def generate_initial_grid(self):
        grid=None
        data = self.dataset
        for variable in self.variables:

            shortname = variable
            print("dataset = ", data)    

            print("shortname = ", shortname)

            lat=data["latitude"][:]
            lon=data["longitude"][:]
            
            dx=float(lon[1]-lon[0])
            dy=float(lat[1]-lat[0])
            
            assert lon[1]>lon[0]                
            assert lat[1]<lat[0]                
            assert dx==-dy
            self.shape=data[shortname][0,:,:].shape
            print("shape = ", self.shape)
            print("dx = ", dx)
            print("dy = ", dy|units.deg)    
            if grid is None:
                grid=new_cartesian_grid(self.shape, 
                                        dx|units.deg,
                                        offset=([lat[-1]-dx/2,lon[0]-dx/2]|units.deg),
                                        axes_names=["lat","lon"])
 
            value=numpy.zeros(self.shape)
            
            setattr(grid, variable, value)

        return grid

    def evolve_model(self, endtime):
        self.tnow=endtime
        self.update_grid()

    def update_grid(self):
        for v in self.variables:
            self.update_variable(v)

    def update_variable(self, var):
        time=self.start_datetime+datetime.timedelta(days=self.tnow.value_in(units.day))
        print(time)
        shortname = var
        dataset=self.dataset
        print("Dataset = \n", dataset)
        index = 15
        _value=dataset[shortname][index, ...]
        print("_value = \n", _value)
        print("self = \n", self)
        assert len(_value.shape)==len(self.shape)

        value=numpy.zeros(self.shape)    
            
        setattr(self.grid, "_"+var, value)

    @property
    def model_time(self):
        return self.tnow
  
    def get_dataset(self):
        if not os.path.isfile(self.source_path):
            raise Exception("this file does not exist")

        dataset = xr.open_dataset(self.source_path)
        return dataset
        

datacached=Data
