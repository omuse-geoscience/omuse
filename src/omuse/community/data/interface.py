import os
import hashlib
import numpy
import datetime
import calendar
import xarray as xr

from omuse.units import units
from omuse.units import quantities
from omuse.community.era5 import era5

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
    nwse_crop: None or [North,West,South,East], optional
        Bounding box of grid to which to crop the original dataset, 
        will be [90,-180,-90,180] | units.deg (whole globe) if None. default: None
    era5_metadata: bool, optional
        Whether the data uses the era5 metadata format
    """

    def __init__(self, source_dir, source_file,
                 start_datetime=datetime.datetime(1979, 1, 2), variables=[], invariate_variables = [],
                 grid_resolution=None, nwse_crop=None):

        self.source_path = os.path.join(source_dir, source_file)

        self.nwse_crop=nwse_crop

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

            shortname=self._get_shortname(data,variable)
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

    def _get_shortname(self, dataset, var):
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
        return shortname


    def update_variable(self, var):
      
        time=self.start_datetime+datetime.timedelta(days=self.tnow.value_in(units.day))
        print(time)
      
        dataset=self.dataset
        print("Dataset = \n", dataset)
        shortname=self._get_shortname(dataset, var)
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

    @staticmethod
    def generate_outputfile(name, request, directory="./"):
        filename="_era5_cache_"+hashlib.sha1((name+repr(request)).encode()).hexdigest() + ".nc"
        return os.path.join(directory, filename)
  
    def get_dataset(self):
        if not os.path.isfile(self.source_path):
            raise Exception("this file does not exist")

        # The Dataset constructor has a "parallel" keyword, this uses mpi4py
        # would be good to investigate whether it would help to load this dataset
        # using the parallel keyword 
        # It should also be considered to move to xarray to load in the data

        # dataset = Dataset(self.source_path, mode="r")

        dataset = xr.open_dataset(self.source_path)
        return dataset
        

datacached=Data
