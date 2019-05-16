import numpy

from netCDF4 import Dataset

from amuse.io import base
from omuse.units import units

from amuse.datamodel import Grid

celsius=units.K # placeholder

units_lookup={ 
    "cm/s" : units.cm/units.s,
    "dyne/cm2" : units.dyn/units.cm**2,
    "dyne/cm3" : units.dyn/units.cm**3,
    "deg C" : celsius, 
    "msu (g/g)" : units.g/units.g,
    "cm2/s" : units.cm**2/units.s
      }

class POPNetCDFFileFormatProcessor(base.FileFormatProcessor):
    def load(self):
      
      dataset=Dataset(self.filename)
      
      shape=()
      order=dict()
      for i,(key,d) in enumerate(dataset.dimensions.items()):
          shape+=(len(d),)
          order[key]=i

      grid=Grid(*shape[0:2])
      
      for s in dataset.ncattrs():
          setattr(grid.collection_attributes,s,getattr(dataset,s))
      
      for name,var in dataset.variables.items():
        perm=[order[x] for x in var.dimensions]
        if hasattr(var,"units"):
          unit=units_lookup[var.units]
        else:
          unit=units.none
        setattr(grid,name,numpy.array(var).transpose(tuple(perm)) | unit)
      
        
      self.dataset=dataset  
        
      return grid
    def close(self):
      self.dataset.close()
