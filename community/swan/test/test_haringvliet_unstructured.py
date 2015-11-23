import numpy

from omuse.community.swan.interface import SwanInterface,Swan
from omuse.units import units

from matplotlib import pyplot

from amuse.test.amusetest import TestWithMPI

from amuse.io import write_set_to_file

from read_triangle_mesh import read_triangle_mesh

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)


def read_bot_data(filename="f31hari.bot",n=88,m=117):
  f=open(filename,"r")
  lines=f.readlines()
  f.close()
  
  dat=[]
  for line in lines:
    for s in line.split():
      dat.append(float(s))
  
  return numpy.transpose(numpy.array(dat).reshape((m,n)))

class TestHaringvliet(object):
    def test2(self):
        bathymetry=read_bot_data()
        bathymetry=numpy.array(bathymetry,dtype="float32")
          
        rt=read_triangle_mesh("f32hari") 
        rt.read_grid()
        nodes,elements=rt.get_sets() 
          
        s=Swan(grid_type="unstructured", redirection="none")

        print s.parameters.coordinates
        print s.parameters.projection_method
        print s.parameters.grid_type
        print s.parameters.input_grid_type
        print s.parameters.calculation_mode
        print s.parameters.number_of_dimensions

        ncells=len(elements)
        nverts=len(nodes)
        msc=32
        mdc=36

        s.parameters.number_of_cells=ncells
        s.parameters.number_of_vertices=nverts
        s.parameters.number_of_directions=mdc
        s.parameters.number_of_frequencies=msc
        s.parameters.lowest_frequency=2*numpy.pi*0.0521 | units.rad/units.s
        s.parameters.highest_frequency=2*numpy.pi | units.rad/units.s

        s.parameters.input_grid_origin_x=0. | units.m
        s.parameters.input_grid_origin_y=0. | units.m
        s.parameters.input_grid_orientation=0. | units.deg
        s.parameters.input_grid_dx=250. | units.m
        s.parameters.input_grid_dy=250. | units.m
        s.parameters.input_grid_nmesh_x=87
        s.parameters.input_grid_nmesh_y=116
      
        s.parameters.constant_water_level=0.3 | units.m
      
        #~ print s.initialize_grids()

        s.parameters.uniform_wind_velocity=12. | units.m/units.s
        s.parameters.uniform_wind_direction=8.8 | units.deg
      
        s.parameters.unstructured_boundary_spec_file="f31har01.bnd"
        s.parameters.boundary_marker=2
      
        s.parameters.use_gen3_parameters=True
        s.parameters.use_breaking_parameters=True
        s.parameters.use_triads_parameters=True
        s.parameters.use_friction_parameters=True
        s.parameters.use_uniform_wind=True
      
        #~ s.commit_parameters()
        print s.parameters
       
        exc=s.get_exc_value(1)
        
        bathymetry[bathymetry==-99.]=exc
        print (bathymetry==-1.e20).sum(),exc
        input_shape=bathymetry.shape
        ii,jj=numpy.mgrid[1:input_shape[0]+1,1:input_shape[1]+1]

        print s.forcings
        print bathymetry.shape
        
        channel=nodes.new_channel_to(s.nodes)
        channel.copy_attributes(["x","y","vmark"])
        channel=elements.new_channel_to(s.elements)
        channel.copy_attributes(["n1","n2","n3"])

        s.forcings.depth=bathymetry | units.m
      
        #~ print s.commit_grids()
      
        #~ print s.initialize_boundary()

        s.evolve_model(0. | units.s)
      
        write_set_to_file(s.nodes,"nodes.amuse","amuse")
        write_set_to_file(s.elements,"elements.amuse","amuse")
        
if __name__=="__main__":
    test=TestHaringvliet()
    test.test2()
