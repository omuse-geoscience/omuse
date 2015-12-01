import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from omuse.community.swan.interface import SwanInterface, Swan

from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

#default_options={}
default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

class TestSwanInterface(TestWithMPI):
    def test1(self):
      s=SwanInterface(**default_options)
      err=s.initialize_code()
      self.assertEqual(err,0)
      err=s.cleanup_code()
      self.assertEqual(err,0)

    def test2(self):
      s=SwanInterface(**default_options)
      err=s.initialize_code()      
      self.assertEqual(err,0)
      
      x,err=s.get_wlev()
      self.assertEqual(x,0.)
      self.assertEqual(err,0)
      err=s.set_wlev(12.)
      self.assertEqual(err,0)
      x,err=s.get_wlev()
      self.assertEqual(err,0)
      self.assertEqual(x,12.)

    def test3(self):
      s=SwanInterface(**default_options)
      err=s.initialize_code()
      self.assertEqual(err,0)

      for name,default in [("calc_mode","stationary"),("grid_type","regular"),
          ("input_grid_type","regular"),("projection_method","quasi-cart."),("number_dimensions",2),
          ("coordinates","cartesian")]:
        x,err=eval("s.get_"+name+"()")
        self.assertEqual(x,default)
        self.assertEqual(err,0)

      for name,default,val in [("wlev",0.,12.),("grav",9.81,12.),
          ("rho",1025.,1234.),("cdcap",99999.,12.5e-3),("rearth",6366197.5,12345.)]:
        x,err=eval("s.get_"+name+"()")
        self.assertAlmostEqual(x,default,6)
        self.assertEqual(err,0)
        err=eval("s.set_"+name+"(val)")
        self.assertEqual(err,0)
        x,err=eval("s.get_"+name+"()")
        self.assertEqual(err,0)
        self.assertAlmostEqual(x,val,6)

    def test4(self):
      s=SwanInterface(**default_options)
      err=s.initialize_code("abc","defg","hij","klm")
      self.assertEqual(err,-1)

      for name,default in [("coordinates","abc"),("calc_mode","defg"),("grid_type","hij"),
          ("input_grid_type","klm")]:
        x,err=eval("s.get_"+name+"()")
        self.assertEqual(x,default)
        self.assertEqual(err,0)

    def test5(self):
      pass

class TestSwan(TestWithMPI):
    def test1(self):
      s=Swan(**default_options)
      print s.parameters

    def test2(self):
      options=dict(coordinates="cartesian")
      options.update(default_options)
      s=Swan(**options)
      u=s.parameters.grid_origin_x
      self.assertEqual(u,0 | units.m)

    def test3(self):
      options=dict(coordinates="spherical")
      options.update(default_options)
      s=Swan(**options)
      u=s.parameters.grid_origin_x
      v=s.parameters.grid_length_x
      self.assertEqual(u,0 | units.deg)
      self.assertEqual(v,0 | units.deg)
