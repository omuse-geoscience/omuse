import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from omuse.community.qgcm.interface import QGCMInterface, QGCM

from amuse.units import units

#default_options={}
default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

class TestQGCMInterface(TestWithMPI):
    def test1(self):
      s=QGCMInterface(**default_options)
      err=s.initialize_code()
      self.assertEqual(err,0)
      #~ err=s.cleanup_code() # cleanup needs changes
      #~ self.assertEqual(err,0)

    def test2(self):
      s=QGCMInterface(**default_options)
      err=s.initialize_code()
      self.assertEqual(err,0)
      err=s.commit_parameters()
      self.assertEqual(err,0)
      #~ err=s.cleanup_code()
      #~ self.assertEqual(err,0)


    def test3(self):
      s=QGCMInterface(**default_options)
      err=s.initialize_code()
      self.assertEqual(err,0)
      err=s.commit_parameters()
      self.assertEqual(err,0)
      err=s.commit_grids()
      self.assertEqual(err,0)      
      err=s.cleanup_code()
      self.assertEqual(err,0)

class TestQGCM(TestWithMPI):
    def test1(self):
      s=QGCM(**default_options)
      s.stop()

    def test2(self):
      s=QGCM(**default_options)
      
      s.evolve_model(1. | units.day)
      
      s.stop()

    def test3(self):
      s=QGCM(**default_options)
      
      s.evolve_model(12. | units.hour)
      s.evolve_model(24. | units.hour)
      
      s.stop()

    def test4(self):
      s=QGCM(**default_options)

      print s.ocean_P_grid
      print s.ocean_P_grid_forcings
      print s.ocean_T_grid_forcings
      print s.atmosphere_T_grid

      s.evolve_model(1. | units.hour)

      print s.ocean_P_grid.pressure[10,10,0]

class TestRestartQGCM(TestWithMPI):
    def test1(self):
      s=QGCM(**default_options)

      s.evolve_model(1. | units.hour)

      pref=s.ocean_P_grid.pressure
  
      s1=QGCM(**default_options)

      s1.evolve_model(0.5 | units.hour)
      s1.evolve_model(1. | units.hour)

      p=s1.ocean_P_grid.pressure

      self.assertEquals(p,pref)

    def test2(self):
      s=QGCM(**default_options)

      dto=s.parameters.atmosphere_timestep*s.parameters.timestep_ratio
      #~ self.assertEquals(dto,540. | units.s)
      
      s.evolve_model(dto)

      self.assertEquals(s.model_time,dto)

      s1=QGCM(**default_options)
      
      s1.parameters.begin_time=s.model_time
      
      self.assertEqual(s1.get_name_of_current_state(), 'UNINITIALIZED')
      
      ch1=s.ocean_P_grid.new_channel_to(s1.ocean_P_grid)
      self.assertEqual(s1.get_name_of_current_state(), 'EDIT')

      ch1.copy_attributes(["pressure","dpressure_dt"])

      ch2=s.ocean_P_grid_forcings.new_channel_to(s1.ocean_P_grid_forcings)
      ch2.copy_attributes(["tau_x","tau_y"])
      ch3=s.ocean_T_grid_forcings.new_channel_to(s1.ocean_T_grid_forcings)
      ch3.copy_attributes(["surface_heat_flux"])
      # mixed_layer_depth

      ch4=s.ocean_T_grid.new_channel_to(s1.ocean_T_grid)
      ch4.copy_attributes(["surface_temperature_anomaly","dsurface_temperature_anomaly_dt"])


      self.assertEquals(s.ocean_P_grid[100,100].pressure,s1.ocean_P_grid[100,100].pressure)
      self.assertEquals(s.ocean_T_grid[100,100].surface_temperature_anomaly,s1.ocean_T_grid[100,100].surface_temperature_anomaly)
      self.assertEquals(s.ocean_T_grid[100,100].dsurface_temperature_anomaly_dt,s1.ocean_T_grid[100,100].dsurface_temperature_anomaly_dt)
      self.assertEquals(s.ocean_P_grid[100,100].dpressure_dt,s1.ocean_P_grid[100,100].dpressure_dt)
      self.assertEquals(s.ocean_P_grid_forcings[100,100].tau_x,s1.ocean_P_grid_forcings[100,100].tau_x)
      self.assertEquals(s.ocean_P_grid_forcings[100,100].tau_y,s1.ocean_P_grid_forcings[100,100].tau_y)
      self.assertEquals(s.ocean_T_grid_forcings[100,100].surface_heat_flux,s1.ocean_T_grid_forcings[100,100].surface_heat_flux)
      self.assertEquals(s.atmosphere_T_grid[100,50].mixed_layer_depth,s1.atmosphere_T_grid[100,50].mixed_layer_depth)
      self.assertEquals(s.atmosphere_T_grid[100,50].dmixed_layer_depth_dt,s1.atmosphere_T_grid[100,50].dmixed_layer_depth_dt)

      s1.evolve_model(2*dto)
      s.evolve_model(2*dto)
      
      self.assertEquals(s.model_time,s1.model_time)

      # diff should be max ~ 5.e-18, however due to error in dpo_dt ~ 1.e-16      
      d=s.ocean_P_grid.pressure-s1.ocean_P_grid.pressure
      print "abs diff in pressure:", abs(d).max()
      self.assertTrue(abs(d).max().number< 2.e-16)
      
    def test3(self):
      s=QGCM(**default_options)

      dto=0.5 | units.day
      
      s.evolve_model(dto)

      self.assertEquals(s.model_time,dto)

      s1=QGCM(**default_options)
      
      s1.parameters.begin_time=s.model_time
      
      self.assertEqual(s1.get_name_of_current_state(), 'UNINITIALIZED')
      
      ch1=s.ocean_P_grid.new_channel_to(s1.ocean_P_grid)
      self.assertEqual(s1.get_name_of_current_state(), 'EDIT')

      ch1.copy_attributes(["pressure","dpressure_dt"])

      ch2=s.ocean_P_grid_forcings.new_channel_to(s1.ocean_P_grid_forcings)
      ch2.copy_attributes(["tau_x","tau_y"])
      ch3=s.ocean_T_grid_forcings.new_channel_to(s1.ocean_T_grid_forcings)
      ch3.copy_attributes(["surface_heat_flux"])
      # mixed_layer_depth

      ch4=s.ocean_T_grid.new_channel_to(s1.ocean_T_grid)
      ch4.copy_attributes(["surface_temperature_anomaly","dsurface_temperature_anomaly_dt"])


      self.assertEquals(s.ocean_P_grid[100,100].pressure,s1.ocean_P_grid[100,100].pressure)
      self.assertEquals(s.ocean_T_grid[100,100].surface_temperature_anomaly,s1.ocean_T_grid[100,100].surface_temperature_anomaly)
      self.assertEquals(s.ocean_T_grid[100,100].dsurface_temperature_anomaly_dt,s1.ocean_T_grid[100,100].dsurface_temperature_anomaly_dt)
      self.assertEquals(s.ocean_P_grid[100,100].dpressure_dt,s1.ocean_P_grid[100,100].dpressure_dt)
      self.assertEquals(s.ocean_P_grid_forcings[100,100].tau_x,s1.ocean_P_grid_forcings[100,100].tau_x)
      self.assertEquals(s.ocean_P_grid_forcings[100,100].tau_y,s1.ocean_P_grid_forcings[100,100].tau_y)
      self.assertEquals(s.ocean_T_grid_forcings[100,100].surface_heat_flux,s1.ocean_T_grid_forcings[100,100].surface_heat_flux)
      self.assertEquals(s.atmosphere_T_grid[100,50].mixed_layer_depth,s1.atmosphere_T_grid[100,50].mixed_layer_depth)
      self.assertEquals(s.atmosphere_T_grid[100,50].dmixed_layer_depth_dt,s1.atmosphere_T_grid[100,50].dmixed_layer_depth_dt)

      s1.evolve_model(2*dto)
      s.evolve_model(2*dto)
      
      self.assertEquals(s.model_time,s1.model_time)

      d=s.ocean_P_grid.pressure-s1.ocean_P_grid.pressure
      print "abs diff in pressure:", abs(d).max()
      self.assertAlmostRelativeEquals(abs(d).max().number,1.19865228854e-12,11)
      # value from stand alone restart!

