import os.path
import numpy
from amuse.test.amusetest import TestWithMPI
from time import sleep

from omuse.community.adcirc.interface import AdcircInterface,Adcirc

from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

default_options={}
#default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from omuse.community.adcirc.read_grid import adcirc_grid_reader, adcirc_parameter_reader
from omuse.community.adcirc.write_grid import adcirc_grid_writer, adcirc_parameter_writer

class TestAdcircInterface(TestWithMPI):

    def test1(self):
        print "Test 1: start"

        instance = AdcircInterface(**default_options)
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        print "Test 2: rootdir"

        instance = AdcircInterface(**default_options)
        instance.set_rootdir("data/test/2d")
        rootdir,err=instance.get_rootdir()
        self.assertEqual(rootdir,"data/test/2d")
        instance.cleanup_code()
        instance.stop()

    def test3(self):
        print "Test 3: depth"

        instance = AdcircInterface(**default_options)
        instance.set_rootdir("data/test/2d")
        instance.commit_parameters()
        depths,err=instance.get_node_depth([1,11])
        self.assertEqual(depths,[3.0480,9.3345])
        instance.cleanup_code()
        instance.stop()

    def test4(self):
        print "Test 4: sigma"

        instance = AdcircInterface(**default_options)
        instance.set_rootdir("data/test/3d")
        instance.commit_parameters()
        zn,err=instance.get_number_of_vertical_nodes()
        self.assertEqual(err,0)
        self.assertEqual(zn,21)
        
        s,z,err=instance.get_node_sigma(1,1)
        self.assertEqual(s,-1.)
        s,z,err=instance.get_node_sigma(1,zn)
        self.assertEqual(s,1.)
                
        instance.cleanup_code()
        instance.stop()


class TestAdcirc(TestWithMPI):

    def test1(self):
        print "Test 1: rootdir"

        instance = Adcirc(**default_options)
        instance.set_rootdir("data/test/2d")
        rootdir=instance.get_rootdir()
        self.assertEqual(rootdir,"data/test/2d")
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        print "Test 1: parameters"

        instance = Adcirc(**default_options)
        instance.initialize_code()
        instance.set_rootdir("data/test/2d")
        instance.commit_parameters()

        self.assertEqual(instance.parameters.use_interface_elevation_boundary, False)
        self.assertEqual(instance.parameters.use_interface_met_forcing, False)

        instance.parameters.use_interface_elevation_boundary=True
        instance.parameters.use_interface_met_forcing=True

        self.assertEqual(instance.parameters.use_interface_elevation_boundary, True)
        self.assertEqual(instance.parameters.use_interface_met_forcing, True)

        self.assertEqual(instance.parameters.bottom_friction_law, "linear")
        self.assertEqual(instance.parameters.linear_bottom_friction_coeff, 0. | units.s**-1)
        self.assertEqual(instance.parameters.quadratic_bottom_friction_coeff, 0. )
        self.assertEqual(instance.parameters.GWCE_weighting_factor, 0.005 )        

        instance.cleanup_code()
        instance.stop()

    def test3(self):
        instance = Adcirc(**default_options)
        instance.initialize_code()
        instance.set_rootdir("data/test/2d")
        instance.commit_parameters()

        self.assertEqual(len(instance.nodes),63)
        self.assertEqual(len(instance.elements),96)
        self.assertEqual(len(instance.forcings),63)
        
        self.assertEquals(instance.forcings.coriolis_f, (63*[0.0])| units.s**-1)
        self.assertEquals(instance.forcings.tau_x, ([0.0])| units.Pa)
        self.assertEquals(instance.forcings.tau_y, ([0.0])| units.Pa)

        instance.forcings.coriolis_f=numpy.arange(63) | units.s**-1
        self.assertEquals(instance.forcings.coriolis_f, range(63)| units.s**-1)
        forcings=instance.forcings.empty_copy()
        forcings.tau_x=(numpy.arange(63)+123) | units.Pa
        forcings.tau_y=numpy.arange(63) | units.Pa
        forcings.new_channel_to(instance.forcings).copy_attributes(["tau_x","tau_y"])
        self.assertEquals(instance.forcings.tau_x, range(123,123+63)| units.Pa)
        self.assertEquals(instance.forcings.tau_y, range(63)| units.Pa)

        instance.cleanup_code()
        instance.stop()


    def test4(self):
        instance = Adcirc(mode="3D",**default_options)
        instance.initialize_code()
        instance.set_rootdir("data/test/3d")
        instance.commit_parameters()
        
        a,b=1.,-1.
        sigma=numpy.arange(21)*0.1-1.
        z=(sigma-a)/(a-b)*(12.1920| units.m)-instance.nodes[11].eta
        self.assertEqual(instance.nodes[11].sigma,sigma)
        self.assertEqual(instance.nodes[11].z,z)
        instance.cleanup_code()
        instance.stop()
        
    def test5(self):
        instance = Adcirc(**default_options)
        instance.initialize_code()
        instance.set_rootdir("data/test/2d")
        instance.commit_parameters()
        depth=instance.nodes.depth[[0,10]]
        self.assertEqual(depth,[3.0480,9.3345] | units.m)
        instance.cleanup_code()
        instance.stop()
        
    def test6(self):
        instance = Adcirc(mode="3D",**default_options)
        instance.initialize_code()
        instance.set_rootdir("data/test/3d")
        instance.commit_parameters()
        vx=instance.nodes[11].wx
        vy=instance.nodes[11].wy
        vz=instance.nodes[11].wz
        self.assertEqual(vx.number,[0.]*21)
        instance.cleanup_code()
        instance.stop()

    def test7(self):
        instance=Adcirc(**default_options)
        instance.parameters.rootdir="data/test/2d"        

        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), "EDIT")
        
    def test8(self):
        instance=Adcirc(**default_options)
        instance.parameters.rootdir="data/test/2d"        

        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), "EDIT")
        self.assertRaises(Exception, lambda : setattr(instance.parameters, "use_interface_grid",True),
         expected_message="While calling before_set_interface_parameter of Adcirc: No transition from current state state 'EDIT' to state 'INITIALIZED' possible")
        
    def test9(self):
        instance=Adcirc(**default_options)
        instance.parameters.rootdir="data/test/2d"        

        instance.nodes.vx
        self.assertEqual(instance.get_name_of_current_state(), "EDIT")
        
    def test10(self):
        instance = Adcirc(**default_options)
        instance.set_rootdir("data/test/2d")

        print instance.parameters
        instance.evolve_model( 1 | units.hour)

        self.assertEqual(instance.model_time, 21*174.656| units.s)
        
        self.assertEquals(instance.nodes.status, u'wet')

        self.assertEquals(instance.elements.status, u'wet')
        
        instance.stop()

    def test11(self):
        instance = Adcirc(**default_options)
        instance.set_rootdir("data/test/2d")

        instance.evolve_model( 1 | units.hour)
        self.assertEqual(instance.model_time, 21*174.656| units.s)
        self.assertEquals(instance.nodes.status, u'wet')
        self.assertEquals(instance.elements.status, u'wet')
        
        instance.nodes[12].eta=0.21 | units.m
        self.assertEqual(instance.nodes[12].eta,0.21 | units.m)

        instance.nodes[12].status=u'dry'
        self.assertEqual(instance.nodes[12].status,u'dry')
        self.assertEqual(instance.nodes[12].eta,-instance.nodes[12].depth)

        instance.nodes[11].deta_dt=0.1 | units.m/units.s
        self.assertEqual(instance.nodes[11].deta_dt,0.1 | units.m/units.s)

        instance.nodes[11].eta=0.21 | units.m
        self.assertEqual(instance.nodes[11].eta,0.21 | units.m)


        instance.nodes[21].vx=0.21 | units.m/units.s
        self.assertEqual(instance.nodes[21].vx,0.21 | units.m/units.s)
        instance.nodes[21].vy=0.21 | units.m/units.s
        self.assertEqual(instance.nodes[21].vy,0.21 | units.m/units.s)


        
        
        instance.stop()

    def test12(self):
        instance = Adcirc(**default_options)
        instance.set_rootdir("data/test/2d")
        instance.commit_parameters()
 
        ref=instance.parameters.atmospheric_reference_pressure

        self.assertEqual(len(instance.nodes),63)
        self.assertEqual(len(instance.elements),96)
        self.assertEqual(len(instance.forcings),63)
        
        self.assertEquals(instance.forcings.coriolis_f, (63*[0.0])| units.s**-1)
        self.assertEquals(instance.forcings.tau_x, ([0.0])| units.Pa)
        self.assertEquals(instance.forcings.tau_y, ([0.0])| units.Pa)
        
        self.assertEquals(instance.forcings.pressure, ref)

        instance.forcings.coriolis_f=numpy.arange(63) | units.s**-1
        self.assertEquals(instance.forcings.coriolis_f, range(63)| units.s**-1)
        forcings=instance.forcings.empty_copy()
        forcings.tau_x=(numpy.arange(63)+123) | units.Pa
        forcings.tau_y=numpy.arange(63) | units.Pa
        forcings.pressure=(numpy.arange(63)+321.) | units.Pa
        forcings.new_channel_to(instance.forcings).copy_attributes(["pressure","tau_x","tau_y"])
        self.assertAlmostEquals(instance.forcings.pressure, range(321,321+63)| units.Pa,10)
        self.assertEquals(instance.forcings.tau_x, range(123,123+63)| units.Pa)
        self.assertEquals(instance.forcings.tau_y, range(63)| units.Pa)

        instance.stop()


AMIG=0.000140525700000 | units.s**-1
PER=2*numpy.pi/AMIG
EMO=0.3048 | units.m
DRAMP=2 | units.day

def ramp(t):
    return numpy.tanh(2*t/DRAMP)

def tidal_force_function(t):
    ncyc=numpy.floor(t/PER)
    aj=AMIG*(t-ncyc*PER)
    return ramp(t)*EMO*numpy.cos(aj)

class TestAdcircLong(TestWithMPI):

    def test1(self):
        tend=5*86400. | units.s

        param=adcirc_parameter_reader("data/test/2d/fort.15")
        param.read_parameters(NETA=9)
        param.parameters['NBFR']=-1

        gr=adcirc_grid_reader("data/test/2d/fort.14")
        gr.read_grid()
        nodes,elements,elev_boundary,flow_boundary=gr.get_sets()

        code=Adcirc()

        code._parameters=param.parameters        
        code.assign_grid_and_boundary(nodes,elements,elev_boundary, flow_boundary)

        code.parameters.use_interface_elevation_boundary=True
        code.parameters.use_interface_parameters=True
        code.parameters.use_interface_grid=True
        code.parameters.A_H=param.parameters["ESLM"] | units.m**2/units.s
        code.parameters.timestep=abs(param.parameters["DTDP"]) | units.s
        code.parameters.bottom_friction_law=["linear","quadratic","hybrid"][param.parameters["NOLIBF"]]
        try:
          code.parameters.linear_bottom_friction_coeff=param.parameters["TAU"]| units.s**-1
        except:
          pass
        try:
          code.parameters.quadratic_bottom_friction_coeff=param.parameters["CF"]
        except:
          pass
        code.parameters.use_predictor_corrector=param.parameters["DTDP"]<0
        code.parameters.use_interface_met_forcing=False

        print code.parameters

        tnow=code.model_time
        dt=code.parameters.timestep
  
        elev_boundaries= list(code.elevation_boundaries())
  
        eta61=[]
        time=[]
        forcing=[]
  
        while tnow<tend-dt/2:
            elev_boundaries[0].eta=tidal_force_function(tnow+dt/2)
            code.evolve_model(tnow+dt)
            tnow=code.get_model_time()
  
            eta=code.nodes[60].eta.number
            time.append(tnow.number)
            eta61.append(eta)  
            forcing.append(elev_boundaries[0].eta[0].number)
      
        code.stop()
        
        from matplotlib import pyplot
          
        pyplot.ion()
        f=pyplot.figure(figsize=(8,6))
        pyplot.show()

        pyplot.clf()
        pyplot.plot(time,eta61,'r+')
        pyplot.plot(time,forcing,'g+')
        pyplot.plot(time,tidal_force_function((time| units.s)).number)
        pyplot.draw()
        sleep(3)

    def test2(self):
        tend=5*86400. | units.s

        param=adcirc_parameter_reader("data/test/2d/fort.15")
        param.read_parameters(NETA=9)

        gr=adcirc_grid_reader("data/test/2d/fort.14")
        gr.read_grid()
        nodes,elements,elev_boundary,flow_boundary=gr.get_sets()

        code=Adcirc()

        code._parameters=param.parameters        
        code.assign_grid_and_boundary(nodes,elements,elev_boundary, flow_boundary)

        code.parameters.use_interface_elevation_boundary=False
        code.parameters.use_interface_parameters=True
        code.parameters.use_interface_grid=True
        code.parameters.A_H=param.parameters["ESLM"] | units.m**2/units.s
        code.parameters.timestep=abs(param.parameters["DTDP"]) | units.s
        code.parameters.bottom_friction_law=["linear","quadratic","hybrid"][param.parameters["NOLIBF"]]
        try:
          code.parameters.linear_bottom_friction_coeff=param.parameters["TAU"]| units.s**-1
        except:
          pass
        try:
          code.parameters.quadratic_bottom_friction_coeff=param.parameters["CF"]
        except:
          pass
        code.parameters.use_predictor_corrector=param.parameters["DTDP"]<0
        code.parameters.use_interface_met_forcing=False

        print code.parameters

        tnow=code.model_time
        dt=code.parameters.timestep
  
        elev_boundaries= list(code.elevation_boundaries())
  
        eta61=[]
        time=[]
        forcing=[]
  
        while tnow<tend-dt/2:
            code.evolve_model(tnow+dt)
            tnow=code.get_model_time()
  
            eta=code.nodes[60].eta.number
            time.append(tnow.number)
            eta61.append(eta)  
            forcing.append(elev_boundaries[0].eta[0].number)
      
        code.stop()
        
        from matplotlib import pyplot
          
        pyplot.ion()
        f=pyplot.figure(figsize=(8,6))
        pyplot.show()

        pyplot.clf()
        pyplot.plot(time,eta61,'r+')
        pyplot.plot(time,forcing,'g+')
        pyplot.plot(time,tidal_force_function((time| units.s)).number)
        pyplot.draw()
        sleep(3)
