import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

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
        
        print instance.nodes.eta,instance.nodes.eta_prev

        instance.stop()
