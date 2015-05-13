import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from omuse.community.adcirc.interface import AdcircInterface,Adcirc

from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

#default_options={}
default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

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
        print "Test 2: rootdir"

        instance = AdcircInterface(**default_options)
        instance.set_rootdir("data/test/3d")
        instance.commit_parameters()
        zn,err=instance.get_number_of_vertical_nodes()
        self.assertEqual(err,0)
        self.assertEqual(zn,21)
        
        s,err=instance.get_node_sigma(1,1)
        self.assertEqual(s,-1.)
        s,err=instance.get_node_sigma(zn,1)
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
        
        print instance.nodes.sigma
        print instance.nodes[11].sigma
        instance.cleanup_code()
        instance.stop()
        
        
