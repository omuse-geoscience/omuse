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
        instance.set_rootdir("data/test/")
        rootdir,err=instance.get_rootdir()
        self.assertEqual(rootdir,"data/test/")
        instance.cleanup_code()
        instance.stop()

class TestAdcirc(TestWithMPI):

    def test1(self):
        print "Test 1: rootdir"

        instance = Adcirc(**default_options)
        instance.set_rootdir("data/test/")
        rootdir=instance.get_rootdir()
        self.assertEqual(rootdir,"data/test/")
        instance.cleanup_code()
        instance.stop()
