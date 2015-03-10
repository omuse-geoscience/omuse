import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from omuse.community.adcirc.interface import AdcircInterface

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
        print "Test 2: test"

        instance = AdcircInterface(**default_options)
        instance.test_amuse()
        instance.cleanup_code()
        instance.stop()
