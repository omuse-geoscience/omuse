from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.dales.interface import DalesInterface
from omuse.community.dales.interface import Dales

from amuse.units import units

default_options=dict(redirection="none")
#default_options=dict(number_of_workers=4)
#default_options=dict(channel_type="sockets",redirection="none")
#default_options=dict(redirection="none",debugger="gdb")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from nose.tools import nottest

class TestDalesInterface(TestWithMPI):

    def test1(self):

        print "Test 1: start"

        instance = Dales(**default_options)
        instance.cleanup_code()

        instance.stop()

    def test2(self):

        print "Test 2: start"

        instance = Dales(redirection="none",number_of_workers=4)
        instance.cleanup_code()

        instance.stop()
    
#    def test2(self):
#
#        print "Test 2: default input"
#
#        instance = Dales(**default_options)
#        inputfile,err=instance.get_input_file()
#
#        self.assertEquals(inputfile,"namoptions.001")
#        instance.cleanup_code()
#        instance.stop()
