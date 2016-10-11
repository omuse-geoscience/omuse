from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.dales.interface import DalesInterface
from omuse.community.dales.interface import Dales

from amuse.units import units

default_options={}
#default_options=dict(number_of_workers=4)
#default_options=dict(channel_type="sockets",redirection="none")
#default_options=dict(redirection="none",debugger="gdb")

import logging
import numpy
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from nose.tools import nottest

class TestDalesInterface(TestWithMPI):

    def test1(self):

        print "Test 1: instantiate and clean up"

        instance = Dales(**default_options)
        instance.cleanup_code()
        instance.stop()

    def test2(self):

        print "Test 2: instantiate with 4 workers and clean up"

        instance = Dales(number_of_workers=4)
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):

        print "Test 3: instantiate, run one minute and clean up"

        instance = Dales(redirection="none")
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (1 | units.minute))
        newtim=instance.get_model_time()
        self.assertTrue(newtim>=tim)
        instance.cleanup_code()
        instance.stop()

    def test4(self):

        print "Test 4: instantiate with 4 workers, run one minute and clean up"

        instance = Dales(number_of_workers=4,redirection="none")
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (1 | units.minute))
        newtim=instance.get_model_time()
        self.assertTrue(newtim>=tim)
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):

        print "Test 5: instantiate, run 10 seconds and retrieve temperature profile"

        instance = Dales(number_of_workers=4,redirection="none")
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (10 | units.s))
        profile=instance.get_profile_field(numpy.arange(1,96))
        print "The retrieved profile is:",profile
        instance.cleanup_code()
        instance.stop()

    def test6(self):
        print "Test 6: instantiate, run 2 seconds and retrieve temperature profile using higher-level interface"

        instance = Dales(number_of_workers=4,redirection="none")
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (2 | units.s))

        temp = numpy.zeros(100)
        
        print "The retrieved U profile is:", instance.get_profile_U(temp)
        print "The retrieved V profile is:", instance.get_profile_V(temp)
        print "The retrieved W profile is:", instance.get_profile_W(temp)
        print 

        print "The retrieved THL profile is:", instance.get_profile_THL(temp)
        print "The retrieved QT profile is:", instance.get_profile_QT(temp)
        

        instance.cleanup_code()
        instance.stop()

            
