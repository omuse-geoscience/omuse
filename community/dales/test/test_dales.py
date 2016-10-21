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
        print "Test 6: instantiate, run 2 seconds and retrieve profiles"

        instance = Dales(number_of_workers=4,redirection="none")
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (2 | units.s))

        
        print "The retrieved U profile is:", instance.get_profile_U()
        print "The retrieved V profile is:", instance.get_profile_V()
        print "The retrieved W profile is:", instance.get_profile_W()
        print 

        print "The retrieved THL profile is:", instance.get_profile_THL()
        print "The retrieved QT profile is:", instance.get_profile_QT()
        print "The retrieved zf levels:", instance.get_zf()
        print "The retrieved zh levels:", instance.get_zh()

        i,j,k,xsize,ysize = instance.get_params_grid()
        print "Grid size", (i, j, k)
        print "Horizontal extent of model", (xsize,ysize)

        instance.cleanup_code()
        instance.stop()

            
    def test7(self):
        print "Test 7: instantiate, set a forcing, run 2 seconds and retrieve profiles"

        instance = Dales(number_of_workers=4,redirection="none")
        tim=instance.get_model_time()
        instance.commit_grid()

        #print "Number of height levels:", instance.k
        #tend = numpy.arange(instance.k)*.01
        #print ('U_tend in python', tend)
        #instance.set_tendency_U(tend)

        tend = numpy.zeros(instance.k)
        tend[5] = -.1
        tend[6] = -.2   # layer 9 sees the most effect of this
        tend[7] = -.1

        instance.set_tendency_U(tend)
        
        instance.evolve_model(tim + (30 | units.s))

        
        print "The retrieved U profile is:", instance.get_profile_U()
        print "The retrieved V profile is:", instance.get_profile_V()
        #print "The retrieved W profile is:", instance.get_profile_W()
        print 

        #print "The retrieved THL profile is:", instance.get_profile_THL()
        print "The retrieved QT profile is:", instance.get_profile_QT()
        #print "The retrieved zf levels:", instance.get_zf()
        #print "The retrieved zh levels:", instance.get_zh()

        i,j,k,xsize,ysize = instance.get_params_grid()
        print "Grid size", (i, j, k)
        print "Horizontal extent of model", (xsize,ysize)

        instance.cleanup_code()
        instance.stop()
