from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.pop.interface import POPInterface
from omuse.community.pop.interface import POP

from amuse.units import units

default_options=dict(channel_type="sockets",redirection="none")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from nose.tools import nottest



class POPTests(TestWithMPI):


    #test the behavior of the state machine
    #@nottest
    def test1(self):
        instance = POP(**default_options)

        self.assertEquals(instance.state_machine._current_state.name, 'UNINITIALIZED')

        #a read of a parameter requires the state to be either EDIT or RUN, which means we pass through INITIALIZED
        fcor = instance.forcings.coriolis_f 

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'RUN')

        #check if we read a sensible value
        print 'fcor[1,1] = ', fcor[1,1]
        self.assertTrue(fcor[1,1] != 0 | units.s**-1, msg='Expected coriolis force to be not equal to zero for node 1,1')

        #proceed to evolve
        time = instance.get_model_time()
        tend = time + instance.get_timestep_next()
        instance.evolve_model(tend)

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'EVOLVED')

        #try to read something again, should cause a state transition to RUN
        fcor = instance.forcings.coriolis_f 

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'RUN')

        #check if we can write to coriolis_f
        instance.forcings.coriolis_f = fcor

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'EDIT')
        
        instance.stop()



    #@nottest
    def test2(self):
        instance = POP(**default_options)

        #proceed to evolve
        #time = instance.get_model_time()
        #tend = time + (1 | units.day)
        #instance.evolve_model(tend)

        #extract the tau_x and tau_y from POP
        tau_x = instance.forcings.tau_x
        tau_y = instance.forcings.tau_y

        print 'tau_x='
        print tau_x[range(1,5),1]
#        print tau_x
        print 'tau_y='
        print tau_y[range(1,5),1]
#        print tau_y

        #check if tau_x and tau_y are not only zeroes, this also fails for things like NaN and Inf
        self.assertTrue(numpy.sum(numpy.abs(tau_x)) > (0.0 | units.Pa), msg='Expected tau_x to contain some actual values')
        self.assertTrue(numpy.sum(numpy.abs(tau_y)) > (0.0 | units.Pa), msg='Expected tau_x to contain some actual values')

        #check to see if we can set and retrieve the same wind stress
        size = instance.get_domain_size()
        tau_x = numpy.random.random(size) | units.Pa
        tau_y = numpy.random.random(size) | units.Pa

        #cannot write to forcings directly unfortunately
        #instance.forcings.tau_x = tau_x
        #instance.forcings.tau_y = tau_y
        forcings=instance.forcings.empty_copy()
        forcings.tau_x=tau_x
        forcings.tau_y=tau_y
        forcings.new_channel_to(instance.forcings).copy_attributes(["tau_x","tau_y"])

        self.assertEquals(instance.state_machine._current_state.name, 'EDIT_FORCINGS')

        #now retrieve the wind stress from the model
        tau_x_returned = instance.forcings.tau_x
        tau_y_returned = instance.forcings.tau_y

        self.assertEquals(instance.state_machine._current_state.name, 'RUN')

        #almost equals
        self.assertAlmostRelativeEqual(tau_x, tau_x_returned, places=10)
        self.assertAlmostRelativeEqual(tau_y, tau_y_returned, places=10)

        #check to see if the forcings we have set are not overwritten by internal routines,
        #evolve the model 
        time = instance.get_model_time()
        tend = time + instance.get_timestep_next()
        instance.evolve_model(tend)

        #check if wind stress is still the same
        tau_x_returned = instance.forcings.tau_x
        tau_y_returned = instance.forcings.tau_y
        self.assertAlmostRelativeEqual(tau_x, tau_x_returned, places=10)
        self.assertAlmostRelativeEqual(tau_y, tau_y_returned, places=10)

        instance.stop()



    def test3(self):
        p = POP(**default_options)

        mystr = ''
        mynum = 0
        bogus_file = '/fake/path/to/file'

        mystr = p.get_ts_option()
        self.assertEquals(mystr, 'internal')
        p.set_ts_option('restart')
        mystr = p.get_ts_option()
        self.assertEquals(mystr, 'restart')
        mystr = p.get_ts_file()
        self.assertEquals(mystr, '')
        p.set_ts_file(bogus_file)
        mystr = p.get_ts_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_ts_file_format()
        self.assertEquals(mystr, 'bin')
        p.set_ts_file_format('nc')
        mystr = p.get_ts_file_format()
        self.assertEquals(mystr, 'nc')

        mystr = p.get_distribution()
        self.assertEquals(mystr, 'cartesian')

        p.set_distribution('predefined')
        mystr = p.get_distribution()
        self.assertEquals(mystr, 'predefined')

        mystr = p.get_distribution_file()
        self.assertEquals(mystr, '')
        p.set_distribution_file(bogus_file)
        mystr = p.get_distribution_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_ew_boundary_type()
        self.assertEquals(mystr, 'cyclic')
        p.set_ew_boundary_type('closed')
        mystr = p.get_ew_boundary_type()
        self.assertEquals(mystr, 'closed')
        mystr = p.get_ns_boundary_type()
        self.assertEquals(mystr, 'closed')
        p.set_ns_boundary_type('tripole')
        mystr = p.get_ns_boundary_type()
        self.assertEquals(mystr, 'tripole')

        mystr = p.get_restart_option()
        self.assertEquals(mystr, 'never')
        p.set_restart_option('nday')
        mystr = p.get_restart_option()
        self.assertEquals(mystr, 'nday')
        mynum = p.get_restart_freq_option()
        self.assertEquals(mynum, 1)
        p.set_restart_freq_option(5)
        mynum = p.get_restart_freq_option()
        self.assertEquals(mynum, 5)
        mystr = p.get_restart_file()
        self.assertEquals(mystr, '')
        p.set_restart_file(bogus_file)
        mystr = p.get_restart_file()
        self.assertEquals(mystr, bogus_file)



