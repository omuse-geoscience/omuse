from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.dales.interface import DalesInterface
from omuse.community.dales.interface import Dales

from amuse.units import units

#default_options=dict(channel_type="sockets",redirection="none")
default_options=dict(redirection="none")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from nose.tools import nottest

class DalesTests(TestWithMPI):

    #test the behavior of the state machine
    def test1(self):
        instance = Dales(**default_options)
        print "executing test..."

        self.assertEquals(instance.state_machine._current_state.name, 'UNINITIALIZED')

        #proceed to evolve
        time = instance.get_model_time()
        tend = time + instance.get_timestep_next()
        instance.evolve_model(tend)

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'EVOLVED')

        instance.stop()

        self.assertEquals(instance.state_machine._current_state.name, 'UNINITIALIZED')
