from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.cdo.interface import CDORemapper

from amuse.units import units

from nose.tools import nottest

default_options=dict(redirection="none", channel="sockets")

class CDORemapperTests(TestWithMPI):
    
    def test1(self):
        r = CDORemapper(**default_options)
        print r
        r.stop()


    #test initialization through weights file and check state
    def test2(self):
        r = CDORemapper(**default_options)
        self.assertEquals(r.state_machine._current_state.name, 'UNINITIALIZED')

        filename = "weights/src_dst_con.nc"
        r.parameters.weights_file = filename
        self.assertEquals(r.state_machine._current_state.name, 'INITIALIZED')

        get_file = r.parameters.weights_file
        self.assertEquals(get_file, filename)

        src_grid_size = r.get_src_grid_size()
        self.assertEquals(src_grid_size, 16200)
        self.assertEquals(r.state_machine._current_state.name, 'RUN')

        r.stop()



