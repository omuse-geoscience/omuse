from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.cdo.interface import CDOInterface

from amuse.units import units

default_options=dict(redirection="none", channel="sockets")

class CDOInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = CDOInterface(**default_options)
        print instance
        instance.stop()

    def test2(self):
        instance = CDOInterface(**default_options)
        instance.set_remap_file("weights/src_dst_con.nc")

        self.assertEquals(instance.get_src_grid_size()['size'], 16200)
        instance.stop()

