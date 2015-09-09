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


    @nottest
    def test3(self):
        cdo = CDOInterface(**default_options)
        cdo.set_remap_file("weights/src_dst_con.nc")

        src_dims = cdo.get_src_grid_dims()['dims']
        dst_dims = cdo.get_dst_grid_dims()['dims']
        
        src = [[0 for x in range(src_dims['x'])] for y in range(src_dims['y'])]  # 180x90
        dst = [[0 for x in range(dst_dims['x'])] for y in range(dst_dims['y'])]  # 192x94

        for i in range (src_dims['x']):
            for j in range (src_dims['y']):
                src[i][j] = (i/10+j/10) % 2

        #still writing this test

        instance.stop()



