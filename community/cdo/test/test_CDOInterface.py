from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.cdo.interface import CDOInterface

from amuse.units import units

from nose.tools import nottest

default_options=dict(redirection="none", channel="sockets")

class CDOInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = CDOInterface(**default_options)
        print instance
        instance.stop()


    #test initialization through weights file
    def test2(self):
        instance = CDOInterface(**default_options)
        instance.initialize_code()
        filename = "weights/src_dst_con.nc"
        instance.set_weights_file(filename)
        instance.commit_parameters()
        get_file = instance.get_weights_file()['filename']

        self.assertEquals(get_file, filename)

        self.assertEquals(instance.get_src_grid_size()['size'], 16200)
        instance.stop()


    #test initialization through grid files
    def test3(self):
        instance = CDOInterface(**default_options)
        instance.initialize_code()
        src_grid_file = "grids/src_grid.nc"
        dst_grid_file = "grids/dst_grid.nc"
        instance.set_src_grid_file(src_grid_file)
        instance.set_dst_grid_file(dst_grid_file)
        instance.commit_parameters()
        get_file = instance.get_src_grid_file()['filename']
        self.assertEquals(get_file, src_grid_file)
        get_file = instance.get_dst_grid_file()['filename']
        self.assertEquals(get_file, dst_grid_file)

        self.assertEquals(instance.get_src_grid_size()['size'], 16200)
        self.assertEquals(instance.get_dst_grid_size()['size'], 18048)
        instance.stop()


    @nottest
    def test4(self):
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



