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


