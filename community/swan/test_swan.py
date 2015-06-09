import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from omuse.community.swan.interface import SwanInterface

from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

#default_options={}
default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

class TestSwanInterface(TestWithMPI):
    def test1(self):
      s=SwanInterface(**default_options)
      s.initialize_code()
      s.cleanup_code()
