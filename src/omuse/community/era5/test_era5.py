# todo:
# test all units are defined
# test all shortnames are defined
# test instantiation
# test basic functionality

import os
import sys
import numpy

from omuse.units import units

from amuse.test.amusetest import TestWithMPI

from omuse.community.era5 import era5
from omuse.community.era5.interface import ERA5, _era5_units_to_omuse

import datetime

class testera5(TestWithMPI):
    def test0(self):
        """ test single level units"""
        for field in era5.SLVARS:
            ustring=era5.UNITS[field]
            if ustring not in _era5_units_to_omuse:
                print("%s unit unknown"%ustring)

class TestERA5(TestWithMPI):
    
    def test0(self):
        instance=ERA5()
