from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from omuse.community.oifs.interface import OpenIFSInterface
from omuse.community.oifs.interface import OpenIFS

import logging
import nose.tools

logging.basicConfig(level = logging.DEBUG)
log = logging.getLogger(__name__)

class TestOpenIFSInterface(TestWithMPI):

    def test_single_proc(self):
        log.info("Testing instantiate OpenIFS with single MPI rank and cleaning up")
        instance = OpenIFS()
        try:
            instance.initialize_code()
            instance.commit_grid()
        finally:
            instance.cleanup_code()
            instance.stop()

    def test_4proc(self):
        log.info("Testing instantiate OpenIFS with 4 MPI ranks and cleaning up")
        instance = OpenIFS(number_of_workers = 4)
        try:
            instance.initialize_code()
            instance.commit_grid()
        finally:
            instance.cleanup_code()
            instance.stop()

    def test_evolve_model(self):
        log.info("Testing running OpenIFS for one hour")
        instance = OpenIFS()
        try:
            instance.initialize_code()
            tim = instance.get_model_time()
            instance.commit_grid()
            instance.evolve_model(tim + (1 | units.hour))
            newtim = instance.get_model_time()
            self.assertTrue(newtim >= tim)
        finally:
            instance.cleanup_code()
            instance.stop()

    def test_evolve_model_4proc(self):
        log.info("Testing running OpenIFS for eight hours with 4 MPI ranks")
        instance = OpenIFS(number_of_workers = 4)
        try:
            instance.initialize_code()
            tim = instance.get_model_time()
            instance.commit_grid()
            instance.evolve_model(tim + (8 | units.hour))
            newtim = instance.get_model_time()
            self.assertTrue(newtim >= tim)
        finally:
            instance.cleanup_code()
            instance.stop()
