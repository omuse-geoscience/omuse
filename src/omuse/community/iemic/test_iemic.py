from amuse.test.amusetest import TestWithMPI

from .interface import iemicInterface
from .interface import iemic

class iemicInterfaceTests(TestWithMPI):

    def test1(self):
        instance = iemicInterface()
        instance.initialize()
        instance.commit_parameters()
        instance.initialize_code()
        instance.step()
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        instance = iemicInterface()
        instance.initialize()
        instance.commit_parameters()
        instance.initialize_code()
        instance.get_u()
        instance.get_v()
        instance.get_w()
        instance.get_p()
        instance.get_t()
        instance.get_s()
        instance.cleanup_code()
        instance.stop()
