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
        instance.get_u([0],[0],[0])
        instance.get_v([0],[0],[0])
        instance.get_w([0],[0],[0])
        instance.get_p([0],[0],[0])
        instance.get_t([0],[0],[0])
        instance.get_s([0],[0],[0])
        instance.cleanup_code()
        instance.stop()
