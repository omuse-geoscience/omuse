from amuse.test.amusetest import TestWithMPI

from .interface import iemicInterface
from .interface import iemic

class iemicInterfaceTests(TestWithMPI):

    def test1(self):
        instance = iemicInterface()
        instance.initialize()
        instance.commit_parameters()
        instance.initialize_code()
        instance.stop()

    def test2(self):
        instance = iemicInterface()
        instance.initialize()
        instance.commit_parameters()
        instance.initialize_code()
        instance.step()
        instance.stop()
