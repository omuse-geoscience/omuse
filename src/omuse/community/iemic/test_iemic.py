from amuse.test.amusetest import TestWithMPI

from .interface import iemicInterface
from .interface import iemic

class iemicInterfaceTests(TestWithMPI):

    def test1(self):
        instance = iemicInterface(redirection="none")
        err = instance.initialize()
        assert not err, "err should be 0"
        err = instance.commit_parameters()
        assert not err, "err should be 0"
        err = instance.initialize_code()
        assert not err, "err should be 0"
        err = instance.step()
        assert not err, "err should be 0"
        err = instance.cleanup_code()
        assert not err, "err should be 0"
        err = instance.stop()
        assert not err, "err should be 0"

    def test2(self):
        instance = iemicInterface(redirection="none")
        err = instance.initialize()
        assert not err, "err should be 0"
        err = instance.commit_parameters()
        assert not err, "err should be 0"
        err = instance.initialize_code()
        assert not err, "err should be 0"
        val, err = instance.get_u([0],[0],[0])
        assert not err, "err should be 0"
        val, err = instance.get_v([0],[0],[0])
        assert not err, "err should be 0"
        val, err = instance.get_w([0],[0],[0])
        assert not err, "err should be 0"
        val, err = instance.get_p([0],[0],[0])
        assert not err, "err should be 0"
        val, err = instance.get_t([0],[0],[0])
        assert not err, "err should be 0"
        val, err = instance.get_s([0],[0],[0])
        assert not err, "err should be 0"
        err = instance.cleanup_code()
        assert not err, "err should be 0"
        err = instance.stop()
        assert not err, "err should be 0"
