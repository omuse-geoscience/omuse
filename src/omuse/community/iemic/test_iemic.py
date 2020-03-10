from amuse.test.amusetest import TestWithMPI

from .interface import iemicInterface
from .interface import iemic

class iemicInterfaceTests(TestWithMPI):

    def test1(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        err = instance.commit_parameters()
        self.assertEqual(err,0)

        err = instance.initialize_code()
        self.assertEqual(err,0)

        err = instance.step()
        self.assertEqual(err,0)

        err = instance.cleanup_code()
        self.assertEqual(err,0)

        instance.stop()

    def test2(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        err = instance.commit_parameters()
        self.assertEqual(err,0)

        err = instance.initialize_code()
        self.assertEqual(err,0)

        val, err = instance.get_u([0],[0],[0])
        self.assertEqual(err,0)

        val, err = instance.get_v([0],[0],[0])
        self.assertEqual(err,0)

        val, err = instance.get_w([0],[0],[0])
        self.assertEqual(err,0)

        val, err = instance.get_p([0],[0],[0])
        self.assertEqual(err,0)

        val, err = instance.get_t([0],[0],[0])
        self.assertEqual(err,0)

        val, err = instance.get_s([0],[0],[0])
        self.assertEqual(err,0)

        err = instance.test_grid("")
        self.assertEqual(err,0)

        err = instance.cleanup_code()
        self.assertEqual(err,0)

        instance.stop()
