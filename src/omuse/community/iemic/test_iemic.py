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

    def paramHelper(self, instance, paramSet, n, sublist, param):
        paramType, err = instance.get_parameter_type(paramSet, sublist + param)
        self.assertEqual(err,0)

        if paramType == "unknown":
            assert 0, "type can't be unknown"

        print ("    " * n) + param + ": " + paramType

        param = sublist + param

        if paramType == "ParameterList":
            subParamCount, err = instance.get_num_parameters(paramSet, param)
            self.assertEqual(err,0)

            for i in range(0, subParamCount):
                subParamName, err = instance.get_parameter_name(paramSet, param, i)
                self.assertEqual(err,0)
                self.paramHelper(instance, paramSet, n+1, param + "?", subParamName)

    def test3(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        paramSetCount, err = instance.get_num_parameter_sets()
        self.assertEqual(err,0)

        for i in range(0, paramSetCount):
            paramSet, err = instance.get_parameter_set_name(i)
            self.assertEqual(err,0)

            paramCount, err = instance.get_num_parameters(paramSet)
            self.assertEqual(err,0)

            print paramSet + "(" + str(paramCount) + "):"

            for j in range(0, paramCount):
                paramName, err = instance.get_parameter_name(paramSet, "", j)
                self.assertEqual(err,0)

                self.paramHelper(instance, paramSet, 0, "", paramName)

        err = instance.cleanup_code()
        self.assertEqual(err,0)

        instance.stop()
