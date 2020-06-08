from math import isnan
from amuse.test.amusetest import TestWithMPI

from omuse.community.iemic.interface import iemicInterface
from omuse.community.iemic.interface import iemic

class iemicInterfaceTests(TestWithMPI):

    def test1(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        err = instance.commit_parameters()
        self.assertEqual(err,0)

        err = instance.set_double_parameter("continuation", "destination 0", 1.0)
        self.assertEqual(err,0)

        err = instance.commit_continuation_parameters()
        self.assertEqual(err,0)

        err = instance.initialize_code()
        self.assertEqual(err,0)

        err = instance.step_continuation()
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

        err = instance.set_double_parameter("continuation", "destination 0", 1.0)
        self.assertEqual(err,0)

        err = instance.commit_continuation_parameters()
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
        print("!!>> |{0}|, |{1}|, |{2}|".format(paramSet, sublist,sublist+param)) 
        paramType, err = instance.get_parameter_type(paramSet, sublist + param)
        self.assertEqual(err,0)

        if paramType == "unknown":
            assert 0, "type can't be unknown"

        print(("    " * n) + param + ": " + paramType, end=' ')

        param = sublist + param

        if paramType == "ParameterList":
            subParamCount, err = instance.get_num_parameters(paramSet, param)
            self.assertEqual(err,0)

            print()

            for i in range(0, subParamCount):
                subParamName, err = instance.get_parameter_name(paramSet, param, i)
                self.assertEqual(err,0)
                self.paramHelper(instance, paramSet, n+1, param + "?", subParamName)

        elif paramType == "bool":
            val, err = instance.get_bool_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance.get_default_bool_parameter(paramSet, param)
            self.assertEqual(err,0)
            self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance.set_bool_parameter(paramSet, param, val)
            self.assertEqual(err,0)

        elif paramType == "char":
            print()

        elif paramType == "double":
            val, err = instance.get_double_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance.get_default_double_parameter(paramSet, param)
            self.assertEqual(err,0)
            if isnan(val):
                self.assertEqual(isnan(val), isnan(default))
            else:
                self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance.set_double_parameter(paramSet, param, val)
            self.assertEqual(err,0)

        elif paramType == "int":
            val, err = instance.get_int_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance.get_default_int_parameter(paramSet, param)
            self.assertEqual(err,0)
            self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance.set_int_parameter(paramSet, param, val)
            self.assertEqual(err,0)

        elif paramType == "string":
            val, err = instance.get_string_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance.get_default_string_parameter(paramSet, param)
            self.assertEqual(err,0)
            self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance.set_string_parameter(paramSet, param, val)

        else:
            assert 0, "unknow type!"

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

            print(paramSet + "(" + str(paramCount) + "):")

            for j in range(0, paramCount):
                paramName, err = instance.get_parameter_name(paramSet, "", j)
                self.assertEqual(err,0)

                self.paramHelper(instance, paramSet, 0, "", paramName)

        err = instance.cleanup_code()
        self.assertEqual(err,0)

        instance.stop()

class iemicTests(TestWithMPI):
    def test1(self):
        instance = iemic()

        sets = instance.parameter_set_names()

        for name in sets:
          print("parameter set: {0}".format(name))
          print(getattr(instance,name))
          print()

        instance.cleanup_code()
        instance.stop()
