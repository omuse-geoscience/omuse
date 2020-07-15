from math import isnan
from amuse.test.amusetest import TestWithMPI
from amuse.support.exceptions import AmuseException

from omuse.community.iemic.interface import iemicInterface
from omuse.community.iemic.interface import iemic

class iemicInterfaceTests(TestWithMPI):

    def test1(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        err = instance.commit_parameters()
        self.assertEqual(err,0)

        err = instance._set_double_parameter("continuation", "destination 0", 1.0)
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

        err = instance._set_double_parameter("continuation", "destination 0", 1.0)
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
        paramType, err = instance._get_parameter_type(paramSet, sublist + param)
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
                self.paramHelper(instance, paramSet, n+1, param + "->", subParamName)

        elif paramType == "bool":
            val, err = instance._get_bool_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance._get_default_bool_parameter(paramSet, param)
            self.assertEqual(err,0)
            self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance._set_bool_parameter(paramSet, param, val)
            self.assertEqual(err,0)

        elif paramType == "char":
            print()

        elif paramType == "double":
            val, err = instance._get_double_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance._get_default_double_parameter(paramSet, param)
            self.assertEqual(err,0)
            if isnan(val):
                self.assertEqual(isnan(val), isnan(default))
            else:
                self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance._set_double_parameter(paramSet, param, val)
            self.assertEqual(err,0)

        elif paramType == "int":
            val, err = instance._get_int_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance._get_default_int_parameter(paramSet, param)
            self.assertEqual(err,0)
            self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance._set_int_parameter(paramSet, param, val)
            self.assertEqual(err,0)

        elif paramType == "string":
            val, err = instance._get_string_parameter(paramSet, param)
            self.assertEqual(err,0)

            default, err = instance._get_default_string_parameter(paramSet, param)
            self.assertEqual(err,0)
            self.assertEqual(val, default)

            print(" (value: " + str(val) + ", default: " + str(default) + ")")

            err = instance._set_string_parameter(paramSet, param, val)

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

    def test4(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        err = instance.save_xml_parameters("ocean", "test.xml")
        self.assertEqual(err,0)

        err = instance._load_xml_parameters("ocean", "test.xml")
        self.assertEqual(err,0)

        err = instance.cleanup_code()
        self.assertEqual(err,0)

        instance.stop()

    def test5(self):
        instance = iemicInterface()
        err = instance.initialize()
        self.assertEqual(err,0)

        err = instance.commit_parameters()
        self.assertEqual(err,0)

        err = instance._set_double_parameter("ocean", "THCM->Global Bound xmin", 5.0)
        self.assertEqual(err,-1)

        err = instance.recommit_parameters()
        self.assertEqual(err,0)

        err = instance._set_double_parameter("ocean", "THCM->Starting Parameters->SPL1", 5.0)
        self.assertEqual(err,0)

        err = instance.recommit_parameters()
        self.assertEqual(err,0)

        instance.cleanup_code()
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

    def test2(self):
        "Test whether all ocean parameters cause the expected state transitions."

        instance = iemic()

        for param_set_name in instance.parameter_set_names():
            params = getattr(instance, param_set_name)
            for param in params.iter_parameters():
                param.set_value(param.get_value())

        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(),
                        "OCEAN-PARAM-CONTINUATION-NOPARAM")

        for param_set_name in instance.parameter_set_names():
            params = getattr(instance, param_set_name)
            isOcean = param_set_name.startswith("ocean")
            for param in params.iter_parameters():
                try:
                    param.set_value(param.get_value())
                    if isOcean:
                        self.assertEqual(instance.get_name_of_current_state(),
                                        "OCEAN-UPDATED-CONTINUATION-NOPARAM")

                        instance.recommit_parameters()
                        self.assertEqual(instance.get_name_of_current_state(),
                                        "OCEAN-PARAM-CONTINUATION-NOPARAM")
                except AmuseException as e:
                    self.assertTrue(isOcean)

        instance.cleanup_code()
        instance.stop()

    def test3(self):
        "Test whether all parameters cause the expected state transitions."
        instance = iemic()

        instance.commit_parameters()
        instance.set_parameter("continuation->destination 0", 1.0)
        instance.commit_continuation_parameters()
        self.assertEqual(instance.get_name_of_current_state(),
                        "OCEAN-PARAM-CONTINUATION-PARAM")

        for param_set_name in instance.parameter_set_names():
            params = getattr(instance, param_set_name)
            isOcean = param_set_name.startswith("ocean")
            for param in params.iter_parameters():
                try:
                    param.set_value(param.get_value())
                    if isOcean:
                        self.assertEqual(instance.get_name_of_current_state(),
                                        "OCEAN-UPDATED-CONTINUATION-PARAM")

                        instance.recommit_parameters()
                        self.assertEqual(instance.get_name_of_current_state(),
                                        "OCEAN-PARAM-CONTINUATION-PARAM")
                    else:
                        self.assertEqual(instance.get_name_of_current_state(),
                                        "OCEAN-PARAM-CONTINUATION-UPDATED")

                        instance.recommit_continuation_parameters()
                        self.assertEqual(instance.get_name_of_current_state(),
                                        "OCEAN-PARAM-CONTINUATION-PARAM")
                except AmuseException as e:
                    pass

        instance.cleanup_code()
        instance.stop()

    def test4(self):
        "Test whether ocean methods trigger state transitions properly."
        simple_methods = {
            "get_u" : (0, 0, 0),
            "get_v" : (0, 0, 0),
            "get_w" : (0, 0, 0),
            "get_p" : (0, 0, 0),
            "get_t" : (0, 0, 0),
            "get_s" : (0, 0, 0),
            "get_nrange" : (),
            "get_mrange" : (),
            "get_lrange" : ()
        }

        for method_name, args in simple_methods.items():
            instance = iemic()
            self.assertEqual(instance.get_name_of_current_state(),
                             "UNINITIALIZED")

            getattr(instance, method_name)(*args)
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-PARAM-CONTINUATION-NOPARAM")

            instance.set_parameter("ocean->THCM->Starting Parameters->SPL1", 2000.0)
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-UPDATED-CONTINUATION-NOPARAM")

            getattr(instance, method_name)(*args)
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-PARAM-CONTINUATION-NOPARAM")

            instance.cleanup_code()
            instance.stop()

    def test5(self):
        "Test whether continuation methods trigger state transitions properly."
        continuation_methods = {
            "step_continuation" : ()
        }

        for method_name, args in continuation_methods.items():
            instance = iemic()

            instance.continuation.destination_0 = 1.0
            self.assertEqual(instance.get_name_of_current_state(), "INITIALIZED")

            getattr(instance, method_name)(*args)
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-PARAM-CONTINUATION-PARAM")

            instance.cleanup_code()
            instance.stop()

    def test6(self):
        "Test whether state vector operations trigger state transitions properly."
        instance = iemic()
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(),
                         "OCEAN-PARAM-CONTINUATION-NOPARAM")

        state_methods = [
            "rhs",
            "solve",
            "jacobian",
        ]

        for method_name in state_methods:
            instance.set_parameter("ocean->THCM->Starting Parameters->SPL1", 2000.0)
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-UPDATED-CONTINUATION-NOPARAM")

            state = instance.new_state()
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-PARAM-CONTINUATION-NOPARAM")


            getattr(instance, method_name)(state)
            self.assertEqual(instance.get_name_of_current_state(),
                             "OCEAN-PARAM-CONTINUATION-NOPARAM")

            del state

        instance.cleanup_code()
        instance.stop()
