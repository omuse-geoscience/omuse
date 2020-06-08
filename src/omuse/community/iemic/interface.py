from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.datamodel import CartesianGrid

from omuse.units import units

from remotestatevector import RemoteStateVector

class iemicInterface(CodeInterface,CommonCodeInterface):
    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="iemic_worker", **keyword_arguments)

    @remote_function
    def initialize():
        returns ()

    @remote_function
    def commit_parameters():
        returns ()

    @remote_function
    def commit_continuation_parameters():
        returns ()

    @remote_function
    def test_grid(logFile="s"):
        returns ()

    @remote_function
    def step_continuation():
        returns ()

    @remote_function
    def cleanup_code():
        returns ()

    @remote_function
    def load_xml_parameters(set_name="", path=""):
        returns ()

    @remote_function
    def save_xml_parameters(set_name="", path=""):
        returns ()

    @remote_function(must_handle_array=True)
    def get_u(i=0, j=0, k=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_v(i=0, j=0, k=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_w(i=0, j=0, k=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_p(i=0, j=0, k=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_t(i=0, j=0, k=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_s(i=0, j=0, k=0):
        returns (var=0.)

    @remote_function
    def get_nrange():
        returns(nmin=0,nmax=0)

    @remote_function
    def get_mrange():
        returns(nmin=0,nmax=0)

    @remote_function
    def get_lrange():
        returns(nmin=0,nmax=0)

    @remote_function
    def get_num_parameter_sets():
        returns(num=0)

    @remote_function
    def get_parameter_set_name(i=0):
        returns(name="")

    @remote_function
    def get_num_parameters(set_name="", param_name=""):
        returns(num=0)

    @remote_function
    def get_parameter_name(set_name="", param_name="", i=0):
        returns(name="")

    @remote_function
    def get_parameter_type(set_name="", param_name=""):
        returns(name="")

    @remote_function
    def get_bool_parameter(set_name="", param_name=""):
        returns(value=False)

    @remote_function
    def set_bool_parameter(set_name="", param_name="", value=False):
        returns()

    @remote_function
    def get_state_norm():
        returns(norm=0.)

    @remote_function
    def _get_state_norm(src=0):
        returns(norm=0.)


    @remote_function
    def get_default_bool_parameter(set_name="", param_name=""):
        returns(value=False)

    #@remote_function
    #def get_char_parameter(set_name="", param_name=""):
    #    returns(value='c')

    #@remote_function
    #def set_char_parameter(set_name="", param_name="", value='c'):
    #    returns()

    #@remote_function
    #def get_default_char_parameter(set_name="", param_name=""):
    #    returns(value='c')

    @remote_function
    def get_double_parameter(set_name="", param_name=""):
        returns(value=0.)

    @remote_function
    def set_double_parameter(set_name="", param_name="", value=0.):
        returns()

    @remote_function
    def get_default_double_parameter(set_name="", param_name=""):
        returns(value=0.)

    @remote_function
    def get_int_parameter(set_name="", param_name=""):
        returns(value=0)

    @remote_function
    def set_int_parameter(set_name="", param_name="", value=0):
        returns()

    @remote_function
    def get_default_int_parameter(set_name="", param_name=""):
        returns(value=0)

    @remote_function
    def get_string_parameter(set_name="", param_name=""):
        returns(value="")

    @remote_function
    def set_string_parameter(set_name="", param_name="", value=""):
        returns()

    @remote_function
    def get_default_string_parameter(set_name="", param_name=""):
        returns(value="")

    @remote_function(can_handle_array=True)
    def _new_state():
        returns(index=0)
    
    @remote_function
    def _remove_state(index=0):
        returns()

    @remote_function
    def _mul_state(src=0, fac=0.):
        returns()

    @remote_function
    def _add_state(s1=0, s2=0):
        returns()

    @remote_function
    def _get_rhs(src=0, target=0):
        returns()

    @remote_function
    def _copy_state(src=0, target=0):
        returns()

    @remote_function
    def _solve(src=0, target=0):
        returns()

    @remote_function
    def _jacobian(src=0):
        returns()
        

class iemic(InCodeComponentImplementation):
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  iemicInterface(**options), **options)

    def define_state(self, handler):
        handler.set_initial_state('UNINITIALIZED')
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize')
        handler.add_transition('INITIALIZED', 'PARAM', 'commit_parameters')
        handler.add_transition('PARAM', 'PARAMC', 'commit_continuation_parameters')
        
        for state in ["PARAM", "PARAMC"]:
            for method in ["get_u", "get_v", "get_w", "get_p", "get_t", "get_s",
                "get_nrange", "get_mrange", "get_lrange", "_new_state"]:
                handler.add_method(state, method)

        for state in ["PARAMC"]:
            for method in ["step", "run_continuation"]:
                handler.add_method(state, method)

    def _grid_range(self):
        return self.get_nrange()+self.get_mrange()+self.get_lrange()
    
    def define_grids(self, handler):
        #~ code_generator.generate_grid_definitions(object)
        
        # ocean dynamic variables p-grid: po, dpo_dt, vorticity
        handler.define_grid('grid',axes_names = ["lon", "lat"], grid_class=CartesianGrid)
        handler.set_grid_range('grid', '_grid_range')
        #~ handler.add_getter('grid', 'get_grid_position', names=["lon", "lat"])
        handler.add_getter('grid', 'get_u', names=["u_velocity"])
        handler.add_getter('grid', 'get_v', names=["v_velocity"])
        handler.add_getter('grid', 'get_w', names=["w_velocity"])


    def _generate_parameter_setter_getters(self, paramSet, sublist, paramName, paramType):
        # define getter (closure or partial..)
        def getter():
          return getattr(self, "get_"+paramType+"_parameter")(paramSet, "?".join(sublist+[paramName]))
        # define setter 
        def setter(val):
          return getattr(self, "set_"+paramType+"_parameter")(paramSet, "?".join(sublist+[paramName]), val)
        return getter,setter


    def define_parameter_set(self, handler, paramSet, sublist=[]):
        paramCount = self.get_num_parameters(paramSet, "?".join(sublist))
        
        allowed_types=["bool", "string", "double", "int"]
        
        for j in range(0, paramCount):
            paramName = self.get_parameter_name(paramSet, "?".join(sublist) , j)
            paramType = self.get_parameter_type(paramSet, "?".join(sublist + [paramName]))
            if paramType=="ParameterList":
                self.define_parameter_set(handler,paramSet, sublist + [paramName])
            else:
                if paramType=="char":
                    # not handled yet?
                    continue
                if paramType not in allowed_types:
                    raise Exception("encountered unknown parameter type")
                # normalized parameter name , set name
                longname="__".join([paramSet] + sublist + [paramName])
                longname=longname.replace(" ", "_").replace("-","_")
                parameter_set_name="__".join([paramSet] + sublist)
              
                getter,setter=self._generate_parameter_setter_getters(paramSet, sublist, paramName, paramType)

                # name of getter
                name_of_getter="get_"+longname                
                # name of setter
                name_of_setter="set_"+longname

                # add getter to interface
                setattr(self, name_of_getter, getter)
                # add setter to interface
                setattr(self, name_of_setter, setter)

                # get default
                default=getattr(self, "get_default_"+paramType+"_parameter")(paramSet, "?".join(sublist+[paramName]))
                # define parameter
                handler.add_method_parameter(
                  name_of_getter,
                  name_of_setter,
                  paramName.replace(" ","_").replace("-","_"),
                  "generated parameter",
                  default_value = default,
                  parameter_set=parameter_set_name.replace(" ", "_").replace("-","_")
               )
 

    def define_parameters(self, handler):
        
        paramSetCount = self.get_num_parameter_sets()

        for i in range(0, paramSetCount):
            paramSet = self.get_parameter_set_name(i)

            self.define_parameter_set(handler,paramSet)

    def new_state(self):
        return RemoteStateVector(self)

    def rhs(self, state):
        result=self.new_state()
        self._get_rhs(state._id, result._id)
        return result

    def solve(self, state):
        result=self.new_state()
        self._solve(state._id, result._id)
        return result
        
    def jacobian(self, state):
        self._jacobian(state._id)
