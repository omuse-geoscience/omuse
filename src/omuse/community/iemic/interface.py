from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.datamodel import CartesianGrid

from omuse.units import units



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
    def test_grid(logFile="s"):
        returns ()

    @remote_function
    def step():
        returns ()

    @remote_function
    def cleanup_code():
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

class iemic(InCodeComponentImplementation):
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  iemicInterface(**options), **options)

    def define_state(self, handler):
        handler.set_initial_state('UNINITIALIZED')
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize')
        handler.add_transition('INITIALIZED', 'PARAM', 'commit_parameters')
        
        for state in ["PARAM"]:
            for method in ["get_u", "get_v", "get_w", "get_p", "get_t", "get_s",
                "get_nrange", "get_mrange", "get_lrange","step"]:
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
