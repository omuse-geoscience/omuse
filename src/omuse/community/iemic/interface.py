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
    def set_default_params():
        returns ()

    @remote_function
    def get_ocean_params():
        returns (ocean_params="s")

    @remote_function
    def set_ocean_params(ocean_params="s"):
        returns ()

    @remote_function
    def get_continuation_params():
        returns (continuation_params="s")

    @remote_function
    def set_continuation_params(continuation_params="s"):
        returns ()

    @remote_function
    def commit_parameters():
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
