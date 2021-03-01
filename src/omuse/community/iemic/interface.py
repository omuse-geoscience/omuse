from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.datamodel import CartesianGrid

from omuse.units import units

from remotestatevector import RemoteStateVector
from implicit_utils import time_stepper

class iemicInterface(CodeInterface,CommonCodeInterface):
    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="iemic_worker", **keyword_arguments)

        self._model_time=0

    def set_ocean_transition(self):
        pass

    def set_continuation_transition(self):
        pass

    @remote_function
    def set_log_file(logFile="s"):
        returns ()

    @remote_function
    def set_output_file(outFile="s"):
        returns ()

    @remote_function
    def initialize():
        returns ()

    @remote_function
    def commit_parameters():
        returns ()

    @remote_function
    def recommit_parameters():
        returns ()

    @remote_function
    def commit_continuation_parameters():
        returns ()

    @remote_function
    def recommit_continuation_parameters():
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
    def _load_xml_parameters(set_name="", path=""):
        returns ()

    @remote_function
    def save_xml_parameters(set_name="", path=""):
        returns ()

    @remote_function(must_handle_array=True)
    def get_land_mask(n=0, m=0, l=0):
        returns (var=0)

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

    @remote_function(must_handle_array=True)
    def get_u_(i=0, j=0, k=0, sindex=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_v_(i=0, j=0, k=0,sindex=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_w_(i=0, j=0, k=0,sindex=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_p_(i=0, j=0, k=0,sindex=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_t_(i=0, j=0, k=0,sindex=0):
        returns (var=0.)

    @remote_function(must_handle_array=True)
    def get_s_(i=0, j=0, k=0,sindex=0):
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
    def get_real_pos(xIn=0,yIn=0,zIn=0):
        returns(xOut=0.0,yOut=0.0,zOut=0.0)

    @remote_function
    def _get_parameter_type(set_name="", param_name=""):
        returns(name="")

    @remote_function
    def _get_bool_parameter(set_name="", param_name=""):
        returns(value=False)

    @remote_function
    def _set_bool_parameter(set_name="", param_name="", value=False):
        returns()

    @remote_function
    def _get_default_bool_parameter(set_name="", param_name=""):
        returns(value=False)

    @remote_function
    def get_state_norm():
        returns(norm=0.)

    @remote_function
    def _get_state_norm(src=0):
        returns(norm=0.)

    @remote_function
    def _get_char_parameter(set_name="", param_name=""):
        returns(value='c')

    @remote_function
    def _set_char_parameter(set_name="", param_name="", value='c'):
        returns()

    @remote_function
    def _get_default_char_parameter(set_name="", param_name=""):
        returns(value='c')

    @remote_function
    def _get_double_parameter(set_name="", param_name=""):
        returns(value=0.)

    @remote_function
    def _set_double_parameter(set_name="", param_name="", value=0.):
        returns()

    @remote_function
    def _get_default_double_parameter(set_name="", param_name=""):
        returns(value=0.)

    @remote_function
    def _get_int_parameter(set_name="", param_name=""):
        returns(value=0)

    @remote_function
    def _set_int_parameter(set_name="", param_name="", value=0):
        returns()

    @remote_function
    def _get_default_int_parameter(set_name="", param_name=""):
        returns(value=0)

    @remote_function
    def _get_string_parameter(set_name="", param_name=""):
        returns(value="")

    @remote_function
    def _set_string_parameter(set_name="", param_name="", value=""):
        returns()

    @remote_function
    def _get_default_string_parameter(set_name="", param_name=""):
        returns(value="")

    @remote_function(can_handle_array=True)
    def _new_state():
        returns(index=0)

    @remote_function
    def _get_model_state(target=0):
        returns()

    @remote_function
    def _set_model_state(src=0):
        returns()

    @remote_function
    def _remove_state(index=0):
        returns()
    
    @remote_function
    def _to_str(src=0):
        returns(out="")

    @remote_function
    def _mul_state(src=0, fac=0.):
        returns()

    @remote_function
    def _update_state(s1=0, s2=0, scal=1.0):
        returns()

    @remote_function
    def _dot(s1=0, s2=0):
        returns(val=0.)

    @remote_function
    def _length(src=0):
        returns(val=0)

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

    @remote_function
    def _jacobian_with_mass_matrix(src=0, sigma=0.0):
        returns()

    @remote_function
    def _apply_mass_matrix(src=0, target=0):
        returns()

    @remote_function
    def _get_psi_m(src=0):
        returns(psi_min=0.0, psi_max=0.0)

class iemic(InCodeComponentImplementation):
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, iemicInterface(**options), **options)

    def define_state(self, handler):
        ocean_sub_states = ["PARAM", "UPDATED"]
        continuation_sub_states = ocean_sub_states + ["NOPARAM"]

        handler.set_initial_state('UNINITIALIZED')
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize')
        handler.add_transition('INITIALIZED', 'OCEAN-PARAM-CONTINUATION-NOPARAM',
                                'commit_parameters')

        handler.add_transition('INITIALIZED', 'INITIALIZED', 'set_ocean_transition')
        handler.add_transition('INITIALIZED', 'INITIALIZED', 'set_continuation_transition')

        for sub_state in ocean_sub_states:
            handler.add_transition('OCEAN-' + sub_state + '-CONTINUATION-NOPARAM',
                                   'OCEAN-' + sub_state + '-CONTINUATION-PARAM',
                                   'commit_continuation_parameters')

            handler.add_transition('OCEAN-' + sub_state + '-CONTINUATION-PARAM',
                                   'OCEAN-' + sub_state + '-CONTINUATION-UPDATED',
                                    'set_continuation_transition')

            handler.add_transition('OCEAN-' + sub_state + '-CONTINUATION-NOPARAM',
                                   'OCEAN-' + sub_state + '-CONTINUATION-NOPARAM',
                                    'set_continuation_transition')

            handler.add_transition('OCEAN-' + sub_state + '-CONTINUATION-UPDATED',
                                   'OCEAN-' + sub_state + '-CONTINUATION-PARAM',
                                   'recommit_continuation_parameters')

        for sub_state in continuation_sub_states:
            handler.add_transition('OCEAN-PARAM-CONTINUATION-' + sub_state,
                                   'OCEAN-UPDATED-CONTINUATION-' + sub_state,
                                   'set_ocean_transition')

            handler.add_transition('OCEAN-UPDATED-CONTINUATION-' + sub_state,
                                   'OCEAN-PARAM-CONTINUATION-' + sub_state,
                                    'recommit_parameters')

        ocean_methods = [
            "get_u", "get_v", "get_w", "get_p", "get_t", "get_s", "get_nrange",
            "get_mrange", "get_lrange", "_new_state", "_get_rhs", "_solve",
            "_jacobian", "_jacobian_with_mass_matrix", "_apply_mass_matrix",
            "_set_model_state", "get_state_norm",
            "_get_model_state", "_get_psi_m", "get_real_pos",
            "get_land_mask"
        ]

        for sub_state in continuation_sub_states:
            for method in ocean_methods:
                handler.add_method('OCEAN-PARAM-CONTINUATION-' + sub_state, method)

        handler.add_method('OCEAN-PARAM-CONTINUATION-PARAM',
                           'step_continuation')

    def _grid_range(self, **kwargs):
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

    def _specify_grid(self, definition, index=0):
        #~ handler.define_grid('grid',axes_names = ["lon", "lat"], grid_class=CartesianGrid)
        definition.set_grid_range('_grid_range')
        #~ handler.add_getter('grid', 'get_grid_position', names=["lon", "lat"])
        definition.add_getter( 'get_u_', names=["u_velocity"])
        definition.add_getter( 'get_v_', names=["v_velocity"])
        definition.add_getter( 'get_w_', names=["w_velocity"])
        definition.define_extra_keywords({'sindex':index})


    def _generate_parameter_setter_getters(self, paramSet, sublist, paramName, paramType):
        # define getter (closure or partial..)
        def getter():
          return getattr(self, "_get_"+paramType+"_parameter")(paramSet, "->".join(sublist+[paramName]))
        # define setter
        def setter(val):
          return getattr(self, "set_"+paramType+"_parameter")(paramSet, "->".join(sublist+[paramName]), val)
        return getter,setter


    def define_parameter_set(self, handler, paramSet, sublist=[]):
        paramCount = self.get_num_parameters(paramSet, "->".join(sublist))

        allowed_types=["bool", "char", "string", "double", "int"]
        
        for j in range(0, paramCount):
            paramName = self.get_parameter_name(paramSet, "->".join(sublist) , j)
            paramType = self._get_parameter_type(paramSet, "->".join(sublist + [paramName]))
            if paramType=="ParameterList":
                self.define_parameter_set(handler,paramSet, sublist + [paramName])
            else:
                if paramType not in allowed_types:
                    raise Exception("encountered unknown parameter type")
                # normalized parameter name , set name
                longname="__".join([paramSet] + sublist + [paramName])
                longname=longname.replace(" ", "_").replace("-","_").replace("/","_")
                parameter_set_name="__".join([paramSet] + sublist)

                getter,setter=self._generate_parameter_setter_getters(paramSet, sublist, paramName, paramType)

                # name of getter
                name_of_getter="_get_"+longname
                # name of setter
                name_of_setter="set_"+longname

                # add getter to interface
                setattr(self, name_of_getter, getter)
                # add setter to interface
                setattr(self, name_of_setter, setter)

                # get default
                default=getattr(self, "_get_default_"+paramType+"_parameter")(paramSet, "->".join(sublist+[paramName]))
                # define parameter
                handler.add_method_parameter(
                  name_of_getter,
                  name_of_setter,
                  paramName.replace(" ","_").replace("-","_").replace("/","_"),
                  "generated parameter",
                  default_value = default,
                  parameter_set=parameter_set_name.replace(" ", "_").replace("-","_").replace("/","_")
                )
                #define also in parameters for flat access 
                handler.add_method_parameter(
                  name_of_getter,
                  name_of_setter,
                  longname,
                  "generated parameter",
                  default_value = default,
                )


    def define_parameters(self, handler):
        paramSetCount = self.get_num_parameter_sets()

        for i in range(0, paramSetCount):
            paramSet = self.get_parameter_set_name(i)

            self.define_parameter_set(handler,paramSet)
            
        handler.add_interface_parameter( "timestepping_theta", "timestepping theta parameter [0-1]", 0.5)
        handler.add_interface_parameter( "timestepping_dt", "timestep for timestepper", 1.)

    def new_state(self):
        return RemoteStateVector(self)

    def get_state(self):
        result = self.new_state()
        self._get_model_state(result._id)
        return result

    def set_state(self, state):
        self._set_model_state(state._id)

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

    def jacobian_with_mass_matrix(self, state, sigma):
        self._jacobian_with_mass_matrix(state._id, sigma)

    def apply_mass_matrix(self, state):
        result=self.new_state()
        self._apply_mass_matrix(state._id, result._id)
        return result

    def get_psi_m(self, state):
        return self._get_psi_m(state._id)

    def get_grid(self, state):
        return self._create_new_grid(self._specify_grid, index=state._id)

    def parameterset_parameters(self, paramSet, sublist=[]):
        paramCount = self.get_num_parameters(paramSet, "->".join(sublist))

        for j in range(0, paramCount):
            paramName = self.get_parameter_name(paramSet, "->".join(sublist) , j)
            paramComponents = sublist + [paramName]
            paramType = self._get_parameter_type(paramSet, "->".join(paramComponents))
            if paramType=="ParameterList":
                for param in self.parameterset_parameters(paramSet, paramComponents):
                    yield param
            else:
                yield "->".join([paramSet] + paramComponents)

    def parameter_names(self):
        paramSetCount = self.get_num_parameter_sets()

        for i in range(0, paramSetCount):
            paramSet = self.get_parameter_set_name(i)

            for param in self.parameterset_parameters(paramSet):
                yield param

    def get_parameter_type(self, param_name):
        set_name, param_name = param_name.split("->", 1)
        return self._get_parameter_type(set_name, param_name)

    def get_parameter(self, name):
        type = self.get_parameter_type(name)
        set_name, param_name = name.split("->", 1)
        return getattr(self, '_get_' + type + '_parameter')(set_name, param_name)

    def get_default_parameter(self, name):
        type = self.get_parameter_type(name)
        set_name, param_name = name.split("->", 1)
        return getattr(self, '_get_default_' + type + '_parameter')(set_name, param_name)

    def set_parameter(self, name, value):
        type_ = self.get_parameter_type(name)
        set_name, param_name = name.split("->", 1)
        getattr(self, 'set_' + type_ + '_parameter')(set_name, param_name, value)

    def set_bool_parameter(self, set_name="", param_name="", value=False):
        result = self._set_bool_parameter(set_name, param_name, value)
        getattr(self, "set_" + set_name.lower() + "_transition")()
        return result

    def set_char_parameter(self, set_name="", param_name="", value='V'):
        if len(value) != 1:
            raise Exception("char parameter must have length 1")

        result = self._set_char_parameter(set_name, param_name, value)
        getattr(self, "set_" + set_name.lower() + "_transition")()
        return result

    def set_double_parameter(self, set_name="", param_name="", value=0.):
        result = self._set_double_parameter(set_name, param_name, value)
        getattr(self, "set_" + set_name.lower() + "_transition")()
        return result

    def set_int_parameter(self, set_name="", param_name="", value=0):
        result = self._set_int_parameter(set_name, param_name, value)
        getattr(self, "set_" + set_name.lower() + "_transition")()
        return result

    def set_string_parameter(self, set_name="", param_name="", value=""):
        result =  self._set_string_parameter(set_name, param_name, value)
        getattr(self, "set_" + set_name.lower() + "_transition")()
        return result

    def load_xml_parameters(set_name="", path=""):
        result =  self._load_xml_parameters(set_name, path)
        getattr(self, "set_" + set_name.lower() + "_transition")()
        return result

    def grid_pos_to_real_pos(self,coordinate):
        return self.get_real_pos(*coordinate)

    def get_surface_mask(self):
        n = self.get_parameter("Ocean->THCM->Global Grid-Size n")
        m = self.get_parameter("Ocean->THCM->Global Grid-Size m")

        indices = [(_n,_m,0) for _m in range(0,m) for _n in range(0,n)]
        ns, ms, ls = zip(*indices)

        return self.get_land_mask(list(ns), list(ms), list(ls)).reshape(n, m)

    def evolve_model(self,tend):
        x=self.get_state()
        t=0
        tmax=tend-self._model_time
        theta=self.parameters.timestepping_theta
        dt=self.parameters.timestepping_dt
        dt=min(dt, tmax)
        xnew=time_stepper(self, x, theta, dt, tmax)
        self.set_state(xnew)
        self._model_time=tmax
