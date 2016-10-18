import os.path
from omuse.units import units
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system
from amuse.community.interface.stopping_conditions import StoppingConditionInterface, StoppingConditions
from amuse import datamodel

from amuse.units import trigo

import numpy

class DalesInterface(CodeInterface,
                     CommonCodeInterface,
                     StoppingConditionInterface,
                     LiteratureReferencesMixIn):
    """

    DALES - Dutch Atmospheric Large Eddy Simulation

    """
    use_modules=["StoppingConditions","dales_interface"]

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return "dales_worker"

    @remote_function
    def get_input_file():
        returns (input_file="s")

    @remote_function
    def set_input_file(input_file="namoptions"):
        pass

    @remote_function
    def get_model_time():
        returns (time=0.| units.s)

    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)

    @remote_function
    def commit_grid():        
        pass

# getter functions for vertical profiles - slab averages
# these take a dummy array as input, and return output of the same length
    @remote_function(must_handle_array=True)
    def get_profile_field(k=0):
        returns (out=0.| units.K)

    @remote_function(must_handle_array=True)
    def get_profile_U_(k=0):
        returns (out=0.| units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_V_(k=0):
        returns (out=0.| units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_W_(k=0):
        returns (out=0.| units.m / units.s)
        
    @remote_function(must_handle_array=True)
    def get_profile_THL_(k=0):
        returns (out=0.| units.K)
        
    @remote_function(must_handle_array=True)
    def get_profile_QT_(k=0):
        returns (out=0.)

    @remote_function(must_handle_array=True)
    def get_zf_(k=0):
        returns (out=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_zh_(k=0):
        returns (out=0. | units.m)

# # #
    @remote_function()
    def get_params_grid():
        returns (i=1, j=1, k=1, xsize=1.0 | units.m, ysize=1.0 | units.m)

    @remote_function(must_handle_array=True)
    def get_layer_field(i=0,j=0,k=0):
        returns (temp=0.| units.K)

    @remote_function(must_handle_array=True)
    def get_volume_field(i=0,j=0,k=0):
        returns (temp=0.| units.K)

    @remote_function
    def evolve_model(tend=0. | units.s):
        pass

    
    
class Dales(CommonCode):
    # these defs match those in daleslib.f90
    FIELDID_U=1
    FIELDID_V=2
    FIELDID_W=3
    FIELDID_THL=4
    FIELDID_QT=5
    
    # grid size
    itot = None
    jtot = None
    k    = None
    xsize = None
    ysize = None

    

    
    def __init__(self,**options):
        CommonCode.__init__(self,  DalesInterface(**options), **options)
        self.stopping_conditions = StoppingConditions(self)

    def commit_parameters(self):
        self.overridden().commit_parameters()

    def define_parameters(self, object):
        object.add_default_form_parameter(
            "input_file",
            "set the input file path",
            "namoptions"
        )


    def commit_grid(self):        
        self.overridden().commit_grid()
        self.get_params()
        
    def define_state(self, object):
        object.set_initial_state("UNINITIALIZED")
        object.add_transition("UNINITIALIZED", "INITIALIZED", "initialize_code")
        object.add_method("UNINITIALIZED", "before_get_parameter")
        object.add_method("!UNINITIALIZED", "before_set_parameter")
        object.add_method("END", "before_get_parameter")
        object.add_transition("!UNINITIALIZED!STOPPED", "END", "cleanup_code")
        object.add_transition("END", "STOPPED", "stop", False)
        object.add_method("STOPPED", 'stop')

        object.add_transition("INITIALIZED","EDIT","commit_parameters")
        object.add_transition("EDIT","RUN","commit_grid")

        object.add_method("INITIALIZED", "before_set_interface_parameter")
        object.add_method("INITIALIZED", "set_input_file")

        for state in ["RUN","EDIT","EVOLVED"]:
          object.add_method(state,"get_model_time")
          object.add_method(state,"get_timestep")
          object.add_method(state,"get_profile_field")
          object.add_method(state,"get_layer_field")
          object.add_method(state,"get_volume_field")

        object.add_transition("RUN", "EVOLVED", "evolve_model", False)
        object.add_method("EVOLVED", "evolve_model")


    # wrapping functions for hiding the dummy array passed to getter functions
    def get_profile_U(self):
        return self.get_profile_U_(numpy.zeros(self.k))

    def get_profile_V(self):
        return self.get_profile_V_(numpy.zeros(self.k))

    def get_profile_W(self):
        return self.get_profile_W_(numpy.zeros(self.k))
        
    def get_profile_THL(self):
        return self.get_profile_THL_(numpy.zeros(self.k))
        
    def get_profile_QT(self):
        return self.get_profile_QT_(numpy.zeros(self.k))

    def get_zf(self):
        return self.get_zf_(numpy.zeros(self.k))

    def get_zh(self):
        return self.get_zh_(numpy.zeros(self.k))


    # get parameters from the fortran code, store them in the Dales interface object
    def get_params(self):
        self.itot,self.jtot,self.k,self.xsize,self.ysize = self.get_params_grid()
        
