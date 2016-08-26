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

class DalesInterface(CodeInterface,
                     CommonCodeInterface,
                     StoppingConditionInterface,
                     LiteratureReferencesMixIn):
    """

    DALES - Dutch Atmospheric Large Eddy Simulation

    """
    use_modules=['StoppingConditions','dales_interface']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'dales_worker'

    @remote_function
    def get_input_file():
        returns (input_file="s")

    @remote_function
    def set_input_file(input_file="namoptions"):
        returns ()

    @remote_function
    def get_model_time():
        returns (time=0.| units.s)

    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)

    @remote_function
    def evolve_model(tend=0. | units.s):
        pass

    @remote_function
    def commit_grid():
        pass


class Dales(CommonCode):

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

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_method('!UNINITIALIZED', 'before_get_parameter')
        object.add_method('!UNINITIALIZED', 'before_set_parameter')
        object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('EDIT','RUN','commit_grid')

        #~ object.set_initial_state('UNINITIALIZED')
        #~ object.add_transition('!STOPPED', 'END', 'cleanup_code')
        #~ object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        #~ object.add_transition('END', 'STOPPED', 'stop', False)
        #~ object.add_method('STOPPED', 'stop')

        object.add_method('INITIALIZED', 'before_set_interface_parameter')
        object.add_method('INITIALIZED', 'set_input_file')

        for state in ["RUN","EDIT","EVOLVED"]:
          object.add_method(state,"get_model_time")
          object.add_method(state,"get_timestep")

        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_method('EVOLVED', 'evolve_model')
