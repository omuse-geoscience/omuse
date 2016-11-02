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

class OpenIFSInterface(CodeInterface,
                     CommonCodeInterface,
                     StoppingConditionInterface,
                     LiteratureReferencesMixIn):
    """

    OpenIFS: The Open Integral Forecasting System

    """
    use_modules=['StoppingConditions','openifs_interface']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'openifs_worker'

    @remote_function
    def evolve_model():
        pass


class OpenIFS(CommonCode):

    def __init__(self,**options):
        CommonCode.__init__(self,  OpenIFSInterface(**options), **options)
        self.stopping_conditions = StoppingConditions(self)

    def commit_parameters(self):
        self.overridden().commit_parameters()

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

        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_method('EVOLVED', 'evolve_model')
