import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option


class AdcircInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
    
    ADCIRC - ADvanced CIRCulation model

    .. [#] http://adcirc.org/
    
    """
    use_modules=['StoppingConditions','adcirc_interface']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'adcirc_worker'

    @remote_function
    def get_model_time():
        returns (time=0.| units.s)
        
    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)
        
    @remote_function
    def evolve_model(tend=0. | units.s):
        pass

    @remote_function(can_handle_array=True)
    def get_node_state(index=0):
        returns (eta=0. | units.m,vx=0.| units.m/units.s,vy=0.| units.m/units.s)

    @remote_function(can_handle_array=True)
    def get_node_position(index=0):
        returns (x=0.| units.m,y=0. | units.m)

    @remote_function(can_handle_array=True)
    def get_element_nodes(index=0):
        returns (n1=0,n2=0,n3=0)
