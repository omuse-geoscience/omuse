from amuse.community import *

class DFlowFMInterface(CodeInterface):

    use_modules=["dflowfm_omuse"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="dflowfm_worker", **keyword_arguments)
    
    @remote_function
    def initialize():
        returns ()
        
    @remote_function
    def evolve_model(tend=0. | units.s):
        returns ()
    
    @remote_function
    def get_model_time():
        returns (model_time = 0. | units.s)
    
class DFlowFM(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  DFlowFMInterface(**options), **options)
    
    def define_properties(self, object):
        object.add_property('get_model_time', public_name = "model_time")
