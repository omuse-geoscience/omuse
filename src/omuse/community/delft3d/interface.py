from amuse.community import *

class DFlowFMInterface(CodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="dflowfm_worker", **keyword_arguments)
    
    @legacy_function
    def initialize():
        function = LegacyFunctionSpecification()  
        return function
        
    
class DFlowFM(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  DFlowFMInterface(**options), **options)
    
