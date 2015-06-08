import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

class SwanInterface(CodeInterface, 
                      CommonCodeInterface,
                      LiteratureReferencesMixIn):
    """
    
    SWAN - 

    .. [#] swanmodel.sf.net
    
    """
    use_modules=['StoppingConditions','swan_interface']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'swan_worker'
