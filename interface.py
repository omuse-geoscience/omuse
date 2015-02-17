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

    @legacy_function
    def test_amuse():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function
