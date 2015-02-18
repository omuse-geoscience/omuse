import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

#~ IN=LegacyFunctionSpecification.IN
#~ OUT=LegacyFunctionSpecification.OUT
#~ 
#~ functions=dict(
 #~ evolve_model={ "parameters": ["tend",dict(dtype="double",direction=IN,unit=units.s),],
                #~ "result_type": "i" } 
#~ )

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

    @legacy_function
    def get_model_time():
        function = LegacyFunctionSpecification()
        function.addParameter("tend", dtype='d', direction=function.OUT ,unit=units.s)
        function.result_type = 'i'
        return function
    @legacy_function
    def get_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter("dt", dtype='d', direction=function.OUT ,unit=units.s)
        function.result_type = 'i'
        return function
    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter("tend", dtype='d', direction=function.IN ,unit=units.s)
        function.result_type = 'i'
        return function


