import numpy
from amuse.community.interface.common import CommonCode, CommonCodeInterface 
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface

from amuse.units import units

import subprocess 

class QGCMInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
        .. [#] Hogg et al. 2006, 2003

    """
    
    use_modules = ['qgcm_interface',]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="q-gcm_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
