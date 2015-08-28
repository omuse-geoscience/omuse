from amuse.community import *

from amuse.community.interface.common import CommonCodeInterface, CommonCode

from omuse.units import units

class CDOInterface(CodeInterface):
    
    """
    CDO - Climate Data Operators

    .. [#] https://code.zmaw.de/projects/cdo

    """
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(), **keyword_arguments)

    def name_of_the_worker(self):
        return 'remapper_worker'

    #example of must_handle_array
    #@remote_function(must_handle_array=True)
    #def get_node_wind_stress(i=0,j=0):
    #    returns (tau_x=0. | units.Pa,tau_y=0. | units.Pa)

    #example without must_handle_array
    #@remote_function
    #def get_number_of_nodes():
    #    returns (n_nodes=0)


    @remote_function()
    def get_src_grid_size():
        returns (size=0)
    @remote_function()
    def get_dst_grid_size():
        returns (size=0)

    @remote_function()
    def get_remap_file():
        returns (filename='s')
    @remote_function()
    def set_remap_file(filename='s'):
        returns ()


    





class CDORemapper(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self, CDOInterface(), **options)


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_remap_file",
            "set_remap_file",
            "remap_file",
            "Specify the filename of weigths file to be used for remapping",
            default_value = ''
        )
 
