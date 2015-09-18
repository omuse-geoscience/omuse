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
    def get_src_grid_dims():
        returns (x=0, y=0)
    @remote_function()
    def get_dst_grid_dims():
        returns (x=0, y=0)

    @remote_function()
    def get_num_links():
        returns (num_links=0)
    @remote_function(must_handle_array=True)
    def get_remap_links(i=0):
        returns (src_address=0, dst_address=0, remap_weights=0.)


    @remote_function()
    def get_remap_file():
        returns (filename='s')
    @remote_function()
    def set_remap_file(filename='s'):
        returns ()



    @remote_function()
    def perform_remap():
        returns ()

    
    @remote_function(must_handle_array=True)
    def set_src_grid_values(i=0, src=0.):
        returns ()
    @remote_function(must_handle_array=True)
    def set_dst_grid_values(i=0, dst=0.):
        returns ()
    @remote_function(must_handle_array=True)
    def get_src_grid_values(i=0):
        returns (src=0.)
    @remote_function(must_handle_array=True)
    def get_dst_grid_values(i=0):
        returns (dst=0.)




    @remote_function(must_handle_array=True)
    def get_src_grid_mask(i=0, j=0):
        returns (mask=0)

    





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
 
