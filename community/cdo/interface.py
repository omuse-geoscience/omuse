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
    def set_weights_file(filename='s'):
        returns ()
    @remote_function()
    def set_src_grid_file(filename='s'):
        returns ()
    @remote_function()
    def set_dst_grid_file(filename='s'):
        returns ()
    @remote_function()
    def get_weights_file():
        returns (filename='s')
    @remote_function()
    def get_src_grid_file():
        returns (filename='s')
    @remote_function()
    def get_dst_grid_file():
        returns (filename='s')


    @remote_function()
    def initialize_code():
        returns ()
    @remote_function()
    def commit_parameters():
        returns ()
    @remote_function()
    def cleanup_code():
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




    # debug functions

    @remote_function()
    def get_num_links():
        returns (num_links=0)
    @remote_function(must_handle_array=True)
    def get_remap_links(i=0):
        returns (src_address=0, dst_address=0, remap_weights=0.)
    @remote_function()
    def print_info():
        returns ()
    





class CDORemapper(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self, CDOInterface(), **options)


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_weights_file",
            "set_weights_file",
            "weights_file",
            "Specify the filename of weigths file to be used for remapping",
            default_value = ''
        )
        object.add_method_parameter(
            "get_src_grid_file",
            "set_src_grid_file",
            "src_grid_file",
            "Specify the filename of src grid file to be used for remapping",
            default_value = ''
        )
        object.add_method_parameter(
            "get_dst_grid_file",
            "set_dst_grid_file",
            "dst_grid_file",
            "Specify the filename of dst grid file to be used for remapping",
            default_value = ''
        )


    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')

        object.add_transition('!UNINITIALIZED!INITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        object.add_method('INITIALIZED', 'set_weights_file')
        object.add_method('INITIALIZED', 'set_src_grid_file')
        object.add_method('INITIALIZED', 'set_dst_grid_file')

        object.add_transition('INITIALIZED', 'RUN', 'commit_parameters')


        for state in ["INITIALIZED","RUN"]:
            object.add_method(state, 'get_weights_file')
            object.add_method(state, 'get_src_grid_file')
            object.add_method(state, 'get_dst_grid_file')

        
        object.add_method("RUN", 'set_src_grid_values')
        object.add_method("RUN", 'get_src_grid_values')
        object.add_method("RUN", 'set_dst_grid_values')
        object.add_method("RUN", 'get_dst_grid_values')
        object.add_method("RUN", 'perform_remap')

        object.add_method("RUN", 'get_src_grid_size')
        object.add_method("RUN", 'get_src_grid_dims')
        object.add_method("RUN", 'get_dst_grid_size')
        object.add_method("RUN", 'get_dst_grid_dims')


