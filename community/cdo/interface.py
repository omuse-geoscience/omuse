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

    @remote_function()
    def get_src_grid_size():
        returns (size=0)
    @remote_function()
    def get_dst_grid_size():
        returns (size=0)
    @remote_function()
    def set_src_grid_size(size=0):
        returns ()
    @remote_function()
    def set_dst_grid_size(size=0):
        returns ()

    @remote_function()
    def get_src_grid_corners():
        returns (corners=0)
    @remote_function()
    def get_dst_grid_corners():
        returns (corners=0)
    @remote_function()
    def set_src_grid_corners(corners=0):
        returns ()
    @remote_function()
    def set_dst_grid_corners(corners=0):
        returns ()

    @remote_function()
    def get_src_grid_dims():
        returns (x=0, y=0)
    @remote_function()
    def get_dst_grid_dims():
        returns (x=0, y=0)
    @remote_function()
    def set_src_grid_dims(x=0, y=0):
        returns ()
    @remote_function()
    def set_dst_grid_dims(x=0, y=0):
        returns ()


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


    @remote_function(must_handle_array=True)
    def set_src_grid_center_lat(i=0, center_lat=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def set_dst_grid_center_lat(i=0, center_lat=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def get_src_grid_center_lat(i=0):
        returns (center_lat=0. | units.rad)
    @remote_function(must_handle_array=True)
    def get_dst_grid_center_lat(i=0):
        returns (center_lat=0. | units.rad)

    @remote_function(must_handle_array=True)
    def set_src_grid_center_lon(i=0, center_lon=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def set_dst_grid_center_lon(i=0, center_lon=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def get_src_grid_center_lon(i=0):
        returns (center_lon=0. | units.rad)
    @remote_function(must_handle_array=True)
    def get_dst_grid_center_lon(i=0):
        returns (center_lon=0. | units.rad)

    @remote_function(must_handle_array=True)
    def set_src_grid_corner_lat(i=0, corner_lat=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def set_dst_grid_corner_lat(i=0, corner_lat=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def get_src_grid_corner_lat(i=0):
        returns (corner_lat=0. | units.rad)
    @remote_function(must_handle_array=True)
    def get_dst_grid_corner_lat(i=0):
        returns (corner_lat=0. | units.rad)

    @remote_function(must_handle_array=True)
    def set_src_grid_corner_lon(i=0, corner_lon=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def set_dst_grid_corner_lon(i=0, corner_lon=0. | units.rad):
        returns ()
    @remote_function(must_handle_array=True)
    def get_src_grid_corner_lon(i=0):
        returns (corner_lon=0. | units.rad)
    @remote_function(must_handle_array=True)
    def get_dst_grid_corner_lon(i=0):
        returns (corner_lon=0. | units.rad)


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


    #while the grids may have rank > 1, their representation here is flat
    def get_src_grid_firstlast_node(self):
        size = self.parameters.src_grid_size
        return 1,size,1,1
    def get_dst_grid_firstlast_node(self):
        size = self.parameters.dst_grid_size
        return 1,size,1,1

    def get_src_grid_corners_firstlast_node(self):
        size = self.parameters.src_grid_size
        corners = self.parameters.src_grid_corners
        return 1,size*corners,1,1
    def get_dst_grid_corners_firstlast_node(self):
        size = self.parameters.dst_grid_size
        corners = self.parameters.dst_grid_corners
        return 1,size*corners,1,1


    def define_particle_sets(self, object):

        object.define_grid('src_grid')
        object.set_grid_range('src_grid', 'get_src_grid_firstlast_node')
        object.add_setter('src_grid', 'set_src_grid_center_lat', names=('lat',))
        object.add_setter('src_grid', 'set_src_grid_center_lon', names=('lon',))
        object.add_getter('src_grid', 'get_src_grid_center_lat', names=('lat',))
        object.add_getter('src_grid', 'get_src_grid_center_lon', names=('lon',))

        object.define_grid('src_grid_corners')
        object.set_grid_range('src_grid_corners', 'get_src_grid_corners_firstlast_node')
        object.add_setter('src_grid_corners', 'set_src_grid_corner_lat', names=('lat',))
        object.add_setter('src_grid_corners', 'set_src_grid_corner_lon', names=('lon',))
        object.add_getter('src_grid_corners', 'get_src_grid_corner_lat', names=('lat',))
        object.add_getter('src_grid_corners', 'get_src_grid_corner_lon', names=('lon',))

        object.define_grid('dst_grid')
        object.set_grid_range('dst_grid', 'get_dst_grid_firstlast_node')
        object.add_setter('dst_grid', 'set_dst_grid_center_lat', names=('lat',))
        object.add_setter('dst_grid', 'set_dst_grid_center_lon', names=('lon',))
        object.add_getter('dst_grid', 'get_dst_grid_center_lat', names=('lat',))
        object.add_getter('dst_grid', 'get_dst_grid_center_lon', names=('lon',))

        object.define_grid('dst_grid_corners')
        object.set_grid_range('dst_grid_corners', 'get_dst_grid_corners_firstlast_node')
        object.add_setter('dst_grid_corners', 'set_dst_grid_corner_lat', names=('lat',))
        object.add_setter('dst_grid_corners', 'set_dst_grid_corner_lon', names=('lon',))
        object.add_getter('dst_grid_corners', 'get_dst_grid_corner_lat', names=('lat',))
        object.add_getter('dst_grid_corners', 'get_dst_grid_corner_lon', names=('lon',))






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



        object.add_method_parameter(
            "get_src_grid_size",
            "set_src_grid_size",
            "src_grid_size",
            "Set source grid size",
            default_value = 0
        )
        object.add_method_parameter(
            "get_dst_grid_size",
            "set_dst_grid_size",
            "dst_grid_size",
            "Set destination grid size",
            default_value = 0
        )
        object.add_method_parameter(
            "get_src_grid_corners",
            "set_src_grid_corners",
            "src_grid_corners",
            "Set source grid number of corners per grid cell",
            default_value = 0
        )
        object.add_method_parameter(
            "get_dst_grid_corners",
            "set_dst_grid_corners",
            "dst_grid_corners",
            "Set destination grid number of corners per grid cell",
            default_value = 0
        )

        object.add_method_parameter(
            "get_src_grid_dims",
            "set_src_grid_dims",
            "src_grid_dims",
            "Set source grid dimensions, use y=0 for unstructured grids",
            default_value = [0,0]
        )
        object.add_method_parameter(
            "get_dst_grid_dims",
            "set_dst_grid_dims",
            "dst_grid_dims",
            "Set destination grid dimensions, use y=0 for unstructured grids",
            default_value = [0,0]
        )


        object.add_method_parameter(
            "get_src_grid_center_lat",
            "set_src_grid_center_lat",
            "src_grid_center_lat",
            "Set source grid cell center latitudes",
            default_value = 0.0 | units.rad
        )
        object.add_method_parameter(
            "get_dst_grid_center_lat",
            "set_dst_grid_center_lat",
            "dst_grid_center_lat",
            "Set destination grid cell center latitudes",
            default_value = []
        )
        object.add_method_parameter(
            "get_src_grid_center_lon",
            "set_src_grid_center_lon",
            "src_grid_center_lon",
            "Set source grid cell center longitudes",
            default_value = 0.0 | units.rad
        )
        object.add_method_parameter(
            "get_dst_grid_center_lon",
            "set_dst_grid_center_lon",
            "dst_grid_center_lon",
            "Set destination grid cell center longitudes",
            default_value = 0.0 | units.rad
        )

        object.add_method_parameter(
            "get_src_grid_corner_lat",
            "set_src_grid_corner_lat",
            "src_grid_corner_lat",
            "Set source grid cell corner latitudes",
            default_value = 0.0 | units.rad
        )
        object.add_method_parameter(
            "get_dst_grid_corner_lat",
            "set_dst_grid_corner_lat",
            "dst_grid_corner_lat",
            "Set destination grid cell corner latitudes",
            default_value = 0.0 | units.rad
        )
        object.add_method_parameter(
            "get_src_grid_corner_lon",
            "set_src_grid_corner_lon",
            "src_grid_corner_lon",
            "Set source grid cell corner longitudes",
            default_value = 0.0 | units.rad
        )
        object.add_method_parameter(
            "get_dst_grid_corner_lon",
            "set_dst_grid_corner_lon",
            "dst_grid_corner_lon",
            "Set destination grid cell corner longitudes",
            default_value = 0.0 | units.rad
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

        object.add_method('INITIALIZED', 'set_src_grid_center_lon')
        object.add_method('INITIALIZED', 'set_src_grid_center_lat')
        object.add_method('INITIALIZED', 'set_src_grid_corner_lon')
        object.add_method('INITIALIZED', 'set_src_grid_corner_lat')
        object.add_method('INITIALIZED', 'set_dst_grid_center_lon')
        object.add_method('INITIALIZED', 'set_dst_grid_center_lat')
        object.add_method('INITIALIZED', 'set_dst_grid_corner_lon')
        object.add_method('INITIALIZED', 'set_dst_grid_corner_lat')
        object.add_method('INITIALIZED', 'set_src_grid_size')
        object.add_method('INITIALIZED', 'set_src_grid_dims')
        object.add_method('INITIALIZED', 'set_src_grid_corners')
        object.add_method('INITIALIZED', 'set_dst_grid_size')
        object.add_method('INITIALIZED', 'set_dst_grid_dims')
        object.add_method('INITIALIZED', 'set_dst_grid_corners')


        object.add_transition('INITIALIZED', 'RUN', 'commit_parameters')


        for state in ["INITIALIZED","RUN"]:
            object.add_method(state, 'get_weights_file')
            object.add_method(state, 'get_src_grid_file')
            object.add_method(state, 'get_dst_grid_file')


            object.add_method(state, 'get_src_grid_center_lon')
            object.add_method(state, 'get_src_grid_center_lat')
            object.add_method(state, 'get_src_grid_corner_lon')
            object.add_method(state, 'get_src_grid_corner_lat')
            object.add_method(state, 'get_dst_grid_center_lon')
            object.add_method(state, 'get_dst_grid_center_lat')
            object.add_method(state, 'get_dst_grid_corner_lon')
            object.add_method(state, 'get_dst_grid_corner_lat')
            object.add_method(state, 'get_src_grid_size')
            object.add_method(state, 'get_src_grid_dims')
            object.add_method(state, 'get_src_grid_corners')
            object.add_method(state, 'get_dst_grid_size')
            object.add_method(state, 'get_dst_grid_dims')
            object.add_method(state, 'get_dst_grid_corners')



        
        object.add_method("RUN", 'set_src_grid_values')
        object.add_method("RUN", 'get_src_grid_values')
        object.add_method("RUN", 'set_dst_grid_values')
        object.add_method("RUN", 'get_dst_grid_values')
        object.add_method("RUN", 'perform_remap')

        object.add_method("RUN", 'get_src_grid_size')
        object.add_method("RUN", 'get_src_grid_dims')
        object.add_method("RUN", 'get_dst_grid_size')
        object.add_method("RUN", 'get_dst_grid_dims')


