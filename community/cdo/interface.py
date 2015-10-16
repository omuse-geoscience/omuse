from amuse.community import *

from amuse.community.interface.common import CommonCodeInterface, CommonCode

from amuse.datamodel.grids import *

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

    _src_grid = []
    _dst_grid = []

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

        object.add_method_parameter(
            "get_src_grid",
            "set_src_grid",
            "src_grid",
            "Specify the source grid",
            default_value = []
        )
        object.add_method_parameter(
            "get_dst_grid",
            "set_dst_grid",
            "dst_grid",
            "Specify the destination grid",
            default_value = []
        )


    def set_src_grid(self, grid):
        if not ((type(grid) is UnstructuredGrid) or (type(grid) is StructuredGrid)):
            raise Exception("expected grid type to be either UnstructuredGrid or StructuredGrid, received {0}".format(type(grid).__name__))

        self._src_grid = grid

        if type(grid) is UnstructuredGrid:
            self.set_src_grid_dims(grid.size, 0)
            self.set_src_grid_corners(grid._num_corners)
            self.set_src_grid_center_lon(range(grid.size), grid.lon)
            self.set_src_grid_center_lat(range(grid.size), grid.lat)
            self.set_src_grid_corner_lon(range(grid.size * grid._num_corners), grid._cell_corners[0].flatten())
            self.set_src_grid_corner_lat(range(grid.size * grid._num_corners), grid._cell_corners[1].flatten())

        if type(grid) is StructuredGrid:
            num_corners = 4 #structured grids have exactly 4 corners

            if len(grid.shape) != 2:
                raise Exception("structured grids with not exactly 2 dimensions are not supported, received {0}".format(grid.shape))

            self.set_src_grid_dims(grid.shape)
            self.set_src_grid_corners(num_corners)
            self.set_src_grid_center_lon(range(grid.size), grid.lon.flatten())
            self.set_src_grid_center_lat(range(grid.size), grid.lat.flatten())

            #construct a cell_corners array that CDO understands
            cell_corners = numpy.zeros((2, grid.size*num_corners), dtype=numpy.double)
            for k in range(len(cell_corners)):
                index = 0
                for j in range(grid.shape[0]):
                    for i in range(grid.shape[1]):
                        cell_corners[k,index+0] = grid._cell_corners[k, j  , i  ] #sw
                        cell_corners[k,index+1] = grid._cell_corners[k, j  , i+1] #se
                        cell_corners[k,index+2] = grid._cell_corners[k, j+1, i+1] #ne
                        cell_corners[k,index+3] = grid._cell_corners[k, j+1, i  ] #nw
                        index += 4

            self.set_src_grid_corner_lon(range(grid.size * num_corners), cell_corners[0])
            self.set_src_grid_corner_lat(range(grid.size * num_corners), cell_corners[1])

    def get_src_grid(self):
        return _src_grid



    def set_dst_grid(self, grid):
        if not ((type(grid) is UnstructuredGrid) or (type(grid) is StructuredGrid)):
            raise Exception("expected grid type to be either UnstructuredGrid or StructuredGrid, received {0}".format(type(grid).__name__))

        self._dst_grid = grid

        if type(grid) is UnstructuredGrid:
            self.set_dst_grid_dims(grid.size, 0)
            self.set_dst_grid_corners(grid._num_corners)
            self.set_dst_grid_center_lon(range(grid.size), grid.lon)
            self.set_dst_grid_center_lat(range(grid.size), grid.lat)
            self.set_dst_grid_corner_lon(range(grid.size * grid._num_corners), grid._cell_corners[0].flatten())
            self.set_dst_grid_corner_lat(range(grid.size * grid._num_corners), grid._cell_corners[1].flatten())

        if type(grid) is StructuredGrid:
            num_corners = 4 #structured grids have exactly 4 corners

            if len(grid.shape) != 2:
                raise Exception("structured grids with not exactly 2 dimensions are not supported, received {0}".format(grid.shape))

            self.set_dst_grid_dims(grid.shape)
            self.set_dst_grid_corners(num_corners)
            self.set_dst_grid_center_lon(range(grid.size), grid.lon.flatten())
            self.set_dst_grid_center_lat(range(grid.size), grid.lat.flatten())

            #construct a cell_corners array that CDO understands
            cell_corners = numpy.zeros((2, grid.size*num_corners), dtype=numpy.double)
            for k in range(len(cell_corners)):
                index = 0
                for j in range(grid.shape[0]):
                    for i in range(grid.shape[1]):
                        cell_corners[k,index+0] = grid._cell_corners[k, j  , i  ] #sw
                        cell_corners[k,index+1] = grid._cell_corners[k, j  , i+1] #se
                        cell_corners[k,index+2] = grid._cell_corners[k, j+1, i+1] #ne
                        cell_corners[k,index+3] = grid._cell_corners[k, j+1, i  ] #nw
                        index += 4

            self.set_dst_grid_corner_lon(range(grid.size * num_corners), cell_corners[0])
            self.set_dst_grid_corner_lat(range(grid.size * num_corners), cell_corners[1])

    def get_dst_grid(self):
        return _dst_grid





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

        state = "RUN"
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


