from amuse.rfi.core import CodeInterface
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.parameter_tools import CodeWithIniFileParameters
from amuse.rfi.core import legacy_function,remote_function
from amuse import datamodel

from omuse.units import units

class DFlowFMInterface(CodeInterface):

    use_modules=["dflowfm_omuse"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="dflowfm_worker", **keyword_arguments)
    
    @remote_function
    def initialize():
        returns ()

    @remote_function
    def commit_parameters():
        returns ()
        
    @remote_function
    def evolve_model(tend=0. | units.s):
        returns ()
    
    @remote_function
    def get_model_time():
        returns (model_time = 0. | units.s)
    
    @remote_function
    def get_2d_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_1d_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_2d_boundary_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_1d_boundary_nodes_range():
        returns (imin=0, imax=0)
    
    @remote_function(must_handle_array=True)
    def get_x_position(index=0):
        returns (x=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_y_position(index=0):
        returns (y=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_water_level(index=0):
        returns (h=0. | units.m)

class DFlowFM(InCodeComponentImplementation, CodeWithIniFileParameters):

    def __init__(self, **options):
        self._ini_file=options.get("ini_file","")
        CodeWithIniFileParameters.__init__(self, options.get("parameters", dict()) ) 
        self._coordinates="cartesian"
        InCodeComponentImplementation.__init__(self,  DFlowFMInterface(**options), **options)
        if self._ini_file:
            self.parameters.ini_file=self._ini_file
        
    def define_properties(self, handler):
        handler.add_property('get_model_time', public_name = "model_time")

    def configuration_file_set(self):
        self.read_inifile_parameters(self.parameters.ini_file)
        handler=self.get_handler('PARAMETER')
        CodeWithIniFileParameters.define_parameters(self,handler)


    def define_parameters(self,handler):
        CodeWithIniFileParameters.define_parameters(self, handler)
      
        handler.add_interface_parameter(
            "ini_file",
            "configuration file with simulation setup",
            self._ini_file,
            state_guard="configuration_file_set"
        )

    def define_grids(self, handler):
        if self._coordinates=="cartesian":
            axes_names=['x','y']
            coordinates="position"
        elif self._coordinates=="spherical":
            axes_names=['lon','lat']
            coordinates="lonlat"

        handler.define_grid('flow_2d_nodes',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('flow_2d_nodes', 'get_2d_flow_nodes_range')
        handler.add_getter('flow_2d_nodes', 'get_x_position', names=axes_names[0:1])
        handler.add_getter('flow_2d_nodes', 'get_y_position', names=axes_names[1:2])
        handler.add_getter('flow_2d_nodes', 'get_water_level', names=["water_level"])

        handler.define_grid('boundary_2d_nodes',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('boundary_2d_nodes', 'get_2d_boundary_nodes_range')
        handler.add_getter('boundary_2d_nodes', 'get_x_position', names=axes_names[0:1])
        handler.add_getter('boundary_2d_nodes', 'get_y_position', names=axes_names[1:2])
        handler.add_getter('boundary_2d_nodes', 'get_water_level', names=["water_level"])

    def commit_parameters(self):
        self.write_inifile_parameters("amuse.mdu")
        self.overridden().commit_parameters()

    def define_state(self, handler):
        handler.set_initial_state('UNINITIALIZED')
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize')
        handler.add_transition('INITIALIZED', 'PARAM', 'commit_parameters')

        for method in ["get_2d_flow_nodes_range",
                        "get_1d_flow_nodes_range",
                        "get_1d_boundary_nodes_range",
                        "get_2d_boundary_nodes_range",
                        "evolve_model",
                      ]:
            handler.add_method('PARAM', method)

