import os.path
from omuse.units import units
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system
from amuse.support.options import option


parameters={
    "constant_water_level": dict(short="wlev", dtype="float64", default=0. | units.m, description="water level that is constant in space and time", ptype="simple"),
    "gravitational_acceleration" : dict(short="grav", dtype="float64", default=9.81 | units.m/units.s**2, description="gravitational acceleration", ptype="simple"),
    "water_density" : dict(short="rho", dtype="float64", default=1025. | units.kg/units.m**3 , description="water density", ptype="simple"),
    "max_wind_drag_coef" : dict(short="cdcap", dtype="float64", default=99999. , description="maximum value for the wind drag coefficient, suggest 0.0025", ptype="simple"),
    "planetary_radius" : dict(short="rearth", dtype="float64", default=6.371e6 | units.m, description="Earth (planetary) radius", ptype="normal"),
    "grid_origin_x" : dict(short="grid_xpc", dtype="float64", default=None , description="x coord. of the origin of the computational grid in the problem coordinate system", ptype="simple"),
    "grid_origin_y" : dict(short="grid_ypc", dtype="float64", default=None , description="y coord. of the origin of the computational grid in the problem coordinate system", ptype="simple"),
    "grid_orientation" : dict(short="grid_alpc", dtype="float64", default=0. | units.deg , description="direction of the positive x-axis of the computational grid (in degrees, Cartesian convention)", ptype="simple"),
    "grid_length_x" : dict(short="grid_xlenc", dtype="float64", default=0. , description="length of the computational grid in x-direction", ptype="simple"),
    "grid_length_y" : dict(short="grid_ylenc", dtype="float64", default=0. , description="length of the computational grid in y-direction", ptype="simple"),
    "grid_nmesh_x" : dict(short="grid_mxc", dtype="int32", default=0 , description="number of meshes in computational grid in x-dir.", ptype="simple"),
    "grid_nmesh_y" : dict(short="grid_myc", dtype="int32", default=0 , description="number of meshes in computational grid in x-dir.", ptype="simple"),
    "numer_of_freq" : dict(short="msc", dtype="int32", default=32 , description="number of frequencies", ptype="simple"),
    "numer_of_directions" : dict(short="mdc", dtype="int32", default=36 , description="number of directional bins", ptype="simple"),
    "lowest_freq" : dict(short="slow", dtype="float64", default=0. | units.Hz, description="lowest angular frequency used in freq. discretization", ptype="simple"),
    "highest_freq" : dict(short="shig", dtype="float64", default=0. | units.Hz , description="highest angular frequency used in freq. discretization", ptype="simple"),
    "input_grid_origin_x" : dict(short="input_xp", dtype="float64", default=0. , description="origin x coord of the input grid", ptype="simple"),
    "input_grid_origin_y" : dict(short="input_yp", dtype="float64", default=0. , description="origin y coord of the input grid", ptype="simple"),
    "input_grid_dx" : dict(short="input_dx", dtype="float64", default=0. , description="input grid x mesh size", ptype="simple"),
    "input_grid_dy" : dict(short="input_dy", dtype="float64", default=0. , description="input grid x mesh size", ptype="simple"),
    "input_grid_orientation" : dict(short="input_alp", dtype="float64", default=0. | units.deg, description="input grid x mesh size", ptype="simple"),
    "input_grid_nmesh_x" : dict(short="input_mx", dtype="int32", default=0, description="input grid x number of mesh cells", ptype="simple"),
    "input_grid_nmesh_y" : dict(short="input_my", dtype="int32", default=0, description="input grid y number of mesh cells", ptype="simple"),
    "uniform_wind_velocity" : dict(short="u10", dtype="float64", default=0. | units.m/units.s , description="wind velocity at 10 m elevation", ptype="simple"),
    "uniform_wind_direction" : dict(short="wdip", dtype="float64", default=0 | units.deg , description="wind direction at 10 m elevation", ptype="normal"),
    "grid_type" : dict(short="grid_type", dtype="string", default="regular" , description="type of grid for computation, [regular, curvilinear or unstructured]", ptype="getter"),
    "input_grid_type" : dict(short="input_grid_type", dtype="string", default="regular" , description="type of grid for input grid, [regular, curvilinear or unstructured]", ptype="getter"),
    "calculation_mode" : dict(short="calc_mode", dtype="string", default="stationary" , description="calculation mode [stationary or dynamic]", ptype="getter"),
    "coordinates" : dict(short="coordinates", dtype="string", default="cartesian" , description="choice of coordinates for input [cartesian, spherical]", ptype="getter"),
    "projection_method" : dict(short="projection_method", dtype="string", default="quasi-cart." , description="projection method (in case of spherical coordinates)", ptype="getter"),
    "number_dimensions" : dict(short="number_dimensions", dtype="int32", default=2 , description="number of dimensions (1 (tbd) or 2)", ptype="getter"),
    "use_uniform_wind" : dict(short="use_uniform_wind", dtype="bool", default=False , description="use constant wind",ptype="simple"),
    "use_input_bottom" : dict(short="use_input_bottom", dtype="bool", default=True , description="use input bathymetry",ptype="simple"),
    "use_input_water_level" : dict(short="use_input_water_level", dtype="bool", default=False , description="use input water level",ptype="simple"),
    "use_input_current" : dict(short="use_input_current", dtype="bool", default=False , description="use input current velocity",ptype="simple"),
    "use_input_air_sea_temp_diff" : dict(short="use_input_air_sea_temp_diff", dtype="bool", default=False , description="use input air sea temp diff.",ptype="simple"),
    "use_input_friction" : dict(short="use_input_friction", dtype="bool", default=False , description="use input bottom friction",ptype="simple"),
    "use_input_wind" : dict(short="use_input_wind", dtype="bool", default=False , description="use input wind velocity",ptype="simple"),
    "use_input_plant_density" : dict(short="use_input_plant_density", dtype="bool", default=False , description="use input plant density",ptype="simple"),
    "use_input_turbulent_visc" : dict(short="use_input_turbulent_visc", dtype="bool", default=False , description="use input turbulent viscosity ",ptype="simple"),
    "use_input_mud_layer" : dict(short="use_input_mud_layer", dtype="bool", default=False , description="use input wave damping mud layrr ",ptype="simple"),
    "use_gen3_parameters" : dict(short="use_gen3", dtype="bool", default=False , description="use third generation wave model parameters set",ptype="simple"),
    "use_friction_parameters" : dict(short="use_friction", dtype="bool", default=False , description="use friction parameter set (change to string parameter)",ptype="simple"),
    "use_breaking_parameters" : dict(short="use_breaking", dtype="bool", default=False , description="use wave breaking",ptype="simple"),
    "use_triads_parameters" : dict(short="use_triads", dtype="bool", default=False , description="use triad interactions",ptype="simple"),
    "north_boundary_spec_file" : dict(short="north_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on north boundary",ptype="simple"),
    "south_boundary_spec_file" : dict(short="south_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on south boundary",ptype="simple"),
    "west_boundary_spec_file" : dict(short="west_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on west boundary",ptype="simple"),
    "east_boundary_spec_file" : dict(short="east_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on east boundary",ptype="simple"),
    "timestep" : dict(short="dt", dtype="float64", default=360. | units.s , description="timestep of computation", ptype="simple"),
    "begin_time" : dict(short="begin_time", dtype="float64", default=0. | units.s , description="start time of computation", ptype="simple"),
    "verbosity" : dict(short="itest", dtype="int32", default=1 , description="verbosity of output (0-200)", ptype="simple"),
    "uniform_air_sea_temp_difference" : dict(short="CASTD", dtype="float64", default=0. | units.Celsius , description="uniform air-sea temp. difference", ptype="simple"),
    "wrap_x_coordinate" : dict(short="wrap_x", dtype="bool", default=False , description="whether grid wraps in x-direction", ptype="simple"),
    "numer_of_vertices" : dict(short="nvertsg", dtype="int32", default=0, description="number of vertices in case of unstructured grid", ptype="simple"),
    "numer_of_cells" : dict(short="ncellsg", dtype="int32", default=0, description="number of cells in case of unstructured grid", ptype="simple"),
#            "parameter_name" : dict(short="abrev.", dtype="float64", default=0 , description=""),
            }

_getter_string="""
  function get_{0}(x) result(ret)
    integer :: ret
    {1} :: x
    x={0}
    ret=0
  end function
               """
_setter_string="""
  function set_{0}(x) result(ret)
  integer :: ret
    {1} :: x
    {0}=x
    ret=0
  end function
               """

def generate_getters_setters(filename="getter_setters.f90"):
    filestring=""
    py_to_f={"string" : "character(len=*) ", "float64" : "real*8", "float32" : "real", "int32" : "integer", "bool" : "logical"}
    for par,d in parameters.iteritems():
      if d["ptype"] in ["simple"]:
        filestring+=_setter_string.format(d["short"],py_to_f[d["dtype"]])
      if d["ptype"] in ["simple","getter"]:
        filestring+=_getter_string.format(d["short"],py_to_f[d["dtype"]])
    with open(filename,"w") as f:
        f.write(filestring)


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

    @remote_function
    def initialize_code(coordinates="cartesian", mode="stationary", grid_type="regular",input_grid_type="regular"):
        returns ()

    @remote_function
    def initialize_grids():
        returns ()

    @remote_function
    def initialize_boundary():
        returns ()


    @remote_function
    def commit_grids():
        returns ()

    @remote_function
    def commit_parameters():
        returns ()

    @remote_function
    def get_time():
        returns (time=0. | units.s)

    @remote_function
    def evolve_model(tend=0. | units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def get_depth_regular(i_index="i",j_index="i"):
        returns (depth="d" | units.m)
    @remote_function(must_handle_array=True)
    def set_depth_regular(i_index="i",j_index="i",depth="d" | units.m):
        returns ()

    @remote_function(must_handle_array=True)
    def get_ac2_regular(i_index='i',j_index='i',k_index='i',l_index='i'):
        returns(ac2='d' | units.m**2*units.s**2/units.rad**2)

    @remote_function(must_handle_array=True)
    def get_grid_position_regular(i_index='i',j_index='i'):
        returns(x='d',y='d')

    #~ @remote_function(must_handle_array=True)
    #~ def get_input_grid_position_regular(i_index='i',j_index='i'):
        #~ returns(x='d',y='d')


    @remote_function
    def get_exc_value(field_index=0):
        returns (exception_value='d')
    @remote_function
    def set_exc_value(field_index="i",exception_value='d'):
        returns ()

    for par,d in parameters.iteritems():
        dtype=d["dtype"]
        if hasattr(d["default"],"unit"):
          unit=d["default"].unit.reference_string()
        else:
          unit="None"
        short=d["short"]
        ptype=d["ptype"]
        exec("@legacy_function\ndef get_"+short+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n"
            "  function.result_type = 'int32'\n  return function")
        if ptype!="getter":
          exec("@legacy_function\ndef set_"+short+"():\n  function = LegacyFunctionSpecification()\n"
              "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.IN, unit="+unit+")\n"
              "  function.result_type = 'int32'\n  return function")    
        
class Swan(InCodeComponentImplementation):
    def __init__(self, coordinates="cartesian", mode="stationary", 
                  grid_type="regular",input_grid_type="regular", **options):
        self._coordinates=coordinates
        self._mode=mode
        self._grid_type=grid_type
        self._input_grid_type=input_grid_type
        InCodeComponentImplementation.__init__(self,  SwanInterface(**options), **options)

    def initialize_code(self):
        self.overridden().initialize_code(self._coordinates,self._mode,
            self._grid_type,self._input_grid_type)

    def commit_grid_and_boundary(self):
        self.overridden().commit_grids()
        self.overridden().initialize_boundary()
        
    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_method('INITIALIZED', 'before_get_parameter')
        object.add_method('INITIALIZED', 'before_set_parameter')
        object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        object.add_transition('INITIALIZED','GRID','initialize_grids')
        object.add_transition('GRID','EDIT','commit_parameters')
        object.add_transition('EDIT','RUN','commit_grid_and_boundary')
        object.add_transition('RUN','EVOLVED','evolve_model')

        for param in ["grid_origin_x","grid_origin_y", "grid_orientation",
          "grid_length_x","grid_length_y","grid_nmesh_x","grid_nmesh_y",
          "numer_of_freq","numer_of_directions","lowest_freq","highest_freq",
          "input_grid_origin_x","input_grid_origin_y","input_grid_dx",
          "input_grid_dy","input_grid_orientation","input_grid_nmesh_x",
          "input_grid_nmesh_y"]:
            short=parameters[param]['short']
            object.add_method('INITIALIZED', 'set_'+short)

        object.add_method('GRID', 'set_depth_regular')
        object.add_method('EDIT', 'set_depth_regular')
        for state in ['EDIT','RUN','EVOLVED']:
            object.add_method(state, 'get_depth_regular')
        object.add_method('RUN', 'get_ac2_regular')
        object.add_method('EVOLVED', 'get_ac2_regular')


    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")

    def define_parameters(self, object):      
        for param in parameters:        
            object.add_default_form_parameter(
                parameters[param]["short"], 
                parameters[param]["description"], 
                parameters[param]["default"] )  
          
    def define_methods(self, object):
        if self._coordinates=="cartesian":
            object.add_method(
                'get_grid_position_regular',
                (object.INDEX,object.INDEX),
                (units.m,units.m,object.ERROR_CODE)
            )
            object.add_method(
                'get_input_grid_position_regular',
                (object.INDEX,object.INDEX),
                (units.m,units.m,object.ERROR_CODE)
            )
            for p,u in [("grid_xpc",units.m),("grid_ypc",units.m),
                        ("grid_xlenc",units.m),("grid_ylenc",units.m),
                        ("input_xp",units.m),("input_yp",units.m),
                        ("input_dx",units.m),("input_dy",units.m)]:
                object.add_method( 'set_'+p, (u), (object.ERROR_CODE) )
                object.add_method( 'get_'+p, (), (u,object.ERROR_CODE))

    def get_grid_range(self):
        return 1,self.get_grid_mxc()+1,1,self.get_grid_myc()+1
    def get_input_grid_range(self):
        return 1,self.get_input_mx()+1,1,self.get_input_my()+1
    def get_dir_freq_range(self):
        return 1,self.get_mdc(),1,self.get_msc()

    def define_particle_sets(self, object):
        if self._coordinates=="cartesian":
            axes_names=['x','y']

        if self._grid_type=="regular":
            object.define_grid('grid',axes_names = axes_names)
            object.set_grid_range('grid', 'get_grid_range')
            object.add_getter('grid', 'get_grid_position_regular', names=axes_names)
            object.add_gridded_getter('grid', 'get_ac2_regular','get_dir_freq_range', names = ["ac2"])

        if self._input_grid_type=="regular":
            object.define_grid('forcings',axes_names = axes_names)
            object.set_grid_range('forcings', 'get_input_grid_range')
            object.add_getter('forcings', 'get_input_grid_position_regular', names=axes_names)
            object.add_getter('forcings', 'get_depth_regular', names=["depth"])
            object.add_setter('forcings', 'set_depth_regular', names=["depth"])
