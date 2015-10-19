import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.units.core import system,no_system
from amuse.support.options import option

parameters={
            "constant_water_level": dict(short="wlev", dtype="float64", default=0. | units.m ,description="water level that is constant in space and time", simple=True),
            "gravitational_acceleration" : dict(short="grav", dtype="float64", default=9.81 | units.m/units.s**2, description="gravitational acceleration", simple=True),
            "water_density" : dict(short="rho", dtype="float64", default=1025. | units.kg/units.m**3 , description="water density", simple=True),
            "max_wind_drag_coef" : dict(short="cdcap", dtype="float64", default=99999. , description="maximum value for the wind drag coefficient, suggest 0.0025", simple=True),
            "planetary_radius" : dict(short="rearth", dtype="float64", default=6.371e6 | units.m, description="Earth (planetary) radius", simple=False),
            "grid_origin_x" : dict(short="grid_xpc", dtype="float64", default=0. , description="x coord. of the origin of the computational grid in the problem coordinate system", simple=True),
            "grid_origin_y" : dict(short="grid_ypc", dtype="float64", default=0. , description="y coord. of the origin of the computational grid in the problem coordinate system", simple=True),
            "grid_orientation" : dict(short="grid_alpc", dtype="float64", default=0. | units.deg , description="direction of the positive x-axis of the computational grid (in degrees, Cartesian convention)", simple=True),
            "grid_length_x" : dict(short="grid_xlenc", dtype="float64", default=0. , description="length of the computational grid in x-direction", simple=True),
            "grid_length_y" : dict(short="grid_ylenc", dtype="float64", default=0. , description="length of the computational grid in y-direction", simple=True),
            "grid_nmesh_x" : dict(short="grid_mxc", dtype="int32", default=0 , description="number of meshes in computational grid in x-dir.", simple=True),
            "grid_nmesh_y" : dict(short="grid_myc", dtype="int32", default=0 , description="number of meshes in computational grid in x-dir.", simple=True),
            "numer_of_freq" : dict(short="nfreq", dtype="int32", default=32 , description="number of frequencies", simple=False),
            "numer_of_directions" : dict(short="mdc", dtype="int32", default=36 , description="number of directional bins", simple=True),
            "lowest_freq" : dict(short="flow", dtype="float64", default=0. | units.Hz, description="lowest frequency used in freq. discretization", simple=False),
            "highest_freq" : dict(short="fhigh", dtype="float64", default=0. | units.Hz , description="highest frequency used in freq. discretization", simple=False),
            "input_grid_origin_x" : dict(short="input_xp", dtype="float64", default=0. , description="origin x coord of the input grid", simple=True),
            "input_grid_origin_y" : dict(short="input_yp", dtype="float64", default=0. , description="origin y coord of the input grid", simple=True),
            "input_grid_mesh_dx" : dict(short="input_dx", dtype="float64", default=0. , description="input grid x mesh size", simple=True),
            "input_grid_mesh_dy" : dict(short="input_dy", dtype="float64", default=0. , description="input grid x mesh size", simple=True),
            "input_grid_mesh_orientation" : dict(short="input_alp", dtype="float64", default=0. | units.deg, description="input grid x mesh size", simple=True),
            "input_grid_nmesh_x" : dict(short="input_mx", dtype="int32", default=0, description="input grid x number of mesh cells", simple=True),
            "input_grid_nmesh_y" : dict(short="input_my", dtype="int32", default=0, description="input grid y number of mesh cells", simple=True),
            "uniform_wind_vel" : dict(short="uniform_wind_vel", dtype="float64", default=0. | units.m/units.s , description="wind velocity at 10 m elevation", simple=False),
            "uniform_wind_direction" : dict(short="uniform_wind_dir", dtype="float64", default=0 | units.deg , description="wind direction at 10 m elevation", simple=False),
            "grid_type" : dict(short="grid_type", dtype="string", default="regular" , description="type of grid for computation, [regular, curvilinear or unstructured]", simple=True),
            "input_grid_type" : dict(short="input_grid_type", dtype="string", default="regular" , description="type of grid for input grid, [regular, curvilinear or unstructured]", simple=True),
            "calculation_mode" : dict(short="calc_mode", dtype="string", default="stationary" , description="calculation mode [stationary or dynamic]", simple=True),
            "coordinates" : dict(short="coordinates", dtype="string", default="cartesian" , description="choice of coordinates for input [cartesian, spherical]", simple=True),
            "projection_method" : dict(short="projection_method", dtype="string", default="quasi-cart." , description="projection method (in case of spherical coordinates)", simple=True),
            "number_dimensions" : dict(short="number_dimensions", dtype="int32", default=2 , description="number of dimensions (1 (tbd) or 2)", simple=True),
            "use_uniform_wind" : dict(short="use_uniform_wind", dtype="bool", default=False , description="use constant wind",simple=True),
            "use_input_bottom" : dict(short="use_input_bottom", dtype="bool", default=True , description="use input bathymetry",simple=True),
            "use_input_water_level" : dict(short="use_input_water_level", dtype="bool", default=False , description="use input water level",simple=True),
            "use_input_current" : dict(short="use_input_current", dtype="bool", default=False , description="use input current velocity",simple=True),
            "use_input_air_sea_temp_diff" : dict(short="use_input_air_sea_temp_diff", dtype="bool", default=False , description="use input air sea temp diff.",simple=True),
            "use_input_friction" : dict(short="use_input_friction", dtype="bool", default=False , description="use input bottom friction",simple=True),
            "use_input_wind" : dict(short="use_input_wind", dtype="bool", default=False , description="use input wind velocity",simple=True),
            "use_input_plant_density" : dict(short="use_input_plant_density", dtype="bool", default=False , description="use input plant density",simple=True),
            "use_input_turbulent_visc" : dict(short="use_input_turbulent_visc", dtype="bool", default=False , description="use input turbulent viscosity ",simple=True),
            "use_input_mud_layer" : dict(short="use_input_mud_layer", dtype="bool", default=False , description="use input wave damping mud layrr ",simple=True),
            "use_gen3_parameters" : dict(short="use_gen3", dtype="bool", default=False , description="use third generation wave model parameters set",simple=True),
            "use_friction_parameters" : dict(short="use_friction", dtype="bool", default=False , description="use friction parameter set (change to string parameter)",simple=True),
            "use_breaking_parameters" : dict(short="use_breaking", dtype="bool", default=False , description="use wave breaking",simple=True),
            "use_triads_parameters" : dict(short="use_triads", dtype="bool", default=False , description="use triad interactions",simple=True),
            "north_boundary_spec_file" : dict(short="north_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on north boundary",simple=True),
            "south_boundary_spec_file" : dict(short="south_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on south boundary",simple=True),
            "west_boundary_spec_file" : dict(short="west_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on west boundary",simple=True),
            "east_boundary_spec_file" : dict(short="east_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on east boundary",simple=True),
            "timestep" : dict(short="dt", dtype="float64", default=360. | units.s , description="timestep of computation", simple=True),
            "begin_time" : dict(short="begin_time", dtype="float64", default=0. | units.s , description="start time of computation", simple=True),
#            "parameter_name" : dict(short="abrev.", dtype="float64", default=0 , description=""),
            }

def generate_getters_setters(filename="getter_setters.f90"):
    filestring=""
    py_to_f={"string" : "character(len=*) ", "float64" : "real*8", "int32" : "integer", "bool" : "logical"}
    for par,d in parameters.iteritems():
      if d["simple"]:
        filestring+="""
                    function set_{0}(x) result(ret)
                    integer :: ret
                      {1} :: x
                      {0}=x
                      ret=0
                    end function
                    function get_{0}(x) result(ret)
                      integer :: ret
                      {1} :: x
                      x={0}
                      ret=0
                    end function
                    """.format(d["short"],py_to_f[d["dtype"]])
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
    def get_depth(i_index="i",j_index="i"):
        returns (depth="f" | units.m)
    @remote_function(must_handle_array=True)
    def set_depth(i_index="i",j_index="i",depth="f" | units.m):
        returns ()

    @remote_function
    def get_exc_value(field_index=0):
        returns (exception_value='f')

    for par,d in parameters.iteritems():
        dtype=d["dtype"]
        if hasattr(d["default"],"unit"):
          unit=d["default"].unit.reference_string()
        else:
          unit="None"
        short=d["short"]
        exec("@legacy_function\ndef get_"+short+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n"
            "  function.result_type = 'int32'\n  return function")
        exec("@legacy_function\ndef set_"+short+"():\n  function = LegacyFunctionSpecification()\n"
            "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.IN, unit="+unit+")\n"
            "  function.result_type = 'int32'\n  return function")    
        
