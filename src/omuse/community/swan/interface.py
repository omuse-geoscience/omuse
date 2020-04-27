import os.path
from omuse.units import units
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system
from amuse.support.options import option
from amuse import datamodel


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
    "number_of_frequencies" : dict(short="msc", dtype="int32", default=32 , description="number of frequencies", ptype="simple"),
    "number_of_directions" : dict(short="mdc", dtype="int32", default=36 , description="number of directional bins", ptype="simple"),
    "lowest_frequency" : dict(short="slow", dtype="float64", default=0. | units.Hz, description="lowest angular frequency used in freq. discretization", ptype="simple"),
    "highest_frequency" : dict(short="shig", dtype="float64", default=0. | units.Hz , description="highest angular frequency used in freq. discretization", ptype="simple"),
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
    "projection_method" : dict(short="projection_method", dtype="string", default="quasi-cart." , description="projection method (in case of spherical coordinates - doesn't actually do anything atm)", ptype="getter"),
    "number_of_dimensions" : dict(short="number_dimensions", dtype="int32", default=2 , description="number of dimensions (1 (tbd) or 2)", ptype="getter"),
    "use_input_depth" : dict(short="use_input_depth", dtype="bool", default=True , description="use input bathymetry",ptype="simple"),
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
    "timestep" : dict(short="dt", dtype="float64", default=0. | units.s , description="timestep of computation", ptype="simple"),
    "begin_time" : dict(short="begin_time", dtype="float64", default=0. | units.s , description="start time of computation", ptype="simple"),
    "verbosity" : dict(short="itest", dtype="int32", default=1 , description="verbosity of output (0-200)", ptype="simple"),
    "uniform_air_sea_temp_difference" : dict(short="CASTD", dtype="float64", default=0. | units.Celsius , description="uniform air-sea temp. difference", ptype="simple"),
    "wrap_x_coordinate" : dict(short="wrap_x", dtype="bool", default=False , description="whether grid wraps in x-direction", ptype="simple"),
    "number_of_vertices" : dict(short="nvertsg", dtype="int32", default=0, description="number of vertices in case of unstructured grid", ptype="simple"),
    "number_of_cells" : dict(short="ncellsg", dtype="int32", default=0, description="number of cells in case of unstructured grid", ptype="simple"),
    "unstructured_boundary_spec_file" : dict(short="unstructured_boundary_spec_file", dtype="string", default="none" , description="file with wave spectrum on unstructured boundary (1 supported)",ptype="simple"),
    "boundary_marker" : dict(short="boundary_marker", dtype="int32", default=0, description="boundary associated with unstructured spec file", ptype="simple"),
    "use_csigma_cfl_limiter" : dict(short="use_csigma_cfl_limiter", dtype="bool", default=False , description="this option prevents an excessive frequency shifting at a single grid point or vertex due to a very coarse bathymetry or current locally",ptype="simple"),
    "use_ctheta_cfl_limiter" : dict(short="use_ctheta_cfl_limiter", dtype="bool", default=False , description="this option prevents an excessive directional turning at a single grid point or vertex due to a very coarse bathymetry or current locally",ptype="simple"),
    "max_iterations_stationary" : dict(short="mxitst", dtype="int32", default=50, description="maximum number of iterations for stationary calc.", ptype="simple"),
    "max_iterations_dynamic" : dict(short="mxitns", dtype="int32", default=1, description="maximum number of iterations for dynamic calc.", ptype="simple"),
    "air_density" : dict(short="rho_air", dtype="float64", default=1.28 | units.kg/units.m**3, description="density of air", ptype="simple"),
    "minimum_wind_speed" : dict(short="umin", dtype="float64", default=1. | units.m/units.s, description="minimum wind speed", ptype="simple"),
    "stationary_propagation_scheme_order" : dict(short="PROPSS", dtype="int32", default=2, description="order of propagation scheme for stationary calc (1=BSBT, 2= SORDUP)", ptype="simple"),    
    "dynamic_propagation_scheme_order" : dict(short="PROPSN", dtype="int32", default=3, description="order of propagation scheme for stationary calc (1=BSBT, 2= S/L)", ptype="simple"),    
    "turbulent_viscosity_factor" : dict(short="ctb", dtype="float64", default=0.01, description="prefactor in turbulent dissipation formula (0.01)", ptype="simple"),    
    "maximum_error_level" : dict(short="maxerr", dtype="int32", default=1, description="maximum internal error level allowed (1)", ptype="simple"),    
    "under_relaxation_factor" : dict(short="under_relaxation_factor", dtype="float64", default=0., description="stationary under relaxation alpha (0 - 0.01)", ptype="simple"),    
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

def parameter_getter_setters(filename="getter_setters.f90"):
    filestring=""
    py_to_f={"string" : "character(len=*) ", "float64" : "real*8", "float32" : "real", "int32" : "integer", "bool" : "logical"}
    for par,d in parameters.items():
      if d["ptype"] in ["simple"]:
        filestring+=_setter_string.format(d["short"],py_to_f[d["dtype"]])
      if d["ptype"] in ["simple","getter"]:
        filestring+=_getter_string.format(d["short"],py_to_f[d["dtype"]])
    return filestring
    
input_grid_variables={
    "depth" : dict(pyvar=["depth"], forvar=["depth"], igrid=1, unit=units.m),
    "water_level" : dict(pyvar=["water_level"], forvar=["wlevl"], igrid=7, unit=units.m),
    "current" : dict(pyvar=["vx","vy"], forvar=["uxb","uyb"], igrid=2, unit=units.m/units.s),
    "wind" : dict(pyvar=["wind_vx","wind_vy"], forvar=["wxi","wyi"], igrid=5, unit=units.m/units.s),
    "turbulent_visc" : dict(pyvar=["visc"], forvar=["turbf"], igrid=12, unit=units.m**2/units.s),
}

_regular_input_grid_template="""
  function {0}et_input_{1}_regular(i,j,{2},n) result(ret)
    integer :: ret,n,i(n),j(n),k,ii,igrid={3}
    real*8 :: {4}
    ret=0
    do k=1,n
      ii=i(k) + (j(k)-1) * MXG(igrid)
      if(ii.LT.1.OR.ii.GT.MXG(igrid)*MYG(igrid)) THEN
        ret=-1
      else
{5}
      endif
    enddo
  end function
"""

_unstructured_input_grid_template="""
  function {0}et_input_{1}_unstructured(i,{2},n) result(ret)
    integer :: ret,n,i(n),k,ii,igrid={3}
    real*8 :: {4}
    ret=0
    do k=1,n
      ii=i(k)
      if(ii.LT.1.OR.ii.GT.MXG(igrid)*MYG(igrid)) THEN
        ret=-1
      else
{5}
      endif
    enddo
  end function
"""

def input_grid_string(template=_unstructured_input_grid_template):
    filestring=""
    for var,d in input_grid_variables.items():
        for getset,getset_template in zip("sg",["{1}={0}","{0}={1}"]):
            args=d["forvar"]
            igrid=d["igrid"]
            n=len(args)
            forvar=','.join(['x'+str(i) for i,x in enumerate(args)])
            forvarn=','.join(['x'+str(i)+'(n)' for i,x in enumerate(args)])
            getset_lines=[]
            for i,a in enumerate(args):
              getset_lines.append(8*" "+getset_template.format('x'+str(i)+"(k)",a+"(ii)"))
            getset_lines='\n'.join(getset_lines)
            filestring+=template.format(getset,var,forvar,igrid,forvarn,getset_lines)
    return filestring

def generate_getters_setters(filename="getter_setters.f90"):
    filestring=""
    filestring+=input_grid_string(_unstructured_input_grid_template)
    filestring+=input_grid_string(_regular_input_grid_template)
    filestring+=parameter_getter_setters()
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
    def initialize_grid():
        returns ()

    @remote_function
    def initialize_input_grids():
        returns ()

    @remote_function
    def initialize_boundary():
        returns ()


    @remote_function
    def commit_grid_positions():
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
    def get_ac2_regular(i_index='i',j_index='i',k_index='i',l_index='i'):
        returns(ac2='d' | units.m**2*units.s**2/units.rad**2)
    @remote_function(must_handle_array=True)
    def get_depth_regular(i_index="i",j_index="i"):
        returns (depth="d" | units.m)
    @remote_function(must_handle_array=True)
    def get_wave_stress_regular(i_index="i",j_index="i"):
        returns (tau_x="d" | units.Pa,tau_y="d" | units.Pa)

    @remote_function(must_handle_array=True)
    def get_ac2_unstructured(i_index='i',k_index='i',l_index='i'):
        returns(ac2='d' | units.m**2*units.s**2/units.rad**2)
    @remote_function(must_handle_array=True)
    def get_depth_unstructured(i_index="i"):
        returns (depth="d" | units.m)
    @remote_function(must_handle_array=True)
    def get_wave_stress_unstructured(i_index="i"):
        returns (tau_x="d" | units.Pa,tau_y="d" | units.Pa)

    @remote_function(must_handle_array=True)
    def get_grid_position_regular(i_index='i',j_index='i'):
        returns(x='d' | units.m,y='d' | units.m)
    @remote_function(must_handle_array=True)
    def get_grid_lonlat_regular(i_index='i',j_index='i'):
        returns(lon='d' | units.deg,lat='d' | units.deg)

    @remote_function(must_handle_array=True)
    def get_input_grid_position_regular(i_index='i',j_index='i'):
        returns(x='d' | units.m,y='d' | units.m)
    @remote_function(must_handle_array=True)
    def get_input_grid_lonlat_regular(i_index='i',j_index='i'):
        returns(lon='d' | units.deg,lat='d' | units.deg)

    @remote_function(must_handle_array=True)
    def set_grid_position_unstructured(i_index='i',x='d' | units.m,y='d' | units.m):
        returns()
    @remote_function(must_handle_array=True)
    def get_grid_position_unstructured(i_index='i'):
        returns(x='d' | units.m,y='d' | units.m)
    @remote_function(must_handle_array=True)
    def set_grid_lonlat_unstructured(i_index='i',lon='d' | units.deg,lat='d' | units.deg):
        returns()
    @remote_function(must_handle_array=True)
    def get_grid_lonlat_unstructured(i_index='i'):
        returns(lon='d' | units.deg,lat='d' | units.deg)
    @remote_function(must_handle_array=True)
    def set_element_nodes(i_index='i', n1='i', n2='i', n3='i'):
        returns()
    @remote_function(must_handle_array=True)
    def get_element_nodes(i_index='i'):
        returns(n1='i', n2='i', n3='i')
    @remote_function(must_handle_array=True)
    def set_grid_vmark_unstructured(i_index='i',vmark='i'):
        returns()
    @remote_function(must_handle_array=True)
    def get_grid_vmark_unstructured(i_index='i'):
        returns(vmark='i')

    @remote_function
    def get_number_of_unstructured_boundary_segments():
        returns (n_boundaries=0)   
    @remote_function
    def get_number_of_nodes_in_unstructured_boundary_segment(index_of_segment=0):
        returns (n_nodes=0)
    @remote_function(can_handle_array=True)
    def get_unstructured_boundary_node(index=0,index_of_segment=0):
        returns (node_index=0)        

    @remote_function
    def get_exc_value(field_index=0):
        returns (exception_value='d')
    @remote_function
    def set_exc_value(field_index="i",exception_value='d'):
        returns ()

    for par,d in parameters.items():
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

    for var,d in input_grid_variables.items():
        pyvars=d["pyvar"]
        forvars=d["forvar"]
        unit=d["unit"].reference_string()
        for getset in "sg":
            args=[]
            for v in pyvars:
              args.append("  function.addParameter('"+v+"', dtype='d', direction=function."+("IN" if getset=="s" else "OUT")+", unit="+unit+")\n")
            args=''.join(args)
            exec( "@legacy_function\n"
                  "def "+getset+"et_input_"+var+"_regular():\n"
                  "  function = LegacyFunctionSpecification()\n"
                  "  function.must_handle_array = True\n"
                  "  function.addParameter('i_index', dtype='i', direction=function.IN)\n"
                  "  function.addParameter('j_index', dtype='i', direction=function.IN)\n" +
                  args +
                  "  function.addParameter('ncells', dtype='i', direction=function.LENGTH)\n"
                  "  function.result_type = 'i'\n"
                  "  return function\n")
        for getset in "sg":
            args=[]
            for v in pyvars:
              args.append("  function.addParameter('"+v+"', dtype='d', direction=function."+("IN" if getset=="s" else "OUT")+", unit="+unit+")\n")
            args=''.join(args)
            exec( "@legacy_function\n"
                  "def "+getset+"et_input_"+var+"_unstructured():\n"
                  "  function = LegacyFunctionSpecification()\n"
                  "  function.must_handle_array = True\n"
                  "  function.addParameter('i_index', dtype='i', direction=function.IN)\n" +
                  args +
                  "  function.addParameter('ncells', dtype='i', direction=function.LENGTH)\n"
                  "  function.result_type = 'i'\n"
                  "  return function\n") 

        
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

    def initialize_input_grids(self):
        self.overridden().commit_grid_positions()
        self.overridden().initialize_input_grids()

    def commit_parameters(self):
        self.overridden().commit_parameters()
        handler=self.get_handler("PARTICLES")
        self.define_additional_grid_attributes(handler)
        
    def initialize_boundary(self):
        self.overridden().commit_grids()
        self.overridden().initialize_boundary()
                        
    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        #~ object.add_method('INITIALIZED', 'before_get_parameter')
        object.add_method('INITIALIZED', 'before_set_parameter')
        #~ object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        object.add_transition('INITIALIZED','PARAM','commit_parameters')
        object.add_transition('PARAM','GRID','initialize_grid')
        object.add_transition('GRID','INPUTGRID','initialize_input_grids')
        object.add_transition('INPUTGRID','RUN','initialize_boundary')
        object.add_transition('RUN','EVOLVED','evolve_model')

        for param in ["grid_origin_x","grid_origin_y", "grid_orientation",
          "grid_length_x","grid_length_y","grid_nmesh_x","grid_nmesh_y",
          "number_of_frequencies","number_of_directions","lowest_frequency","highest_frequency",
          "input_grid_origin_x","input_grid_origin_y","input_grid_dx",
          "input_grid_dy","input_grid_orientation","input_grid_nmesh_x",
          "input_grid_nmesh_y","number_of_vertices","number_of_cells"]:
            short=parameters[param]['short']
            object.add_method('INITIALIZED', 'set_'+short)

        for state in ['INPUTGRID','RUN','EVOLVED']:
            object.add_method(state, 'before_new_set_instance')
        for state in ['GRID','INPUTGRID','RUN','EVOLVED']:
            object.add_method(state, 'before_new_set_instance2')

        object.add_method('GRID', 'set_grid_position_unstructured')
        object.add_method('GRID', 'set_grid_lonlat_unstructured')
        for var,d in input_grid_variables.items():
            for state in ['INPUTGRID','RUN','EVOLVED']:
                object.add_method(state, 'set_input_'+var+'_regular')
                object.add_method(state, 'set_input_'+var+'_unstructured')
        for state in ['RUN','EVOLVED']:
            object.add_method(state, 'get_input_depth_regular')
            object.add_method(state, 'get_grid_position_unstructured')
            object.add_method(state, 'get_wave_stress_unstructured')
            object.add_method(state, 'get_wave_stress_regular')
            object.add_method(state, 'get_ac2_regular')
        object.add_method('EVOLVED', 'get_depth_regular')
        object.add_method('EVOLVED', 'get_depth_unstructured')
        object.add_method('EVOLVED', 'evolve_model')
    
    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")

    def define_parameters(self, object):      
        for param in parameters:
            short=parameters[param]["short"]
            ptype=parameters[param]["ptype"]
            dtype=parameters[param]["dtype"]
            getter="get_"+short
            if ptype in ["simple","normal"]:
              setter="set_"+short
            else:
              setter=None
            if dtype!='bool':
                object.add_method_parameter(
                    getter,
                    setter,
                    param,
                    parameters[param]["description"], 
                    parameters[param]["default"]
                )
            else:
                object.add_boolean_parameter(
                    getter,
                    setter,
                    param,
                    parameters[param]["description"], 
                    parameters[param]["default"]
                )            
          
    def define_methods(self, object):
        if self._coordinates=="cartesian":
            for p,u in [("grid_xpc",units.m),("grid_ypc",units.m),
                        ("grid_xlenc",units.m),("grid_ylenc",units.m),
                        ("input_xp",units.m),("input_yp",units.m),
                        ("input_dx",units.m),("input_dy",units.m)]:
                object.add_method( 'set_'+p, (u), (object.ERROR_CODE) )
                object.add_method( 'get_'+p, (), (u,object.ERROR_CODE))
        if self._coordinates=="spherical":
            for p,u in [("grid_xpc",units.deg),("grid_ypc",units.deg),
                        ("grid_xlenc",units.deg),("grid_ylenc",units.deg),
                        ("input_xp",units.deg),("input_yp",units.deg),
                        ("input_dx",units.deg),("input_dy",units.deg)]:
                object.add_method( 'set_'+p, (u), (object.ERROR_CODE) )
                object.add_method( 'get_'+p, (), (u,object.ERROR_CODE))

    def get_grid_range(self):
        return 1,self.get_grid_mxc()+1,1,self.get_grid_myc()+1
    def get_grid_range_unstructured(self):
        return 1,self.get_nvertsg()
    def get_element_range_unstructured(self):
        return 1,self.get_ncellsg()
    def get_input_grid_range_regular(self):
        return 1,self.get_input_mx()+1,1,self.get_input_my()+1
    def get_dir_freq_range(self):
        return 1,self.get_mdc(),1,self.get_msc()

    def define_particle_sets(self, object):
        if self._coordinates=="cartesian":
            axes_names=['x','y']
            coordinates="position"
        elif self._coordinates=="spherical":
            axes_names=['lon','lat']
            coordinates="lonlat"

        if self._grid_type=="regular":
            object.define_grid('grid',axes_names = axes_names, state_guard="before_new_set_instance")
            object.set_grid_range('grid', 'get_grid_range')
            object.add_getter('grid', 'get_grid_'+coordinates+'_regular', names=axes_names)
            object.add_getter('grid', 'get_depth_regular', names=["depth"])
#            object.add_getter('grid', 'get_wave_stress_regular',names = ["wave_tau_x","wave_tau_y"])
            object.add_gridded_getter('grid', 'get_ac2_regular','get_dir_freq_range', names = ["ac2"])

        if self._grid_type=="unstructured":
            object.define_grid('nodes',axes_names = axes_names, 
                state_guard="before_new_set_instance2", grid_class=datamodel.UnstructuredGrid)
            object.set_grid_range('nodes', 'get_grid_range_unstructured')
            object.add_getter('nodes', 'get_grid_'+coordinates+'_unstructured', names=axes_names)
            object.add_setter('nodes', 'set_grid_'+coordinates+'_unstructured', names=axes_names)
            object.add_getter('nodes', 'get_grid_vmark_unstructured', names=("vmark",))
            object.add_setter('nodes', 'set_grid_vmark_unstructured', names=("vmark",))
            object.add_gridded_getter('nodes', 'get_ac2_unstructured','get_dir_freq_range', names = ["ac2"])
            object.add_getter('nodes', 'get_depth_unstructured',names = ["depth"])
            object.add_getter('nodes', 'get_wave_stress_unstructured',names = ["wave_tau_x","wave_tau_y"])

            object.define_grid('elements', grid_class=datamodel.UnstructuredGrid)
            object.set_grid_range('elements', 'get_element_range_unstructured')
            object.add_getter('elements', 'get_element_nodes', names=["n1","n2","n3"])
            object.add_setter('elements', 'set_element_nodes', names=["n1","n2","n3"])

        if self._input_grid_type=="regular":
            object.define_grid('forcings',axes_names = axes_names, state_guard="before_new_set_instance")
            object.set_grid_range('forcings', 'get_input_grid_range_regular')

        if self._input_grid_type=="unstructured":
            object.define_grid('forcings',axes_names = axes_names, state_guard="before_new_set_instance",
                                grid_class=datamodel.UnstructuredGrid)
            object.set_grid_range('forcings', 'get_grid_range_unstructured')

    def before_new_set_instance2(self):
        pass

    def define_additional_grid_attributes(self,object):
        if self._coordinates=="cartesian":
            axes_names=['x','y']
            coordinates="position"
        elif self._coordinates=="spherical":
            axes_names=['lon','lat']
            coordinates="lonlat"
            
        if self._input_grid_type=="regular":
            #~ object.define_grid('forcings',axes_names = axes_names)
            #~ object.set_grid_range('forcings', 'get_input_grid_range_regular')
            object.add_getter('forcings', 'get_input_grid_'+coordinates+'_regular', names=axes_names)

        if self._input_grid_type=="unstructured":
            #~ object.define_grid('forcings',axes_names = axes_names)
            #~ object.set_grid_range('forcings', 'get_grid_range_unstructured')
            object.add_getter('forcings', 'get_grid_'+coordinates+'_unstructured', names=axes_names)


        for var,d in input_grid_variables.items():
            #~ if eval("self.get_use_input_"+var+"()"):
            if getattr(self.parameters,"use_input_"+var):
                object.add_getter('forcings', 'get_input_'+var+'_'+self._input_grid_type, names=d["pyvar"])
                object.add_setter('forcings', 'set_input_'+var+'_'+self._input_grid_type, names=d["pyvar"])

