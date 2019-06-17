import numpy
from amuse.community.interface.common import CommonCode, CommonCodeInterface 
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system
from amuse.datamodel import CartesianGrid

from fortran_tools import FortranCodeGenerator

from omuse.units import units

import subprocess 

parameters={
    "coriolis_f" : dict(short="fnot", dtype="float64", default=0. | units.rad/units.s, description="coriolis frequency", ptype="ro"),
    "coriolis_beta" : dict(short="beta", dtype="float64", default=0. | units.rad/units.s/units.m, description="coriolis beta parameter", ptype="ro"),
# grid parameters
    "number_atmospheric_layers": dict(short="nla", dtype="int32", default="None", description="number of atmospheric layers", ptype="ro"),
    "number_oceanic_layers": dict(short="nlo", dtype="int32", default="None", description="number of oceanic layers", ptype="ro"),
    "number_ocean_grid_points_x" : dict(short="nxpo", dtype="int32", default="None", description="number of oceanic grid points in x direction", ptype="ro"),
    "number_ocean_grid_points_y" : dict(short="nypo", dtype="int32", default="None", description="number of oceanic grid points in y direction", ptype="ro"),
    "number_atmosphere_grid_points_x" : dict(short="nxpa", dtype="int32", default="None", description="number of atm. grid points in x direction", ptype="ro"),
    "number_atmosphere_grid_points_y" : dict(short="nypa", dtype="int32", default="None", description="number of atm, grid points in y direction", ptype="ro"),
# timestepping
    "begin_time" : dict(short="begin_time", dtype="float64", default=0. | units.s, description="start time of simulation", ptype="simple"),
    "atmosphere_timestep": dict(short="dta", dtype="float64", default=180. | units.s, description="atmospheric timestep", ptype="simple"),
    "timestep_ratio": dict(short="nstr", dtype="int32", default=3, description="timestep ratio oceanic dt/ atmospheric dt", ptype="simple"),
# Physical parameters
    "ocean_grid_spacing": dict(short="dxo", dtype="float64", default=5.e3 | units.m, description="Ocean grid spacing", ptype="simple"),
    "ocean_bottom_ekman_layer_thickness": dict(short="delek", dtype="float64", default=2. | units.m, description="Ocean bottom Ekman layer thickness", ptype="simple"),
    "air_sea_momentum_drag_coefficient": dict(short="cdat", dtype="float64", default=1.3e-3 , description="Air-Sea momentum drag coefficient (quad)", ptype="simple"),
    "atmosphere_density": dict(short="rhoat", dtype="float64", default=1 | units.kg/units.m**3 , description="Atmospheric density ", ptype="simple"),
    "ocean_density": dict(short="rhooc", dtype="float64", default=1000. | units.kg/units.m**3 , description="Ocean density ", ptype="simple"),
    "atmosphere_heat_capacity": dict(short="cpat", dtype="float64", default=1000. | units.J/units.kg/units.K , description="Atmos. specific heat capacity", ptype="simple"),
    "ocean_heat_capacity": dict(short="cpoc", dtype="float64", default=4000. | units.J/units.kg/units.K , description="Ocean specific heat capacity", ptype="simple"),
    "atmosphere_mixed_bc_coeff": dict(short="bccoat", dtype="float64", default=1. , description="Mixed BC coefficient for atmos. (0=freeslip, >2 ~ no-slip)", ptype="simple"),
    "ocean_mixed_bc_coeff": dict(short="bccooc", dtype="float64", default=0.2 , description="Mixed BC coefficient for ocean (0=freeslip, >2 ~ no-slip)", ptype="simple"),
    "coupling_coeff_x": dict(short="xcexp", dtype="float64", default=1. , description="coupling coefficient x", ptype="simple"),
    "coupling_coeff_y": dict(short="ycexp", dtype="float64", default=1. , description="coupling coefficient y", ptype="simple"),
# Mixed layer parameters
    "sensible_latent_transfer": dict(short="xlamda", dtype="float64", default=35. | units.W/units.m**2/units.K , description="Sensible and latent transfer", ptype="simple"),
    "ocean_mixed_layer_depth": dict(short="hmoc", dtype="float64", default=100. |units.m, description="Fixed ocean  mixed layer depth", ptype="simple"),
    "ocean_temperature_2nd_diffusivity": dict(short="st2d", dtype="float64", default=100.| units.m**2/units.s , description="sst grad sqd diffusivity", ptype="simple"),
    "ocean_temperature_4th_diffusivity": dict(short="st4d", dtype="float64", default=2.e9| units.m**4/units.s , description="sst grad 4th diffusivity", ptype="simple"),
    "atmosphere_mixed_layer_depth": dict(short="hmat", dtype="float64", default=1000. | units.m , description="Fixed atmosphere mixed layer depth", ptype="simple"),
    "atmosphere_mixed_layer_depth_min": dict(short="hmamin", dtype="float64", default=100. | units.m , description="minimal allowed  atmosphere mixed layer depth", ptype="simple"),
    "atmosphere_hmix_diffusivity": dict(short="ahmd", dtype="float64", default=2.e5 | units.m**2/units.s, description="atmosphere hmix diffusivity", ptype="simple"),
    "atmosphere_temperature_2nd_diffusivity": dict(short="at2d", dtype="float64", default=2.5e4 | units.m**2/units.s , description="horizontal temperature 2nd order coeff", ptype="simple"),
    "atmosphere_temperature_4th_diffusivity": dict(short="at4d", dtype="float64", default=2.e14 | units.m**4/units.s , description="horizontal temperature 4th order coeff", ptype="simple"),
    "atmosphere_mixed_layer_relaxation_coeff": dict(short="hmadmp", dtype="float64", default=0.15 | units.W/units.m**3 , description="Atmos. m.l. damping constant", ptype="simple"),
#Radiation parameters
    "radiative_forcing": dict(short="fsbar", dtype="float64", default=-210. | units.W/units.m**2 , description="mean radiative forcing", ptype="simple"),
    "radiation_perturbation": dict(short="fspamp", dtype="float64", default=80 | units.W/units.m**2 , description="Radiation perturbation magnitude (W m^-2) (i.e. peak-trough variation, always +ve)", ptype="simple"),
    "optical_depth": dict(short="zm", dtype="float64", default=200. | units.m , description="Optical depth in a.m.l.  (m)", ptype="simple"),
    "adiabatic_lapse_rate": dict(short="gamma", dtype="float64", default=0.01 | units.K /units.m , description="Adiabatic lapse rate", ptype="simple"),
    "optical_depth": dict(short="zopt", dtype="float64", default=None | units.m, ptype="vector", description="Optical depth in layers",length="nla"),
# Oceanic QG layer parameters
    "ocean_viscosity" : dict(short="ah2oc", dtype="float64", default=None | units.m**2/units.s, ptype="vector", description="horizontal viscosity coefficient of ocean layers ",length="nlo"),
    "ocean_biharmonic_viscosity" : dict(short="ah4oc", dtype="float64", default=None | units.m**4/units.s, ptype="vector", description="horizontal biharmonic viscosity coefficient of ocean layers",length="nlo"),
    "ocean_potential_temperature" : dict(short="tabsoc", dtype="float64", default=None | units.K, ptype="vector", description="potential temperature of ocean layers",length="nlo"),
    "ocean_layer_thickness" : dict(short="hoc", dtype="float64", default=None | units.m, ptype="vector", description="thickness of ocean layers",length="nlo"),
    "ocean_reduced_gravity" : dict(short="gpoc", dtype="float64", default=None | units.m/units.s**2, ptype="vector", description="reduced gravity for ocean interfaces",length="nlo1"),
# atmospheric QG layer parameters
    "atmosphere_biharmonic_viscosity" : dict(short="ah4at", dtype="float64", default=None | units.m**4/units.s, ptype="vector", description="horizontal biharmonic viscosity coefficient of atm. layers",length="nla"),
    "atmosphere_temperature" : dict(short="tabsat", dtype="float64", default=None | units.K, ptype="vector", description="temperature of atm. layers",length="nla"),
    "atmosphere_layer_thickness" : dict(short="hat", dtype="float64", default=None | units.m, ptype="vector", description="thickness of atm. layers",length="nla"),
    "atmosphere_reduced_gravity" : dict(short="gpat", dtype="float64", default=None | units.m/units.s**2, ptype="vector", description="reduced gravity for atmospheric  interfaces",length="nla1"),
    #~ "Name": dict(short="fname", dtype="float64", default= , description=" ", ptype="simple"),
# other parameters
    "output_directory" : dict(short="outdir", dtype="string", default="outdata", ptype="simple", description="output directory"),
    "ocean_topography_option" : dict(short="topocname", dtype="string", default="flat", ptype="simple", description="ocean topography option: flat, extant(=set in interface) or filename"),
    "atmosphere_topography_option" : dict(short="topatname", dtype="string", default="flat", ptype="simple", description="atm. topography option: flat, extant(=set in interface) or filename"),
    "use_interface_forcings" : dict( short="use_interface_forcings", dtype="bool", default=False, ptype="simple", description="flag to use interface for forcings, otherwise proper avges.nc file must be present"),
    "check_interface_heights" : dict( short="hcheck", dtype="bool", default=True, ptype="simple", description="whether to check consistency of solution with layer heights"),
    "ocean_ra_alpha" : dict( short="oalpha", dtype="float64", default=0.0, ptype="simple", description="Robert-Asselin filter alpha parameter for ocean (0.01-1) "),
 }    

grid_variables={
    "ocean_dynamic_pressure" : dict(pyvar=["pressure"], forvar=["po"], ndim=3, unit=units.m**2/units.s**2, index_ranges=[(1,"nxpo"),(1,"nypo"),(1,"nlo")]),
    "ocean_dynamic_pressure_tendency" : dict(pyvar=["dpressure_dt"], forvar=["dpo_dt"], ndim=3, unit=units.m**2/units.s**3, index_ranges=[(1,"nxpo"),(1,"nypo"),(1,"nlo")]),
    "ocean_vorticity" : dict(pyvar=["vorticity"], forvar=["qo"], ndim=3, unit=1/units.s, index_ranges=[(1,"nxpo"),(1,"nypo"),(1,"nlo")], vartype="ro"),
    "ocean_topography" : dict(pyvar=["topography"], forvar=["dtopoc"], ndim=2, unit=units.m, index_ranges=[(1,"nxpo"),(1,"nypo")]),
    "ocean_wind_stress" : dict(pyvar=["tau_x","tau_y"], forvar=["tauxo","tauyo"], ndim=2, unit=units.m**2/units.s**2, index_ranges=[(1,"nxpo"),(1,"nypo")]),
    "ocean_surface_heat_flux" : dict(pyvar=["surface_heat_flux"], forvar=["fnetoc"], ndim=2, unit=units.W/units.m**2, index_ranges=[(1,"nxto"),(1,"nyto")]),
    "ocean_surface_temperature_anomaly" : dict(pyvar=["surface_temperature_anomaly"], forvar=["sst"], ndim=2, unit=units.Celsius, index_ranges=[(1,"nxto"),(1,"nyto")]),
    "ocean_surface_temperature_anomaly_tendency" : dict(pyvar=["dsurface_temperature_dt"], forvar=["dsst_dt"], ndim=2, unit=units.Celsius/units.s, index_ranges=[(1,"nxto"),(1,"nyto")]),
    "atmosphere_dynamic_pressure" : dict(pyvar=["pressure"], forvar=["pa"], ndim=3, unit=units.m**2/units.s**2, index_ranges=[(1,"nxpa"),(1,"nypa"),(1,"nla")]),
    "atmosphere_dynamic_pressure_tendency" : dict(pyvar=["dpressure_dt"], forvar=["dpa_dt"], ndim=3, unit=units.m**2/units.s**3, index_ranges=[(1,"nxpa"),(1,"nypa"),(1,"nla")]),
    "atmosphere_vorticity" : dict(pyvar=["vorticity"], forvar=["qa"], ndim=3, unit=1/units.s, index_ranges=[(1,"nxpa"),(1,"nypa"),(1,"nla")], vartype="ro"),
    "atmosphere_surface_temperature_anomaly" : dict(pyvar=["surface_temperature_anomaly"], forvar=["ast"], ndim=2, unit=units.Celsius, index_ranges=[(1,"nxta"),(1,"nyta")]),
    "atmosphere_surface_temperature_anomaly_tendency" : dict(pyvar=["dsurface_temperature_dt"], forvar=["dast_dt"], ndim=2, unit=units.Celsius/units.s, index_ranges=[(1,"nxta"),(1,"nyta")]),
    "atmosphere_mixed_layer_depth" : dict(pyvar=["mixed_layer_depth"], forvar=["hmixa"], ndim=2, unit=units.m, index_ranges=[(1,"nxta"),(1,"nyta")], vartype="ro"),
    "atmosphere_mixed_layer_depth_tendency" : dict(pyvar=["dmixed_layer_depth_dt"], forvar=["dhmixa_dt"], ndim=2, unit=units.m/units.s, index_ranges=[(1,"nxta"),(1,"nyta")], vartype="ro"),
# variables needed for (exact, consistent) restarts
    "ocean_dynamic_pressure_prev" : dict(pyvar=["pressure_prev"], forvar=["pom"], ndim=3, unit=units.m**2/units.s**2, index_ranges=[(1,"nxpo"),(1,"nypo"),(1,"nlo")]),
    "ocean_surface_temperature_anomaly_prev" : dict(pyvar=["surface_temperature_anomaly_prev"], forvar=["sstm"], ndim=2, unit=units.Celsius, index_ranges=[(1,"nxto"),(1,"nyto")]),
    "atmosphere_dynamic_pressure_prev" : dict(pyvar=["pressure_prev"], forvar=["pam"], ndim=3, unit=units.m**2/units.s**2, index_ranges=[(1,"nxpa"),(1,"nypa"),(1,"nla")]),
    "atmosphere_surface_temperature_anomaly_prev" : dict(pyvar=["surface_temperature_anomaly_prev"], forvar=["astm"], ndim=2, unit=units.Celsius, index_ranges=[(1,"nxta"),(1,"nyta")]),
    "atmosphere_mixed_layer_depth_prev" : dict(pyvar=["mixed_layer_depth_prev"], forvar=["hmixam"], ndim=2, unit=units.m, index_ranges=[(1,"nxta"),(1,"nyta")], vartype="ro"),
}



code_generator=FortranCodeGenerator(parameters, grid_variables)

class QGCMInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
        .. [#] Hogg et al. 2006, 2003

    """
    
    use_modules = ['qgcm_interface',]

    def name_of_the_worker(self, mode):
        return "qgcm_worker_"+mode

    def __init__(self, **keyword_arguments):
        mode=keyword_arguments.pop("mode", "ocean_only")
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    exec(code_generator.generate_interface_functions())

    @remote_function
    def commit_grids():
      returns ()

    @remote_function
    def initialize_grids():
      returns ()

    @remote_function
    def get_model_time():
      returns (model_time = 0. | units.day)

    @remote_function
    def get_ocean_time():
      returns (model_time = 0. | units.s)
    @remote_function
    def get_ocean_prev_time():
      returns (model_time = 0. | units.s)

    @remote_function
    def evolve_model(tend = 0. | units.day):
      returns ()

    def get_nlo1(self):
      result=self.get_nlo()
      result["nlo1"]=result["nlo"]-1
      return result

    def get_nla1(self):
      result=self.get_nla()
      result["nla1"]=result["nla"]-1
      return result
      
    
class QGCM(InCodeComponentImplementation):
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  QGCMInterface(**options), **options)

    def get_ocean_P_grid_range(self):
        return 1,self.get_nxpo(),1,self.get_nypo()
    def get_ocean_T_grid_range(self):
        return 1,self.get_nxpo()-1,1,self.get_nypo()-1
    def get_atmosphere_T_grid_range(self):
        return 1,self.get_nxpa()-1,1,self.get_nypa()-1

    def get_ocean_P_grid_position(self,i,j):
        dxo=self.get_dxo()
        return dxo*(i-1),dxo*(j-1)
    def get_ocean_T_grid_position(self,i,j):
        dxo=self.get_dxo()
        return dxo*(i-0.5),dxo*(j-0.5)
# fix offset!
    def get_atmosphere_T_grid_position(self,i,j):
        dxa=self.get_dxa()
        return dxa*(i-0.5),dxa*(j-0.5)

    def get_nlo_range(self):
        return 1, self.get_nlo()

    def define_parameters(self, object):
        code_generator.generate_parameter_definitions(object)

    def define_properties(self, object):
        object.add_property('get_model_time', public_name = "model_time")

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_transition('INITIALIZED', 'PARAM', 'commit_parameters')
        object.add_transition('PARAM', 'EDIT', 'initialize_grids')
        object.add_transition('EDIT', 'RUN', 'commit_grids')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model')
        object.add_method('!UNINITIALIZED', 'before_get_parameter')
        object.add_method('INITIALIZED', 'before_set_parameter')
        #~ object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED!INITIALIZED!PARAM!EDIT', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')
        object.add_method('EVOLVED', 'evolve_model')

        for state in ["EDIT",'RUN','EVOLVED']:
            object.add_method(state, 'before_new_set_instance')
        for state in ['RUN','EVOLVED']:
            object.add_method(state, "get_tauxo")
            object.add_method(state, "get_tauyo")

        object.add_method("EDIT","set_po")
        object.add_method("EDIT","set_dpo_dt")
        # if grids var po, dpo_dt etc are touched then go back to state EDIT


    def define_grids(self, object):
        #~ code_generator.generate_grid_definitions(object)
        
        # ocean dynamic variables p-grid: po, dpo_dt, vorticity
        object.define_grid('ocean_P_grid',axes_names = "xy", grid_class=CartesianGrid,state_guard="before_new_set_instance")
        object.set_grid_range('ocean_P_grid', 'get_ocean_P_grid_range')
        object.add_getter('ocean_P_grid', 'get_ocean_P_grid_position', names="xy")
        object.add_gridded_getter('ocean_P_grid', 'get_po', "get_nlo_range", names=["pressure"])
        object.add_gridded_setter('ocean_P_grid', 'set_po', "get_nlo_range", names=["pressure"])
        object.add_gridded_getter('ocean_P_grid', 'get_dpo_dt', "get_nlo_range", names=["dpressure_dt"])
        object.add_gridded_setter('ocean_P_grid', 'set_dpo_dt', "get_nlo_range", names=["dpressure_dt"])
        object.add_gridded_getter('ocean_P_grid', 'get_qo', "get_nlo_range", names=["vorticity"])
        #~ object.add_getter('ocean_P_grid', 'get_dtopoc', names=["topography"])

        object.define_grid('ocean_T_grid',axes_names = "xy", grid_class=CartesianGrid,state_guard="before_new_set_instance")
        object.set_grid_range('ocean_T_grid', 'get_ocean_T_grid_range')
        object.add_getter('ocean_T_grid', 'get_ocean_T_grid_position', names="xy")
        object.add_getter('ocean_T_grid', 'get_sst', names=["surface_temperature_anomaly"])
        object.add_setter('ocean_T_grid', 'set_sst', names=["surface_temperature_anomaly"])
        object.add_getter('ocean_T_grid', 'get_dsst_dt', names=["dsurface_temperature_anomaly_dt"])
        object.add_setter('ocean_T_grid', 'set_dsst_dt', names=["dsurface_temperature_anomaly_dt"])


        
        object.define_grid('ocean_P_grid_forcings',axes_names = "xy", grid_class=CartesianGrid,state_guard="before_new_set_instance")
        object.set_grid_range('ocean_P_grid_forcings', 'get_ocean_P_grid_range')
        object.add_getter('ocean_P_grid_forcings', 'get_ocean_P_grid_position', names="xy")

        object.add_getter('ocean_P_grid_forcings', 'get_tauxo', names=["tau_x"])
        object.add_getter('ocean_P_grid_forcings', 'get_tauyo', names=["tau_y"])
        object.add_setter('ocean_P_grid_forcings', 'set_tauxo', names=["tau_x"])
        object.add_setter('ocean_P_grid_forcings', 'set_tauyo', names=["tau_y"])
        object.add_getter('ocean_P_grid_forcings', 'get_dtopoc', names=["topography"])
        object.add_setter('ocean_P_grid_forcings', 'set_dtopoc', names=["topography"])
        
        object.define_grid('ocean_T_grid_forcings',axes_names = "xy", grid_class=CartesianGrid,state_guard="before_new_set_instance")
        object.set_grid_range('ocean_T_grid_forcings', 'get_ocean_T_grid_range')
        object.add_getter('ocean_T_grid_forcings', 'get_ocean_T_grid_position', names="xy")

        object.add_getter('ocean_T_grid_forcings', 'get_fnetoc', names=["surface_heat_flux"])
        object.add_setter('ocean_T_grid_forcings', 'set_fnetoc', names=["surface_heat_flux"])

        object.define_grid('atmosphere_T_grid',axes_names = "xy", grid_class=CartesianGrid,state_guard="before_new_set_instance")
        object.set_grid_range('atmosphere_T_grid', 'get_atmosphere_T_grid_range')
        object.add_getter('atmosphere_T_grid', 'get_atmosphere_T_grid_position', names="xy")

        object.add_getter('atmosphere_T_grid', 'get_hmixa', names=["mixed_layer_depth"])
        object.add_getter('atmosphere_T_grid', 'get_dhmixa_dt', names=["dmixed_layer_depth_dt"])
        #~ object.add_setter('atmosphere_T_grid', 'set_hmixa', names=["mixed_layer_depth"])
        #~ object.add_setter('atmosphere_T_grid', 'set_dhmixa_dt', names=["dmixed_layer_depth_dt"])

        
        
    def commit_parameters(self):
        self.overridden().commit_parameters()
        handler=self.get_handler("DATASETS")
        self.define_additional_grid_attributes(handler)

    def define_additional_grid_attributes(self,object):
        pass

