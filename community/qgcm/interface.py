import numpy
from amuse.community.interface.common import CommonCode, CommonCodeInterface 
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system

from fortran_tools import FortranCodeGenerator

from amuse.units import units

import subprocess 

parameters={
# grid parameters
    "number_atmospheric_layers": dict(short="nla", dtype="int32", default="None", description="number of atmospheric layers", ptype="ro"),
    "number_oceanic_layers": dict(short="nlo", dtype="int32", default="None", description="number of oceanic layers", ptype="ro"),
# timestepping
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
    "atmosphere_mixed_bc_coeff": dict(short="bccoat", dtype="float64", default=1. , description="Mixed BC coefficient for atmos.", ptype="simple"),
    "ocean_mixed_bc_coeff": dict(short="bccooc", dtype="float64", default=0.2 , description="Mixed BC coefficient for ocean", ptype="simple"),
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
    "ocean_reduced_gravity" : dict(short="gpoc", dtype="float64", default=None | units.m, ptype="vector", description="reduced gravity for ocean interfaces",length="nlo1"),
# atmospheric QG layer parameters
    "atmosphere_biharmonic_viscosity" : dict(short="ah4at", dtype="float64", default=None | units.m**4/units.s, ptype="vector", description="horizontal biharmonic viscosity coefficient of atm. layers",length="nla"),
    "atmosphere_temperature" : dict(short="tabsat", dtype="float64", default=None | units.K, ptype="vector", description="temperature of atm. layers",length="nla"),
    "atmosphere_layer_thickness" : dict(short="hat", dtype="float64", default=None | units.m, ptype="vector", description="thickness of atm. layers",length="nla"),
    "atmosphere_reduced_gravity" : dict(short="gpat", dtype="float64", default=None | units.m, ptype="vector", description="reduced gravity for atmospheric  interfaces",length="nla1"),
    #~ "Name": dict(short="fname", dtype="float64", default= , description=" ", ptype="simple"),
# other parameters
    "output_directory" : dict(short="outdir", dtype="string", default="outdata", ptype="simple", description="output directory"),
    "ocean_topography_option" : dict(short="topocname", dtype="string", default="flat", ptype="simple", description="ocean topography option: flat, extant(=set in interface) or filename"),
    "atmosphere_topography_option" : dict(short="topatname", dtype="string", default="flat", ptype="simple", description="atm. topography option: flat, extant(=set in interface) or filename"),
}    

code_generator=FortranCodeGenerator(parameters)

class QGCMInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
        .. [#] Hogg et al. 2006, 2003

    """
    
    use_modules = ['qgcm_interface',]

    def name_of_the_worker(self, mode):
        return "q-gcm_worker_"+mode

    def __init__(self, **keyword_arguments):
        mode=keyword_arguments.pop("mode", "ocean_only")
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    exec(code_generator.generate_interface_functions())

    @remote_function
    def commit_grids():
      returns ()

    @remote_function
    def get_model_time():
      returns (model_time = 0. | units.day)

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

    def define_parameters(self, object):
        code_generator.generate_parameter_definitions(object)

    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_transition('INITIALIZED', 'EDIT', 'commit_parameters')
        object.add_transition('EDIT', 'RUN', 'commit_grids')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model')
        #~ object.add_method('INITIALIZED', 'before_get_parameter')
        object.add_method('INITIALIZED', 'before_set_parameter')
        #~ object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED!INITIALIZED!EDIT', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')
        object.add_method('EVOLVED', 'evolve_model')



