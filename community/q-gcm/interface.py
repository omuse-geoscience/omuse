import numpy
from amuse.community.interface.common import CommonCode, CommonCodeInterface 
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system


from amuse.units import units

import subprocess 

parameters={
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
    #~ "Name": dict(short="fname", dtype="float64", default= , description=" ", ptype="simple"),
}    
"""  
 2.0d4  2.0d4  3.0d4        !! zopt(k) = Optical depth in layer k (m)
!!
!! Oceanic QG layer parameters
!! ---------------------------
  0.0d0    0.0d0    0.0d0   !! ah2oc(k)  = Del-sqd coefft for ocean layer k (m^2 s^-1)
!5.0d+9   5.0d+9   5.0d+9   !! ah4oc(k)  = Del-4th coefft for ocean layer k (m^4 s^-1)
 2.0d+9   2.0d+9   2.0d+9   !! ah4oc(k)  = Del-4th coefft for ocean layer k (m^4 s^-1)
 287.0d0  282.0d0  276.0d0  !! tabsoc(k) = Potential temp. of ocean layer k (K)
 350.0d0  750.0d0 2900.0d0  !! hoc(k)    = Thickness of ocean layer k (m)
!! 0.0500d0  0.0250d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!! 0.0400d0  0.0200d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!! 0.0350d0  0.0175d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!! 0.0300d0  0.0150d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!! 0.0250d0  0.0125d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!! 0.0200d0  0.0100d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
 0.0150d0  0.0075d0         !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!!
!! Atmospheric QG layer parameters
!! -------------------------------
 1.5d+14  1.5d+14  1.5d+14  !! ah4at(k)  = Del-4th coefft for atmos. layer k (m^4 s^-1)
 330.0d0  340.0d0  350.0d0  !! tabsat(k) = Temperature for atmos. layer k (K)
 2000.d0  3000.d0  4000.d0  !! hat(k)    = Thickness of atmos. layer k (m)
 1.2d0     0.4d0            !! gpat(i)   = Reduced gravity for atmos. interface i (m s^-2)
"""

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
    for par,d in parameters.iteritems():
      if d["ptype"] in ["simple"]:
        filestring+=_setter_string.format(d["short"],py_to_f[d["dtype"]])
      if d["ptype"] in ["simple","getter"]:
        filestring+=_getter_string.format(d["short"],py_to_f[d["dtype"]])
    return filestring

def generate_getters_setters(filename="getter_setters.f90"):
    filestring=""
    #~ filestring+=input_grid_string(_unstructured_input_grid_template)
    #~ filestring+=input_grid_string(_regular_input_grid_template)
    filestring+=parameter_getter_setters()
    with open(filename,"w") as f:
        f.write(filestring)

class QGCMInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
        .. [#] Hogg et al. 2006, 2003

    """
    
    use_modules = ['qgcm_interface',]

    def name_of_the_worker(self, mode):
        return "q-gcm_worker_"+mode

    def __init__(self, **keyword_arguments):
        mode=keyword_arguments["mode"]
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(mode), **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

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


class QGCM(InCodeComponentImplementation):
    def __init__(self, mode="ocean_only", **options):
        InCodeComponentImplementation.__init__(self,  QGCMInterface(mode=mode,**options), **options)

