from omuse.units import units

# dict of parameters: keys become the OMUSE parameter names, short is the name in the namelist file
namelist_parameters=dict(
#
itot  =  dict(group_name="DOMAIN", short="itot", dtype="int32", default=64, description="Number of horizontal grid points in x-direction", ptype="nml"),
jtot  =  dict(group_name="DOMAIN", short="jtot", dtype="int32", default=64, description="Number of horizontal grid points in y-direction", ptype="nml"),
kmax  =  dict(group_name="DOMAIN", short="kmax", dtype="int32", default=96, description="Number of vertical grid points", ptype="nml"),
xsize =  dict(group_name="DOMAIN", short="xsize", dtype="float64", default=6400. | units.m, description="Horizontal size of the simulated domain in x", ptype="nml"),
ysize =  dict(group_name="DOMAIN", short="ysize", dtype="float64", default=6400. | units.m, description="Horizontal size of the simulated domain in y", ptype="nml"),
xlat  =  dict(group_name="DOMAIN", short="xlat", dtype="float64", default=52 | units.deg, description="latitude of domain", ptype="nml"),
xlon  =  dict(group_name="DOMAIN", short="xlon", dtype="float64", default=0 | units.deg, description="longitude of domain", ptype="nml"),
xyear =  dict(group_name="DOMAIN", short="xyear", dtype="float64", default=0 | units.yr, description="year, only for time units in netcdf", ptype="nml"),
xday  =  dict(group_name="DOMAIN", short="xday", dtype="float64", default=1 | units.day, description="start day, number of the day", ptype="nml"),
xtime =  dict(group_name="DOMAIN", short="xtime", dtype="float64", default=1 | units.hour, description="start time, UTC time of the day", ptype="nml"),
ksp   =  dict(group_name="DOMAIN", short="ksp", dtype="float64", default=None, description="lower height of sponge layer in grid cells", ptype="nml"),
#
llsadv     =  dict(group_name="DYNAMICS", short="llsadv", dtype="bool", default=False, description="Switch for large scale forcings", ptype="nml"),
lqlnr      =  dict(group_name="DYNAMICS", short="lqlnr", dtype="bool", default=True, description="Switch for Newton-Raphson approximation of the liquid water content", ptype="nml"),
lnoclouds  =  dict(group_name="DYNAMICS", short="lnoclouds", dtype="bool", default=False, description="Switch to disable q_l calculations", ptype="nml"),
cu         =  dict(group_name="DYNAMICS", short="cu", dtype="float64", default=0. | units.m/units.s, description="Transformation velocity of the Galilei transformation in x-direction",ptype="nml"),
cv         =  dict(group_name="DYNAMICS", short="cv", dtype="float64", default=0. | units.m/units.s, description="Transformation velocity of the Galilei transformation in y-direction",ptype="nml"),
iadv_mom   =  dict(group_name="DYNAMICS", short="iadv_mom", dtype="int32", default=5, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
iadv_tke   =  dict(group_name="DYNAMICS", short="iadv_tke", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
iadv_thl   =  dict(group_name="DYNAMICS", short="iadv_thl", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
iadv_qt    =  dict(group_name="DYNAMICS", short="iadv_qt", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
iadv_sv    =  dict(group_name="DYNAMICS", short="iadv_sv", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
ibas_prf   =  dict(group_name="DYNAMICS", short="ibas_prf", dtype="int32", default=3, description="Flag for density calculations, -1 = copy iadv_mom, 1= based constant theta_0, 2 = Boussinesq-like, 3=Standard lapse rate, based on surface temp, 4= standard lr, tsurf=15C, 5=user defined",ptype="nml"),
lambda_crit=  dict(group_name="DYNAMICS", short="lambda_crit", dtype="float64", default=100, description="Maximum value for the smoothness.",ptype="nml"),
#
z0         =  dict(group_name="PHYSICS", short="z0", dtype="float64", default=1.6e-4 | units.m, description ="surface roughness (negative = ", ptype="nml"),
ustin      =  dict(group_name="PHYSICS", short="ustin", dtype="float64", default=0.32 | units.m/units.s, description ="Prescribed friction velocity (negative = ", ptype="nml"),
wtsurf     =  dict(group_name="PHYSICS", short="wtsurf", dtype="float64", default=0. | units.K*units.m/units.s, description ="Flux of liq. water pot. temp. at the surface (negative = ", ptype="nml"),
wqsurf     =  dict(group_name="PHYSICS", short="wqsurf", dtype="float64", default=0. | units.K*units.m/units.s, description ="Flux of total water content at the surface (negative = ", ptype="nml"),
wsvsurf    =  dict(group_name="PHYSICS", short="wsvsurf", dtype="float64", default=0. | units.ppb*units.m/units.s, description ="Flux of scalars at the surface", ptype="nml"),
thls       =  dict(group_name="PHYSICS",short="thls",dtype="float64",default=298.5 | units.K, description="Liquid water potential temperature at the surface", ptype="nml"),
ps         =  dict(group_name="PHYSICS",short="ps",dtype="float64",default=101540. | units.Pa, description="surface pressure", ptype="nml"),
isurf      =  dict(group_name="PHYSICS",short="isurf",dtype="int32",default=4, description="surface parametrization (1 [interactive, use rad.],2 [forced surf T, flux calculated],3 [forced mom, moisture and heat flux],4 [forced moisture and heat flux] or 10 [user def]", ptype="nml"),
ltimedep   =  dict(group_name="PHYSICS",short="ltimedep", dtype="bool", default=False, description="switch for timedependent fluxes", ptype="nml"),
lcoriol    =  dict(group_name="PHYSICS",short="lcoriol", dtype="bool", default=True, description="switch for coriolis force", ptype="nml"),
igrw_damp  =  dict(group_name="PHYSICS",short="igrw_damp", dtype="int32", default=2, description="flag for grav wave damping ( -1 [nudge to lscale.inp], 0 [no damping], 1 [fast damping to av. and slow damping of av to geowind], 2[fast damping to geowind] or 3[fast damping to av wind]", ptype="nml"), 
geodamptime=  dict(group_name="PHYSICS",short="geodamptime", dtype="float64", default=7200. | units.s, description="timescale for nudging to geowind", ptype="nml"),
lmomsubs   =  dict(group_name="PHYSICS",short="lmomsubs", dtype="bool", default=False, description="switch to apply subsidence on momentum", ptype="nml"),
lmoist     =  dict(group_name="PHYSICS",short="lmoist", dtype="bool", default=True, description="switch for calculation of moisture field", ptype="nml"),
chi_half   =  dict(group_name="PHYSICS",short="chi_half", dtype="float64", default=0.5, description="Wet, dry or intermediate (default) mixing over the cloud edge", ptype="nml"),
timerad   =   dict(group_name="PHYSICS",short="timerad", dtype="float64", default=0. | units.s, description="Sampling interval of radiation scheme", ptype="nml"),
iradiation =  dict(group_name="PHYSICS",short="iradiation", dtype="int32", default=0, description="flag for rad. calc. (0 [no rad], 1 [full], 2 [parametrized], 3 [simple for land surface], 4 [rapid RTM] or 10 [user]", ptype="nml"),
useMcICA   =  dict(group_name="PHYSICS",short="useMcICA", dtype="bool", default=True, description="switch for the Monte Carlo Indep. Column approach (only for iradiation=1)", ptype="nml"), 
rad_ls     =  dict(group_name="PHYSICS",short="rad_ls", dtype="bool", default=True, description="switch for prescribed radiation forcing", ptype="nml"), 
rad_longw  =  dict(group_name="PHYSICS",short="rad_longw", dtype="bool", default=True, description="switch for parameterized long wave radiation forcing", ptype="nml"), 
rad_shortw =  dict(group_name="PHYSICS",short="rad_shortw", dtype="bool", default=True, description="switch for parameterized short wave radiation forcing", ptype="nml"), 
rad_smoke  =  dict(group_name="PHYSICS",short="rad_smoke", dtype="bool", default=False, description="switch for longwave divergence for smoke cloud", ptype="nml"), 
rka        =  dict(group_name="PHYSICS",short="rka",dtype="float64", default=130. | units.m**2/units.kg, description="Extinction coefficient (if iradiation=2)", ptype="nml"),
dlwbot     =  dict(group_name="PHYSICS",short="dlwbot",dtype="float64", default=0. | units.W/units.m**2, description="long wave rad flux jump at cloud bottom", ptype="nml"),
dlwtop     =  dict(group_name="PHYSICS",short="dlwtop",dtype="float64", default=74. | units.W/units.m**2, description="long wave rad flux jump at cloud top", ptype="nml"),
sw0        =  dict(group_name="PHYSICS",short="sw0",dtype="float64", default=1100. | units.W/units.m**2, description="Direct solar radiative component at cloud top (assumes zero diffusive contr.)", ptype="nml"),
gc         =  dict(group_name="PHYSICS",short="gc",dtype="float64", default=0.85, description="Asymmetry factor for droplet scattering", ptype="nml"),
reff       =  dict(group_name="PHYSICS",short="reff",dtype="float64", default=1.e-5 | units.m, description="cloud drop effective radius", ptype="nml"),
isvsmoke   =  dict(group_name="PHYSICS",short="isvsmoke",dtype="int32", default=1, description="scalar field to be used for optical depth calc. (when rad_smoke=True)", ptype="nml"),
# omitted lforce_user
lcloudshading=  dict(group_name="PHYSICS",short="lcloudshading", dtype="bool", default=False, description="switch to let clouds shade the surface for rad_lsm", ptype="nml"),
lrigidlid  =  dict(group_name="PHYSICS",short="lrigidlid", dtype="bool", default=False, description="Switch to enable simulations with a rigid lid", ptype="nml"),
unudge     =  dict(group_name="PHYSICS",short="unudge", dtype="float64", default=1., description="Nudging factor if igrw_damp is -1", ptype="nml"),
#
iexpnr     =  dict(group_name="RUN",short="iexpnr", dtype="int32", default=0, description="Experiment number", ptype="nml"),
dtmax      =  dict(group_name="RUN",short="dtmax", dtype="float64", default=20. | units.s, description="Maximum timestep", ptype="nml"),
runtime    =  dict(group_name="RUN",short="runtime", dtype="float64", default=1.e6 | units.s, description="total maximum simulation time", ptype="nml"),
lwarmstart =  dict(group_name="RUN",short="lwarmstart", dtype="bool", default=False, description="flag for cold or warm start", ptype="nml"),
startfile  =  dict(group_name="RUN",short="startfile", dtype="string", default="initdlatestx000y000.001", description="restart file name", ptype="nml"),
trestart   =  dict(group_name="RUN",short="trestart", dtype="float64", default=3600. | units.s, description="restart file write interval", ptype="nml"),
dtav_glob  =  dict(group_name="RUN",short="dtav_glob", dtype="float64", default=60. | units.s, description="sampling interval for statistical routines (multiple of dtmax)", ptype="nml"),
timeav_glob=  dict(group_name="RUN",short="timeav_glob", dtype="float64", default=3600. | units.s, description="write interval for statistical routines (multiple of dtav_glob)", ptype="nml"),
irandom    =  dict(group_name="RUN",short="irandom", dtype="int32", default=0, description="random seed", ptype="nml"),
krand      =  dict(group_name="RUN",short="krand", dtype="int32", default=96, description="max level for random (0<krand<kmax)", ptype="nml"),
randqt     =  dict(group_name="RUN",short="randqt", dtype="float64", default=1.e-5 | units.shu, description="amplitude for random qt", ptype="nml"),
randthl    =  dict(group_name="RUN",short="randthl", dtype="float64", default=0.1 | units.K, description="amplitude for random thl", ptype="nml"),
nsv        =  dict(group_name="RUN",short="nsv", dtype="int32", default=0, description="number of additional passive scalars", ptype="nml"),
ladaptive  =  dict(group_name="RUN",short="ladaptive", dtype="bool", default=False, description="flag for adaptive timestep dep. on stab criteria", ptype="nml"),
courant    =  dict(group_name="RUN",short="courant", dtype="float64", default=0.7, description="Courant number", ptype="nml"),
peclet     =  dict(group_name="RUN",short="peclet", dtype="float64", default=0.15, description="Peclet number", ptype="nml"),
author     =  dict(group_name="RUN",short="author", dtype="string", default="OMUSE", description="Name of the author", ptype="nml"),
krandumin  =  dict(group_name="RUN",short="krandumin", dtype="int32", default=1, description="Bottom vertical full level of wind randomization", ptype="nml"),
krandumax  =  dict(group_name="RUN",short="krandumax", dtype="int32", default=0, description="Top vertical full level of wind randomization", ptype="nml"),
randu      =  dict(group_name="RUN",short="randu", dtype="float64", default=0.5 | units.m/units.s, description="amplitude of wind randomization", ptype="nml"),
nprocx     =  dict(group_name="RUN",short="nprocx", dtype="int32", default=0 , description="Number of processes in the x direction (0=means let MPI decide)", ptype="nml"),
nprocy     =  dict(group_name="RUN",short="nprocy", dtype="int32", default=0 , description="Number of processes in the y direction (0=means let MPI decide)", ptype="nml"),
#
ldelta     =  dict(group_name="NAMSUBGRID", short="ldelta", dtype="bool", default=False, description="Switch for diminished sfs in stable flow", ptype="nml"),
lmason     =  dict(group_name="NAMSUBGRID", short="lmason", dtype="bool", default=False, description="Switch for decreased length scale near the surface", ptype="nml"),
cf         =  dict(group_name="NAMSUBGRID", short="cf", dtype="float64", default=2.5, description="Filter constant", ptype="nml"),
cn         =  dict(group_name="NAMSUBGRID", short="cn", dtype="float64", default=0.76, description="Subfilter scale parameter", ptype="nml"),
Rigc       =  dict(group_name="NAMSUBGRID", short="Rigc", dtype="float64", default=0.25, description="Critical Richardson number", ptype="nml"),
Prandtl    =  dict(group_name="NAMSUBGRID", short="Prandtl", dtype="float64", default=1./3., description="Prandtl number", ptype="nml"),
lsmagorinsky=  dict(group_name="NAMSUBGRID", short="lsmagorinsky", dtype="bool", default=False, description="Switch for smagorinsky subgrid scheme ", ptype="nml"),
cs         =  dict(group_name="NAMSUBGRID", short="cs", dtype="float64", default=-1., description="Smagorinsky constant", ptype="nml"),
nmason     =  dict(group_name="NAMSUBGRID", short="nmason", dtype="float64", default=2., description="Exponent in Mason correction function", ptype="nml"),
sgs_surface_fix=  dict(group_name="NAMSUBGRID", short="sgs_surface_fix", dtype="bool", default=True, description="Switch to apply a fix to the coupling of SFS TKE to the surface (experimental)", ptype="nml"),
#
# skipped NAMAGScross
# skipped NAMBUDGET
# skipped NAMBULKMICROSTAT
#
lcanopy    =  dict(group_name="NAMCANOPY", short="lcanopy", dtype="bool", default=False, description="Switch to represent canopy drag", ptype="nml"),
ncanopy    =  dict(group_name="NAMCANOPY", short="ncanopy", dtype="int32", default=10, description="Amount of layers that contain canopy", ptype="nml"),
cd         =  dict(group_name="NAMCANOPY", short="cd", dtype="float64", default=0.15, description="canopy drag coefficient", ptype="nml"),
lai        =  dict(group_name="NAMCANOPY", short="lai", dtype="float64", default=2., description="One-sided plant area index of the canopy", ptype="nml"),
lpaddistr  =  dict(group_name="NAMCANOPY", short="lpaddistr", dtype="bool", default=False, description="Switch to make use of customized plant area density (prescribed at half levels in paddistr.inp)", ptype="nml"),
npaddistr  =  dict(group_name="NAMCANOPY", short="npaddistr", dtype="int32", default=11, description="number of half-levels in paddistr.inp", ptype="nml"),
wth_can    =  dict(group_name="NAMCANOPY", short="wth_can", dtype="float64", default=0. | units.K*units.m/units.s, description="prescribed SH canopy flux(at top)", ptype="nml"),
wqt_can    =  dict(group_name="NAMCANOPY", short="wqt_can", dtype="float64", default=0. | units.shu*units.m/units.s, description="Prescribed LE canopy flux (at top)", ptype="nml"),
wsv_can    =  dict(group_name="NAMCANOPY", short="wsv_can", dtype="float64", default=0. | units.ppb*units.m/units.s, description="Prescribed scalar flux (at top)", ptype="nml"),
wth_total  =  dict(group_name="NAMCANOPY", short="wth_total", dtype="bool", default=False, description="If true, wth_can includes the surface flux", ptype="nml"),
wqt_total  =  dict(group_name="NAMCANOPY", short="wqt_total", dtype="bool", default=False, description="If true, wqt_can includes the surface flux", ptype="nml"),
wsv_total  =  dict(group_name="NAMCANOPY", short="wsv_total", dtype="bool", default=False, description="If true, wsv_can includes the surface flux", ptype="nml"),
wth_alph   =  dict(group_name="NAMCANOPY", short="wth_alph", dtype="float64", default=0., description="Decay constant for SH with integrated PAD", ptype="nml"),
wqt_alph   =  dict(group_name="NAMCANOPY", short="wqt_alph", dtype="float64", default=0., description="Decay constant for LE with integrated PAD", ptype="nml"),
wsv_alph   =  dict(group_name="NAMCANOPY", short="wsv_alph", dtype="float64", default=0., description="Decay constant for scalar fluxes with integrated PAD", ptype="nml"),
#
# skipped NAMCAPE
# skipped NAMCHECKSIM
# skipped NAMCHEM
# skipped NAMCLOUDFIELD
#
lcross     =  dict(group_name="NAMCROSSSECTION", short="lcross", dtype="bool", default=False, description="Switch for dumping of crossections", ptype="nml"),
lbinary    =  dict(group_name="NAMCROSSSECTION", short="lbinary", dtype="bool", default=False, description="Switch for dumping of crossections in binary files", ptype="nml"),
dtav       =  dict(group_name="NAMCROSSSECTION", short="dtav", dtype="int32", default=60 | units.s, description="Time interval for dumping crosssections", ptype="nml"),
crossheight=  dict(group_name="NAMCROSSSECTION", short="crossheight", dtype="int32", default=2, description="indices of layers of the horizontal crossections (max 100)", ptype="nml"),
crossplane =  dict(group_name="NAMCROSSSECTION", short="crossplane", dtype="int32", default=2, description="location of the  certical xz plane", ptype="nml"),
crossortho =  dict(group_name="NAMCROSSSECTION", short="crossortho", dtype="int32", default=2, description="location of the vertical yz plane", ptype="nml"),
#
# skipped NAMDE
# skipped NAMFIELDDUMP
# skipped NAMGENSTAT
# skipped NAMHETEROSTATS
# skipped NAMLSMCROSSSECTION
# skipped NAMLSMSTAT
#
imicro    =  dict(group_name="NAMMICROPHYSICS", short="imicro", dtype="int32", default=0, description="Flag for the microphysics scheme: (0 [no microphysics], 1 [drizzle], 2[bulk], 3[bin (inactive)], 5[simple ice] or  10[user])", ptype="nml"),
l_sb      =  dict(group_name="NAMMICROPHYSICS", short="l_sb", dtype="bool", default=False, description="toggle for (False) Khairoutdinov and Kogan, 2000 or (True) Seifert & Beheng 2001,2006 scheme",ptype="nml"),
l_sedc    =  dict(group_name="NAMMICROPHYSICS", short="l_sedc", dtype="bool", default=True, description="switch for cloud droplet sedimentation", ptype="nml"),
l_rain    =  dict(group_name="NAMMICROPHYSICS", short="l_rain", dtype="bool", default=True, description="switch for rain formation and evolution", ptype="nml"),
l_mur_cst =  dict(group_name="NAMMICROPHYSICS", short="l_mur_cst", dtype="bool", default=False, description="switch to use constant mu_r (in raindrop gamma dist)", ptype="nml"),
l_berry   =  dict(group_name="NAMMICROPHYSICS", short="l_berry", dtype="bool", default=True, description="Berry-Hsie autoconversion instead of Kessler-Lin", ptype="nml"),
l_graupel =  dict(group_name="NAMMICROPHYSICS", short="l_graupel", dtype="bool", default=True, description="switch for graupel", ptype="nml"),
l_warm    =  dict(group_name="NAMMICROPHYSICS", short="l_warm", dtype="bool", default=False, description="Check: run ice micro in warm mode", ptype="nml"),
mur_cst   =  dict(group_name="NAMMICROPHYSICS", short="mur_cst", dtype="float64", default=5, description="Value for mu_r, a shape parameter for the rain drop dens distr.", ptype="nml"),
Nc_0      =  dict(group_name="NAMMICROPHYSICS", short="Nc_0", dtype="float64", default=7.e7, description="Initial number of cloud droplets", ptype="nml"),
sig_g     =  dict(group_name="NAMMICROPHYSICS", short="sig_g", dtype="float64", default=1.34, description="Geometric standard deviation of the cloud droplet drop size distribution", ptype="nml"),
sig_gr    =  dict(group_name="NAMMICROPHYSICS", short="sig_gr", dtype="float64", default=1.5, description="Geometric standard deviation of the rain droplet drop size distribution", ptype="nml"),
#
# skipped NAMNETCDFSTATS
# skipped NAMNUDGE
# skipped NAMPARTICLES
# skipped NAMprojection (obsolete)
# skipped NAMquadrant
# skipped NAMRADIATION
# skipped NAMRADSTAT
# skipped NAMSAMPLING
# skipped NAMSIMPLEICESTAT
# skipped NAMSTATTEND
# skipped NAMSTRESS
#
# NAMSURFACE
#
# skipped NAMTILT
# skipped NAMTIMESTAT
# 
)

import f90nml
from collections import defaultdict
from omuse.units.quantities import new_quantity, to_quantity, is_quantity


# CodeWithNamelistParameters
# 
# namelist_parameters=dict(
#   parametername  =  dict(group_name="name", short="codename", dtype="int32", default=64, description="description", ptype="nml" [, set_name="name"]),
# )

class CodeWithNamelistParameters(object):
    def __init__(self, namelist_parameters):
        self._namelist_parameters=namelist_parameters
    
    def define_parameters(self,object):
        for name,p in self._namelist_parameters.iteritems():
            if p["ptype"] in ["nml", "nml+normal"]:
                parameter_set_name=p.get("set_name", "parameters_"+p["group_name"])
                object.add_interface_parameter( name, p["description"], p["default"], "before_set_interface_parameter", parameter_set=parameter_set_name)

    def read_namelist_parameters(self, inputfile):

        self._nml_file=inputfile
        self._nml_params = f90nml.read(inputfile)

        for group, d in self._nml_params.iteritems():
            for name, val in d.iteritems():
                if name in self._namelist_parameters:
                    group_name=self._namelist_parameters[name]["group_name"]
                    parameter_set_name=self._namelist_parameters[name].get("set_name", "parameters_"+group_name)
                    parameter_set=getattr(self, parameter_set_name)
                    if is_quantity(self._namelist_parameters[name]["default"]):
                        setattr(parameter_set, name, new_quantity(val, to_quantity(self._namelist_parameters[name]["default"]).unit) )
                    else:
                        setattr(parameter_set, name, val )
                else:
                    print "'%s' of group '%s' not in the namelist_parameters"%(name, group)

    def write_namelist_parameters(self, outputfile, do_patch=False, nml_file=None):
        patch=defaultdict( dict )
        for name, v in self._namelist_parameters.iteritems():
            group_name=v["group_name"]
            group=patch[group_name]
            short=v["short"]
            parameter_set_name=v.get("set_name", "parameters_"+group_name)
            parameter_set=getattr(self, parameter_set_name)
            if getattr(parameter_set, name) is None:  # omit if value is None
                continue
            if is_quantity(self._namelist_parameters[name]["default"]):
                group[short]=to_quantity(getattr(parameter_set, name)).value_in(self._namelist_parameters[name]["default"].unit)
            else:
                group[short]=getattr(parameter_set, name)
        
        if do_patch:
            f90nml.patch(nml_file or self._nml_file,patch,outputfile)
        else:
            f90nml.write(patch, outputfile, force=True)      
if __name__=="__main__":
    pass
    
    
  
  
