from omuse.units import units

import f90nml

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
xtime =  dict(group_name="DOMAIN", short="xday", dtype="float64", default=1 | units.hour, description="start time, UTC time of the day", ptype="nml"),
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
wsvsurf    =  dict(group_name="PHYSICS", short="wqsurf", dtype="float64", default=0. | units.ppb*units.m/units.s, description ="Flux of scalars at the surface", ptype="nml"),
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
timerad   =   dict(group_name="PHYSICS",short="time_rad", dtype="float64", default=0. | units.s, description="Sampling interval of radiation scheme", ptype="nml"),
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
runtime    =  dict(group_name="RUN",short="runtime", dtype="float64", default=300. | units.s, description="total maximum simulation time", ptype="nml"),
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


)


#~ class namelist_parameters(object):
    #~ def __init__(self, parameters):
        #~ self.parameter_dict=dict()
        #~ self.namelist_dict=dict()

    #~ def read_namelist_file(self, filename):
        #~ self.namelist_dict=f90nml.read(filename)
        #~ for group in self.namelist_dict:
            #~ for key in self.namelist_dict:
                #~ self.parameter_dict[key][group_name]=group
                #~ self.parameter_dict[key][short]=key
                #~ self.parameter_dict[key][value]=
      
    #~ def write_namelist_file(self, filename):
      
      
      
if __name__=="__main__":
    pass
    #~ code=Dales(input_file="namoptions.001")
    
    #~ print code.parameters
    
    
  
  
