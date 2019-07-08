from omuse.units import units

# dict of parameters: keys become the OMUSE parameter names, short is the name in the namelist file
namelist_parameters = (
#
    dict(name="itot", group_name="DOMAIN", short="itot", dtype="int32", default=64, description="Number of horizontal grid points in x-direction", ptype="nml"),
    dict(name="jtot", group_name="DOMAIN", short="jtot", dtype="int32", default=64, description="Number of horizontal grid points in y-direction", ptype="nml"),
    dict(name="kmax", group_name="DOMAIN", short="kmax", dtype="int32", default=96, description="Number of vertical grid points", ptype="nml"),
    dict(name="xsize", group_name="DOMAIN", short="xsize", dtype="float64", default=6400. | units.m, description="Horizontal size of the simulated domain in x", ptype="nml"),
    dict(name="ysize", group_name="DOMAIN", short="ysize", dtype="float64", default=6400. | units.m, description="Horizontal size of the simulated domain in y", ptype="nml"),
    dict(name="xlat", group_name="DOMAIN", short="xlat", dtype="float64", default=52 | units.deg, description="latitude of domain", ptype="nml"),
    dict(name="xlon", group_name="DOMAIN", short="xlon", dtype="float64", default=0 | units.deg, description="longitude of domain", ptype="nml"),
#    dict(name="xyear", group_name="DOMAIN", short="xyear", dtype="float64", default=0 | units.yr, description="year, only for time units in netcdf", ptype="nml"),
    dict(name="xday", group_name="DOMAIN", short="xday", dtype="float64", default=1 | units.day, description="start day, number of the day", ptype="nml"),
    dict(name="xtime", group_name="DOMAIN", short="xtime", dtype="float64", default=1 | units.hour, description="start time, UTC time of the day", ptype="nml"),
    dict(name="ksp", group_name="DOMAIN", short="ksp", dtype="float64", default=None, description="lower height of sponge layer in grid cells", ptype="nml"),
#
    dict(name="iexpnr", group_name="RUN", short="iexpnr", dtype="int32", default=1, description="Experiment number", ptype="nml"),
    dict(name="dtmax", group_name="RUN", short="dtmax", dtype="float64", default=20. | units.s, description="Maximum timestep", ptype="nml"),
    dict(name="runtime", group_name="RUN", short="runtime", dtype="float64", default=1.e6 | units.s, description="total maximum simulation time", ptype="nml"),
    dict(name="lwarmstart", group_name="RUN", short="lwarmstart", dtype="bool", default=False, description="flag for cold or warm start", ptype="nml"),
    dict(name="startfile", group_name="RUN", short="startfile", dtype="string", default="initdlatestx000y000.001", description="restart file name", ptype="nml"),
    dict(name="trestart", group_name="RUN", short="trestart", dtype="float64", default=3600. | units.s, description="restart file write interval", ptype="nml"),
    dict(name="dtav_glob", group_name="RUN", short="dtav_glob", dtype="float64", default=60. | units.s, description="sampling interval for statistical routines (multiple of dtmax)", ptype="nml"),
    dict(name="timeav_glob", group_name="RUN", short="timeav_glob", dtype="float64", default=3600. | units.s, description="write interval for statistical routines (multiple of dtav_glob)", ptype="nml"),
    dict(name="irandom", group_name="RUN", short="irandom", dtype="int32", default=43, description="random seed", ptype="nml"),
    dict(name="krand", group_name="RUN", short="krand", dtype="int32", default=96, description="max level for random (0<krand<kmax)", ptype="nml"),
    dict(name="randqt", group_name="RUN", short="randqt", dtype="float64", default=1.e-5 | units.shu, description="amplitude for random qt", ptype="nml"),
    dict(name="randthl", group_name="RUN", short="randthl", dtype="float64", default=0.1 | units.K, description="amplitude for random thl", ptype="nml"),
    dict(name="nsv", group_name="RUN", short="nsv", dtype="int32", default=0, description="number of additional passive scalars", ptype="nml"),
    dict(name="ladaptive", group_name="RUN", short="ladaptive", dtype="bool", default=False, description="flag for adaptive timestep dep. on stab criteria", ptype="nml"),
    dict(name="courant", group_name="RUN", short="courant", dtype="float64", default=0.7, description="Courant number", ptype="nml"),
    dict(name="peclet", group_name="RUN", short="peclet", dtype="float64", default=0.15, description="Peclet number", ptype="nml"),
    dict(name="author", group_name="RUN", short="author", dtype="string", default="OMUSE", description="Name of the author", ptype="nml"),
    dict(name="krandumin", group_name="RUN", short="krandumin", dtype="int32", default=1, description="Bottom vertical full level of wind randomization", ptype="nml"),
    dict(name="krandumax", group_name="RUN", short="krandumax", dtype="int32", default=0, description="Top vertical full level of wind randomization", ptype="nml"),
    dict(name="randu", group_name="RUN", short="randu", dtype="float64", default=0.5 | units.m/units.s, description="amplitude of wind randomization", ptype="nml"),
    dict(name="nprocx", group_name="RUN", short="nprocx", dtype="int32", default=0, description="Number of processes in the x direction (0=means let MPI decide)", ptype="nml"),
    dict(name="nprocy", group_name="RUN", short="nprocy", dtype="int32", default=0, description="Number of processes in the y direction (0=means let MPI decide)", ptype="nml"),
#
    dict(name="llsadv", group_name="DYNAMICS", short="llsadv", dtype="bool", default=False, description="Switch for large scale forcings", ptype="nml"),
    dict(name="lqlnr", group_name="DYNAMICS", short="lqlnr", dtype="bool", default=True, description="Switch for Newton-Raphson approximation of the liquid water content", ptype="nml"),
    dict(name="lnoclouds", group_name="DYNAMICS", short="lnoclouds", dtype="bool", default=False, description="Switch to disable q_l calculations", ptype="nml"),
    dict(name="cu", group_name="DYNAMICS", short="cu", dtype="float64", default=0. | units.m/units.s, description="Transformation velocity of the Galilei transformation in x-direction",ptype="nml"),
    dict(name="cv", group_name="DYNAMICS", short="cv", dtype="float64", default=0. | units.m/units.s, description="Transformation velocity of the Galilei transformation in y-direction",ptype="nml"),
    dict(name="iadv_mom", group_name="DYNAMICS", short="iadv_mom", dtype="int32", default=5, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
    dict(name="iadv_tke", group_name="DYNAMICS", short="iadv_tke", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
    dict(name="iadv_thl", group_name="DYNAMICS", short="iadv_thl", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
    dict(name="iadv_qt", group_name="DYNAMICS", short="iadv_qt", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
    dict(name="iadv_sv", group_name="DYNAMICS", short="iadv_sv", dtype="int32", default=-1, description="Advection scheme for momentum (1 [1st order upwind],2 [2nd order central diff],5 [5th order upwind],52 [horiz 5th, vert 2nd], 55 [hybrid],6 [6th order cebntral diff],62 [6th horiz. , 2nd vert.] or 7 [Kappa])",ptype="nml"),
    dict(name="ibas_prf", group_name="DYNAMICS", short="ibas_prf", dtype="int32", default=3, description="Flag for density calculations, -1 = copy iadv_mom, 1= based constant theta_0, 2 = Boussinesq-like, 3=Standard lapse rate, based on surface temp, 4= standard lr, tsurf=15C, 5=user defined",ptype="nml"),
    dict(name="lambda_crit", group_name="DYNAMICS", short="lambda_crit", dtype="float64", default=100, description="Maximum value for the smoothness.",ptype="nml"),
#
    dict(name="z0", group_name="PHYSICS", short="z0", dtype="float64", default=1.6e-4 | units.m, description ="surface roughness (negative = ", ptype="nml"),
    dict(name="ustin", group_name="PHYSICS", short="ustin", dtype="float64", default=0.32 | units.m/units.s, description ="Prescribed friction velocity (negative = ", ptype="nml"),
    dict(name="wtsurf", group_name="PHYSICS", short="wtsurf", dtype="float64", default=0. | units.K*units.m/units.s, description ="Flux of liq. water pot. temp. at the surface (negative = ", ptype="nml"),
    dict(name="wqsurf", group_name="PHYSICS", short="wqsurf", dtype="float64", default=0. | units.K*units.m/units.s, description ="Flux of total water content at the surface (negative = ", ptype="nml"),
    dict(name="wsvsurf", group_name="PHYSICS", short="wsvsurf", dtype="float64", default=0. | units.ppb*units.m/units.s, description ="Flux of scalars at the surface", ptype="nml"),
    dict(name="thls", group_name="PHYSICS", short="thls",dtype="float64",default=298.5 | units.K, description="Liquid water potential temperature at the surface", ptype="nml"),
    dict(name="ps", group_name="PHYSICS", short="ps",dtype="float64",default=101540. | units.Pa, description="surface pressure", ptype="nml"),
    dict(name="isurf", group_name="PHYSICS", short="isurf",dtype="int32",default=4, description="surface parametrization (1 [interactive, use rad.],2 [forced surf T, flux calculated],3 [forced mom, moisture and heat flux],4 [forced moisture and heat flux] or 10 [user def]", ptype="nml"),
    dict(name="ltimedep", group_name="PHYSICS", short="ltimedep", dtype="bool", default=False, description="switch for timedependent fluxes", ptype="nml"),
    dict(name="lcoriol", group_name="PHYSICS", short="lcoriol", dtype="bool", default=True, description="switch for coriolis force", ptype="nml"),
    dict(name="igrw_damp", group_name="PHYSICS", short="igrw_damp", dtype="int32", default=2, description="flag for grav wave damping ( -1 [nudge to lscale.inp], 0 [no damping], 1 [fast damping to av. and slow damping of av to geowind], 2[fast damping to geowind] or 3[fast damping to av wind]", ptype="nml"),
    dict(name="geodamptime", group_name="PHYSICS", short="geodamptime", dtype="float64", default=7200. | units.s, description="timescale for nudging to geowind", ptype="nml"),
    dict(name="lmomsubs", group_name="PHYSICS", short="lmomsubs", dtype="bool", default=False, description="switch to apply subsidence on momentum", ptype="nml"),
    dict(name="lmoist", group_name="PHYSICS", short="lmoist", dtype="bool", default=True, description="switch for calculation of moisture field", ptype="nml"),
    dict(name="chi_half", group_name="PHYSICS", short="chi_half", dtype="float64", default=0.5, description="Wet, dry or intermediate (default) mixing over the cloud edge", ptype="nml"),
    dict(name="timerad", group_name="PHYSICS", short="timerad", dtype="float64", default=0. | units.s, description="Sampling interval of radiation scheme", ptype="nml"),
    dict(name="iradiation", group_name="PHYSICS", short="iradiation", dtype="int32", default=0, description="flag for rad. calc. (0 [no rad], 1 [full], 2 [parametrized], 3 [simple for land surface], 4 [rapid RTM] or 10 [user]", ptype="nml"),
    dict(name="useMcICA", group_name="PHYSICS", short="useMcICA", dtype="bool", default=True, description="switch for the Monte Carlo Indep. Column approach (only for iradiation=1)", ptype="nml"),
    dict(name="rad_ls", group_name="PHYSICS", short="rad_ls", dtype="bool", default=True, description="switch for prescribed radiation forcing", ptype="nml"),
    dict(name="rad_longw", group_name="PHYSICS", short="rad_longw", dtype="bool", default=True, description="switch for parameterized long wave radiation forcing", ptype="nml"),
    dict(name="rad_shortw", group_name="PHYSICS", short="rad_shortw", dtype="bool", default=True, description="switch for parameterized short wave radiation forcing", ptype="nml"),
    dict(name="rad_smoke", group_name="PHYSICS", short="rad_smoke", dtype="bool", default=False, description="switch for longwave divergence for smoke cloud", ptype="nml"),
    dict(name="rka", group_name="PHYSICS", short="rka",dtype="float64", default=130. | units.m**2/units.kg, description="Extinction coefficient (if iradiation=2)", ptype="nml"),
    dict(name="dlwbot", group_name="PHYSICS", short="dlwbot",dtype="float64", default=0. | units.W/units.m**2, description="long wave rad flux jump at cloud bottom", ptype="nml"),
    dict(name="dlwtop", group_name="PHYSICS", short="dlwtop",dtype="float64", default=74. | units.W/units.m**2, description="long wave rad flux jump at cloud top", ptype="nml"),
    dict(name="sw0", group_name="PHYSICS", short="sw0",dtype="float64", default=1100. | units.W/units.m**2, description="Direct solar radiative component at cloud top (assumes zero diffusive contr.)", ptype="nml"),
    dict(name="gc", group_name="PHYSICS", short="gc",dtype="float64", default=0.85, description="Asymmetry factor for droplet scattering", ptype="nml"),
    dict(name="reff", group_name="PHYSICS", short="reff",dtype="float64", default=1.e-5 | units.m, description="cloud drop effective radius", ptype="nml"),
    dict(name="isvsmoke", group_name="PHYSICS", short="isvsmoke",dtype="int32", default=1, description="scalar field to be used for optical depth calc. (when rad_smoke=True)", ptype="nml"),
# omitted lforce_user
    dict(name="lcloudshading", group_name="PHYSICS", short="lcloudshading", dtype="bool", default=False, description="switch to let clouds shade the surface for rad_lsm", ptype="nml"),
    dict(name="lrigidlid", group_name="PHYSICS", short="lrigidlid", dtype="bool", default=False, description="Switch to enable simulations with a rigid lid", ptype="nml"),
    dict(name="unudge", group_name="PHYSICS", short="unudge", dtype="float64", default=1., description="Nudging factor if igrw_damp is -1", ptype="nml"),
#
    dict(name="ldelta", group_name="NAMSUBGRID", short="ldelta", dtype="bool", default=False, description="Switch for diminished sfs in stable flow", ptype="nml"),
    dict(name="lmason", group_name="NAMSUBGRID", short="lmason", dtype="bool", default=False, description="Switch for decreased length scale near the surface", ptype="nml"),
    dict(name="cf", group_name="NAMSUBGRID", short="cf", dtype="float64", default=2.5, description="Filter constant", ptype="nml"),
    dict(name="cn", group_name="NAMSUBGRID", short="cn", dtype="float64", default=0.76, description="Subfilter scale parameter", ptype="nml"),
    dict(name="Rigc", group_name="NAMSUBGRID", short="Rigc", dtype="float64", default=0.25, description="Critical Richardson number", ptype="nml"),
    dict(name="Prandtl", group_name="NAMSUBGRID", short="Prandtl", dtype="float64", default=1./3., description="Prandtl number", ptype="nml"),
    dict(name="lsmagorinsky", group_name="NAMSUBGRID", short="lsmagorinsky", dtype="bool", default=False, description="Switch for smagorinsky subgrid scheme ", ptype="nml"),
    dict(name="cs", group_name="NAMSUBGRID", short="cs", dtype="float64", default=-1., description="Smagorinsky constant", ptype="nml"),
    dict(name="nmason", group_name="NAMSUBGRID", short="nmason", dtype="float64", default=2., description="Exponent in Mason correction function", ptype="nml"),
#    dict(name="sgs_surface_fix", group_name="NAMSUBGRID", short="sgs_surface_fix", dtype="bool", default=True, description="Switch to apply a fix to the coupling of SFS TKE to the surface (experimental)", ptype="nml"), # removed in dales 4.2

#
# skipped NAMAGScross
# skipped NAMBUDGET
# skipped NAMBULKMICROSTAT
#
    dict(name="lcanopy", group_name="NAMCANOPY", short="lcanopy", dtype="bool", default=False, description="Switch to represent canopy drag", ptype="nml"),
    dict(name="ncanopy", group_name="NAMCANOPY", short="ncanopy", dtype="int32", default=10, description="Amount of layers that contain canopy", ptype="nml"),
    dict(name="cd", group_name="NAMCANOPY", short="cd", dtype="float64", default=0.15, description="canopy drag coefficient", ptype="nml"),
    dict(name="lai", group_name="NAMCANOPY", short="lai", dtype="float64", default=2., description="One-sided plant area index of the canopy", ptype="nml"),
    dict(name="lpaddistr", group_name="NAMCANOPY", short="lpaddistr", dtype="bool", default=False, description="Switch to make use of customized plant area density (prescribed at half levels in paddistr.inp)", ptype="nml"),
    dict(name="npaddistr", group_name="NAMCANOPY", short="npaddistr", dtype="int32", default=11, description="number of half-levels in paddistr.inp", ptype="nml"),
    dict(name="wth_can", group_name="NAMCANOPY", short="wth_can", dtype="float64", default=0. | units.K*units.m/units.s, description="prescribed SH canopy flux(at top)", ptype="nml"),
    dict(name="wqt_can", group_name="NAMCANOPY", short="wqt_can", dtype="float64", default=0. | units.shu*units.m/units.s, description="Prescribed LE canopy flux (at top)", ptype="nml"),
    dict(name="wsv_can", group_name="NAMCANOPY", short="wsv_can", dtype="float64", default=0. | units.ppb*units.m/units.s, description="Prescribed scalar flux (at top)", ptype="nml"),
    dict(name="wth_total", group_name="NAMCANOPY", short="wth_total", dtype="bool", default=False, description="If true, wth_can includes the surface flux", ptype="nml"),
    dict(name="wqt_total", group_name="NAMCANOPY", short="wqt_total", dtype="bool", default=False, description="If true, wqt_can includes the surface flux", ptype="nml"),
    dict(name="wsv_total", group_name="NAMCANOPY", short="wsv_total", dtype="bool", default=False, description="If true, wsv_can includes the surface flux", ptype="nml"),
    dict(name="wth_alph", group_name="NAMCANOPY", short="wth_alph", dtype="float64", default=0., description="Decay constant for SH with integrated PAD", ptype="nml"),
    dict(name="wqt_alph", group_name="NAMCANOPY", short="wqt_alph", dtype="float64", default=0., description="Decay constant for LE with integrated PAD", ptype="nml"),
    dict(name="wsv_alph", group_name="NAMCANOPY", short="wsv_alph", dtype="float64", default=0., description="Decay constant for scalar fluxes with integrated PAD", ptype="nml"),
#
# skipped NAMCAPE
# skipped NAMCHECKSIM
# skipped NAMCHEM
# skipped NAMCLOUDFIELD
#
    dict(name="lcross", group_name="NAMCROSSSECTION", short="lcross", dtype="bool", default=False, description="Switch for dumping of crossections", ptype="nml"),
    dict(name="lbinary", group_name="NAMCROSSSECTION", short="lbinary", dtype="bool", default=False, description="Switch for dumping of crossections in binary files", ptype="nml"),
    dict(name="dtav", group_name="NAMCROSSSECTION", short="dtav", dtype="int32", default=60 | units.s, description="Time interval for dumping crosssections", ptype="nml"),
    dict(name="crossheight", group_name="NAMCROSSSECTION", short="crossheight", dtype="int32", default=2, description="indices of layers of the horizontal crossections (max 100)", ptype="nml"),
    dict(name="crossplane", group_name="NAMCROSSSECTION", short="crossplane", dtype="int32", default=2, description="location of the  certical xz plane", ptype="nml"),
    dict(name="crossortho", group_name="NAMCROSSSECTION", short="crossortho", dtype="int32", default=2, description="location of the vertical yz plane", ptype="nml"),
#
# skipped NAMDE
# skipped NAMFIELDDUMP
# skipped NAMGENSTAT
# skipped NAMHETEROSTATS
# skipped NAMLSMCROSSSECTION
# skipped NAMLSMSTAT
#
    dict(name="imicro", group_name="NAMMICROPHYSICS", short="imicro", dtype="int32", default=0, description="Flag for the microphysics scheme: (0 [no microphysics], 1 [drizzle], 2[bulk], 3[bin (inactive)], 5[simple ice] or  10[user])", ptype="nml"),
    dict(name="l_sb", group_name="NAMMICROPHYSICS", short="l_sb", dtype="bool", default=False, description="toggle for (False) Khairoutdinov and Kogan, 2000 or (True) Seifert & Beheng 2001,2006 scheme",ptype="nml"),
    dict(name="l_sedc", group_name="NAMMICROPHYSICS", short="l_sedc", dtype="bool", default=True, description="switch for cloud droplet sedimentation", ptype="nml"),
    dict(name="l_rain", group_name="NAMMICROPHYSICS", short="l_rain", dtype="bool", default=True, description="switch for rain formation and evolution", ptype="nml"),
    dict(name="l_mur_cst", group_name="NAMMICROPHYSICS", short="l_mur_cst", dtype="bool", default=False, description="switch to use constant mu_r (in raindrop gamma dist)", ptype="nml"),
    dict(name="l_berry", group_name="NAMMICROPHYSICS", short="l_berry", dtype="bool", default=True, description="Berry-Hsie autoconversion instead of Kessler-Lin", ptype="nml"),
    dict(name="l_graupel", group_name="NAMMICROPHYSICS", short="l_graupel", dtype="bool", default=True, description="switch for graupel", ptype="nml"),
    dict(name="l_warm", group_name="NAMMICROPHYSICS", short="l_warm", dtype="bool", default=False, description="Check: run ice micro in warm mode", ptype="nml"),
    dict(name="mur_cst", group_name="NAMMICROPHYSICS", short="mur_cst", dtype="float64", default=5, description="Value for mu_r, a shape parameter for the rain drop dens distr.", ptype="nml"),
    dict(name="Nc_0", group_name="NAMMICROPHYSICS", short="Nc_0", dtype="float64", default=7.e7, description="Initial number of cloud droplets", ptype="nml"),
    dict(name="sig_g", group_name="NAMMICROPHYSICS", short="sig_g", dtype="float64", default=1.34, description="Geometric standard deviation of the cloud droplet drop size distribution", ptype="nml"),
    dict(name="sig_gr", group_name="NAMMICROPHYSICS", short="sig_gr", dtype="float64", default=1.5, description="Geometric standard deviation of the rain droplet drop size distribution", ptype="nml"),
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
    dict(name="isurf", group_name="NAMSURFACE", short="isurf", dtype="int32", default=4, description="Overrides isurf flag of PHYSICS if used", ptype="nml"),
    dict(name="lmostlocal", group_name="NAMSURFACE", short="lmostlocal", dtype="bool", default=False, description="use local Obukhov length", ptype="nml"),
    dict(name="lsmoothflux", group_name="NAMSURFACE", short="lsmoothflux", dtype="bool", default=False, description="Switch to create uniform sensible and latent heat flux over domain", ptype="nml"),
    dict(name="lneutral", group_name="NAMSURFACE", short="lneutral", dtype="bool", default=False, description="Switch to disable stability corrections", ptype="nml"),
    dict(name="z0mav", group_name="NAMSURFACE", short="z0mav", dtype="float64", default=-1. | units.m, description="Roughness length of momentum", ptype="nml"),
    dict(name="z0hav", group_name="NAMSURFACE", short="z0hav", dtype="float64", default=-1. | units.m, description="Roughness length of heat", ptype="nml"),
    #~ dict(name="thls", group_name="NAMSURFACE", short="thls", dtype="float64", default=-1. | units.K, description="surface liquid water potential temperature", ptype="nml"),
    #~ dict(name="ps", group_name="NAMSURFACE", short="ps", dtype="float64", default=-1. | units.Pa, description="surface pressure", ptype="nml"),
    #~ dict(name="ustin", group_name="NAMSURFACE", short="ustin", dtype="float64", default=-1. | units.m/units.s, description="prescribed friction vel.", ptype="nml"),
    #~ dict(name="wtsurf", group_name="NAMSURFACE", short="wtsurf", dtype="float64", default=-1. | units.K*units.m/units.s, description="prescribed kinematic temp flux", ptype="nml"),
    #~ dict(name="wqsurf", group_name="NAMSURFACE", short="wqsurf", dtype="float64", default=-1. | units.shu*units.m/units.s, description="prescribed moisture flux", ptype="nml"),
    #~ dict(name="wsvsurf", group_name="NAMSURFACE", short="wsvsurf", dtype="float64", default=-1. | units.ppb*units.m/units.s, description="prescribed surface scalar fluxes", ptype="nml"),
    dict(name="tsoilav", group_name="NAMSURFACE", short="tsoilav", dtype="float64", default=0. | units.K, description="initial soil temperature for 4 layers only for isurf=1", ptype="nml"),
    dict(name="tsoildeepav", group_name="NAMSURFACE", short="tsoildeepav", dtype="float64", default=0. | units.K, description="soil bottom temperature only for isurf=1", ptype="nml"),
    dict(name="phiwav", group_name="NAMSURFACE", short="phiwav", dtype="float64", default=0. | units.vsmc, description="soil moisture for 4 layers for isurf=1 <0.472, prefer <0.323", ptype="nml"),
    dict(name="rootfav", group_name="NAMSURFACE", short="rootfav", dtype="float64", default=0, description="root fraction, should sum to 1.", ptype="nml"),
    dict(name="Cskinav", group_name="NAMSURFACE", short="Cskinav", dtype="float64", default=-1. | units.J/units.K/units.m**2, description="heat capacity skin layer (isurf=1)", ptype="nml"),
    dict(name="lambdaskinav", group_name="NAMSURFACE", short="lambdaskinav", dtype="float64", default=-1. | units.J/units.s/units.K/units.m**2, description="heat conductivity skin layer (isurf=1)", ptype="nml"),
    dict(name="albedoav", group_name="NAMSURFACE", short="albedoav", dtype="float64", default=-1. , description="albedo (isurf=1)", ptype="nml"),
    dict(name="Qnetav", group_name="NAMSURFACE", short="Qnetav", dtype="float64", default=-1. | units.J/units.s/units.m**2, description="net radiation (iradiation ne 1, isurf=1)", ptype="nml"),
    dict(name="cvegav", group_name="NAMSURFACE", short="cvegav", dtype="float64", default=-1. , description="vegetation cover", ptype="nml"),
    dict(name="Wlav", group_name="NAMSURFACE", short="Wlav", dtype="float64", default=-1. | units.m, description="init. water cover on veg.", ptype="nml"),
    dict(name="rsminav", group_name="NAMSURFACE", short="rsminav", dtype="float64", default=-1. | units.s/units.m, description="minimal vegetation resistance (isurf=1)", ptype="nml"),
    dict(name="rssoilminav", group_name="NAMSURFACE", short="rssoilminav", dtype="float64", default=-1. , description="minimum soil evap. resistance (isurf=1)", ptype="nml"),
    dict(name="LAIav", group_name="NAMSURFACE", short="LAIav", dtype="float64", default=-1. | units.m**2/units.m**2, description="leaf area index (isurf=1)", ptype="nml"),
    dict(name="gDav", group_name="NAMSURFACE", short="gDav", dtype="float64", default=-1. , description="correction for evap. of tall veg. (isurf=1)", ptype="nml"),
    dict(name="rsisurf2", group_name="NAMSURFACE", short="rsisurf2", dtype="float64", default=0. | units.s/units.m, description="vegatiation resistance (isurf=2)", ptype="nml"),
    dict(name="lhetero", group_name="NAMSURFACE", short="lhetero", dtype="bool", default=False, description="switch to apply heterogeneous surfaces", ptype="nml"),
    dict(name="xpatches", group_name="NAMSURFACE", short="xpatches", dtype="int32", default=2, description="number of patches in the x-direction", ptype="nml"),
    dict(name="ypatches", group_name="NAMSURFACE", short="ypatches", dtype="int32", default=1, description="number of patches in the y-direction", ptype="nml"),
    dict(name="land_use", group_name="NAMSURFACE", short="land_use", dtype="int32", default=0, description="indicator for land type (1..10) ", ptype="nml"),
    dict(name="loldtable", group_name="NAMSURFACE", short="loldtable", dtype="bool", default=False, description="flag to use old type surface files ", ptype="nml"),
    dict(name="lrsAgs", group_name="NAMSURFACE", short="lrsAgs", dtype="bool", default=False, description="switch to use A-g_s for resistance calc.", ptype="nml"),
    dict(name="lCO2Ags", group_name="NAMSURFACE", short="lCO2Ags", dtype="bool", default=False, description="calc, CO2 fluxes with A-g_s (if lrsags True)", ptype="nml"),
    dict(name="planttype", group_name="NAMSURFACE", short="planttype", dtype="int32", default=3, description="use (C)3 or (C)4 plants for A-g_s (3 or 4)", ptype="nml"),
    dict(name="lrelaxgc", group_name="NAMSURFACE", short="lrelaxgc", dtype="bool", default=False, description="switch to nudge towards calculated conductivity", ptype="nml"),
    dict(name="kgc", group_name="NAMSURFACE", short="kgc", dtype="float64", default=0.00113 | units.s**-1, description="response rate for stomatal conductivity (<dtmax**-1)", ptype="nml"),
    dict(name="lrelaxci", group_name="NAMSURFACE", short="lrelaxci", dtype="bool", default=False, description="switch to nudge towards calc. internal CO2 concentration", ptype="nml"),
    dict(name="kci", group_name="NAMSURFACE", short="kci", dtype="float64", default=0.00113 | units.s**-1, description="response rate for internal CO2 conc. (<dtmax**-1)", ptype="nml"),
    dict(name="phi", group_name="NAMSURFACE", short="phi", dtype="float64", default=0.472, description="volumetric soil porosity", ptype="nml"),
    dict(name="phifc", group_name="NAMSURFACE", short="phifc", dtype="float64", default=0.323, description="volumetric moisture at field capacity", ptype="nml"),
    dict(name="phiwp", group_name="NAMSURFACE", short="phiwp", dtype="float64", default=0.171, description="volumetric moisture at wilting point", ptype="nml"),
    dict(name="R10", group_name="NAMSURFACE", short="R10", dtype="float64", default=0.23 | units.milli(units.g)/ units.m**2/units.s , description="respiration at 10 degrees C", ptype="nml"),

#
# skipped NAMTILT
# skipped NAMTIMESTAT
# 
)

if __name__=="__main__":
    pass
    
    
  
  
