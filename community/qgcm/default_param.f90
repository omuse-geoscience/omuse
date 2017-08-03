!!
!! Input parameters for Southern Ocean configuration
!! -------------------------------------------------
!! N.B. At most the first 80 columns of this file are significant
!! All lines with a "!" in column 1 are ignored; there must be
!! no blank lines until all the parameter values have been read
!!
!! Timing parameters
!! -----------------
trun=30./365.          !! timestep in "years" 
dta= 180.0D0          !! dta     = Atmos. timestep in seconds
nstr= 3               !! nstr    = Timestep ratio dto/dta
!!
!! Physical parameters
!! -------------------
dxo= 5.0D3           !! dxo     = Ocean grid spacing (m)
delek= 2.0d0           !! delek   = Ocean bottom Ekman layer thickness (m)
cdat= 1.3d-3          !! cdat    = Air-Sea momentum drag coefficient (quad)
rhoat= 1.0d0           !! rhoat   = Atmospheric density (kg m^-3)
rhooc= 1.0d3           !! rhooc   = Oceanic density     (kg m^-3)
cpat= 1.0d3           !! cpat    = Atmos. specific heat capacity (J kg^-1 K^-1)
cpoc= 4.0d3           !! cpoc    = Ocean  specific heat capacity (J kg^-1 K^-1)
bccoat= 1.0d0           !! bccoat  = Mixed BC coefficient for atmos. (nondim.)
bccooc= 0.20d0          !! bccooc  = Mixed BC coefficient for ocean  (nondim.)
xcexp= 1.0d0           !! xcexp   = coupling coefficient x
ycexp= 1.0d0           !! ycexp   = coupling coefficient y
!!
!! Data dumping/averaging parameters
!! ---------------------------------
valday= 0.250d0         !! valday  = Validity checking interval (days)
odiday= 5.00d0          !! odiday  = Ocean  data dump interval  (days)
adiday= 5.00d0          !! adiday  = Atmos. data dump interval  (days)
dgnday= 1.00d0          !! dgnday  = Diagnostics data dump interval (days)
prtday= 5.00d0          !! prtday  = Print to std output interval (days)
resday= 5.0d0           !! resday  = Restart dump interval  (days) (zero => off)
nsko= 2               !! nsko    = Subsampling interval for ocean output
nska= 1               !! nska    = Subsampling interval for atmos. output
dtavat= 0.25d0          !! dtavat  = Atmos. averaging interval (days) (zero => off)
dtavoc= 1.0d0           !! dtavoc  = Ocean  averaging interval (days) (zero => off)
dtcovat= 1.0d0           !! dtcovat = Atmos. covar interval (days) (zero => off)
dtcovoc= 5.0d0           !! dtcovoc = Ocean  covar interval (days) (zero => off)
!!
!! Mixed layer parameters
!! ----------------------
xlamda= 35.0d0          !! xlamda  = Sensible and latent transfer
hmoc= 100.0d0         !! hmoc    = Fixed ocean  mixed layer depth  (m)
st2d= 100.0d0         !! st2d    = sst grad-sqd diffusivity (m^2 s^-1)
st4d= 2.0d+9          !! st4d    = sst grad-4th diffusivity (m^4 s^-1)
hmat= 1000.0d0        !! hmat    = Fixed atmos. mixed layer depth  (m)
hmamin= 100.0d0         !! hmamin  = Min. allowed atmos. m. l. depth (m)
ahmd= 2.0d5           !! ahmd    = atmos. hmix diffusivity  (m^2 s^-1)
at2d= 2.5d4           !! at2d    = ast grad-sqd diffusivity (m^2 s^-1)
at4d= 2.0d14          !! at4d    = ast grad-4th diffusivity (m^4 s^-1)
hmadmp= 0.15d0          !! hmadmp  = Atmos. m.l. damping constant
!!
!! Radiation parameters
!! --------------------
fsbar= -210.0d0        !! fsbar   = Mean radiative forcing  (W m^-2)
fspamp= 80.0d0          !! fspamp  = Radiation perturbation magnitude (W m^-2) (i.e. peak-trough variation, always +ve)
zm= 2.0d2           !! zm      = Optical depth in a.m.l.  (m)
zopt=(/ 2.0d4,  2.0d4,  3.0d4/)        !! zopt(k) = Optical depth in layer k (m)
gamma= 1.0d-2          !! gamma   = Adiabatic lapse rate (K m^-1)
!!
!! Oceanic QG layer parameters
!! ---------------------------
ah2oc=(/  0.0d0 ,   0.0d0 ,   0.0d0 /)  !! ah2oc(k)  = Del-sqd coefft for ocean layer k (m^2 s^-1)
ah4oc=(/ 2.0d+9 ,  2.0d+9 ,  2.0d+9 /)  !! ah4oc(k)  = Del-4th coefft for ocean layer k (m^4 s^-1)
tabsoc=(/ 287.0d0 , 282.0d0 , 276.0d0 /) !! tabsoc(k) = Potential temp. of ocean layer k (K)
hoc=(/ 350.0d0 , 750.0d0, 2900.0d0 /) !! hoc(k)    = Thickness of ocean layer k (m)
gpoc=(/ 0.0150d0 , 0.0075d0 /)        !! gpoc(i)   = Reduced gravity for ocean interface i (m s^-2)
!!
!! Atmospheric QG layer parameters
!! -------------------------------
ah4at= (/1.5d+14 , 1.5d+14 , 1.5d+14 /)  !! ah4at(k)  = Del-4th coefft for atmos. layer k (m^4 s^-1)
tabsat= (/ 330.0d0 , 340.0d0 , 350.0d0 /)  !! tabsat(k) = Temperature for atmos. layer k (K)
hat=(/ 2000.d0 , 3000.d0 , 4000.d0 /) !! hat(k)    = Thickness of atmos. layer k (m)
gpat=(/ 1.2d0 ,    0.4d0 /)           !! gpat(i)   = Reduced gravity for atmos. interface i (m s^-2)
!!
!! Configuration/control parameters
!! --------------------------------
!! Initial state. Options are: zero, rbal or the name of a restart file
name="zero"
outdir="./outdata"
!!./lastday.nc
!! Ocean  topography. Options are: flat, define or the name of a file
topocname="flat"
!! Atmos. topography. Options are: flat, define or the name of a file
topatname="flat"
outfloc=(/ 1, 1, 1, 1, 1 ,1 ,0 /)              !! outfloc = output flags for ocean
outflat=(/ 1, 1, 1, 1, 1, 1, 1 /)              !! outflat = output flags for atmos.
                            !! outfloc/at(i): An integer vector containing a set
                            !! of flags specifying which output fields are required.
                            !! 1 implies variable is output to netcdf
                            !! file, 0 means not. Variables (in order) are:
                            !! ml temp, p, q, Ekman vel. at T pts,
                            !! interface height, windstress.
                            !! mixed layer thickness
