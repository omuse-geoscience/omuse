DALES interface
===============

This is the OMUSE interface to DALES, the Dutch Atmospheric Large-Eddy Simulation.

Example::
  
  from omuse.community.dales.interface import Dales
  from omuse.units import units

  d = Dales(workdir='dales-workdir', channel_type='sockets', number_of_workers=1)
  
  # set parameters
  d.parameters_DYNAMICS.iadv_mom = 6  # 6th order advection for momentum
  d.parameters_DYNAMICS.iadv_thl = 5  # 5th order advection for scalars, less overshoots than 6th order
  d.parameters_DYNAMICS.iadv_qt = 5
  d.parameters_DYNAMICS.iadv_tke = 5

  # access model variables
  d.fields[:, :, :].U = 0 | units.m / units.s

  # evolve the model
  d.evolve_model(120 | units.s, exactEnd=True)
  

OMUSE models use variables with units, see :ref:`units`.

The parameters are doumented in the DALES documentation: options_. Setting
parameters is possible before the model is time stepped or the
model variables have been accessed, but not after.

Variable grids
--------------

various variables of the DALES model can be accessed from the OMUSE interface. The variables
are grouped in grids, according to their dimensionality.
For example, the U-velocity component in the model component in the model can
be accessed as ``d.fields[i,j,k].U``.

The following grids are defined:

================ ===================================================================
Grid             Description
================ ===================================================================
fields           3D grid
profiles         horizontal averages of the 3D grid (vertical profiles)
nudging_profiles vertical profiles of U, V, THL, QT to nudge the model towards.
forcing_profiles external forcings of U, V, THL, QT - vertical profiles
scalars          surface fluxes and roughness lengths, scalars.
surface_fields   2D surface fields, e.g. liquid water path.
================ ===================================================================

The grids contain the following variables:

================ ===================================================================
Grid             Variables
================ ===================================================================
fields           U, V, W, THL, QT, QL, QL_ice, QR, E12, T, pi, rswd, rswdir, rswdif,
\                rswu, rlwd, rlwu, rswdcs, rswucs, rlwdcs, rlwucs
profiles         U, V, W, THL, QT, QL, QL_ice, QR, E12, T, P, rho
nudging_profiles U, V, THL, QT
forcing_profiles U, V, THL, QT
scalars          wt, wq, z0m, z0h, Ps, QR
surface_fields   LWP, RWP, TWP, ustar, z0m, z0h, tskin, qskin, LE, H, obl
================ ===================================================================


Description of the variables:

============== ====== ===========================================================================
Variable       Unit   Description
============== ====== ===========================================================================
U,V,W          m/s    velocity components in x, y, and z directions
THL            K      liquid water potential temperature
QT             kg/kg  specific total humidity (i.e. vapor and cloud condensate but excluding precipitation)
QL             kg/kg  specific cloud condensate
QL_ice         kg/kg  specific cloud condensate in the form of ice
QR             kg/kg  precipitation
E12            m/s    sqrt(subgrid-scale turbulent kinetic energy)
T              K      temperature
P              Pa     pressure
Ps             Pa     surface pressure
rho            kg/m^3 density
pi             Pa     modified air pressure
rswu,rswd      W/m^2  up-/down-welling short-wave radiative flux
rlwu,rlwd      W/m^2  up-/down-welling long-wave radiative flux
r{s,l}w{u,d}cs W/m^2  up/down-welling long/short-wave radiative flux for clear sky
rswdir         W/m^2  downwelling shortwave direct radiative flux
rswdif         W/m^2  downwelling shortwave diffusive radiative flux
wt             K m/s  surface flux of heat
wq             m/s    surface flux of humidity
z0m, z0h       m      surface roughness length
ustar          m/s    friction velocity
LWP, RWP, TWP  kg/m^2 liquid-, rain-, total water path
tskin          K      skin temperature
qskin          kg/kg  skin humidity
H, LE          W/m^2  sensible and latent heat fluxes   
obl            m      Obukhov length
============== ====== ===========================================================================


The forcing and nudging profiles, and the 3D grids for prognostic variables can be read and written.
Diagnosed quantities, horizontal averages, and water paths can only be read.

Note: the indexing brackets can also be placed on the variable name instead of
on the grid name, e.g. ``d.fields.U[:, :, :]`` vs ``d.fields[:, :, :].U``.
Avoid using the first form, since assigning new values to the grid this way
does not work.



Example scripts
---------------

Some examples of using Dales with OMUSE are bundled included in
omuse/src/omuse/community/dales/example/ .

bubble.py
^^^^^^^^^

A Dales experiment with a warm air bubble. Shows model initialization, setting model parameters, setting initial conditions, evolving the model, and retreiving the state.

async.py
^^^^^^^^

Shows how to use OMUSE's asynchronous function calls, to let several model instances run concurrently.
See also :ref:`asynchronous`.


.. autoclass:: omuse.community.dales.interface.Dales
   :members:

   
      
Links and references
--------------------

* The official DALES git repository_
* DALES manual_
* DALES namelist options_
* Dales model description paper: Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications, T. Heus et al, `Geosci. Model Dev., 3, 415-444, 2010 <https://doi.org/10.5194/gmd-3-415-2010>`_

.. _repository: https://github.com/dalesteam/dales
.. _manual: https://github.com/dalesteam/dales/blob/master/utils/doc/input/dales-manual.pdf
.. _options: https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf
