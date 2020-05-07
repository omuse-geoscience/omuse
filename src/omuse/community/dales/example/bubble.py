"""
A warm air bubble experiment using the DALES Python interface

This script performs the following steps:
* create a Dales object
* set parameters
* create an initial state with a warm air bubble on a constant liquid-water-potential temperature (thl) background
* perform time evolution
* plot a time series of the model using matplotlib.

"""

import numpy
import matplotlib.pyplot as plt
from omuse.community.dales.interface import Dales
from omuse.units import units

# create Dales object
d=Dales(workdir='daleswork', channel_type='sockets', number_of_workers=1)
# add redirection='none' to see the model log messages

# Set parameters: domain size and resolution, advection scheme
d.parameters_DOMAIN.itot = 32  # number of grid cells in x
d.parameters_DOMAIN.jtot = 32  # number of grid cells in y
d.parameters_DOMAIN.xsize = 6400 | units.m
d.parameters_DOMAIN.ysize = 6400 | units.m
d.parameters_DYNAMICS.iadv_mom = 6 # 6th order advection for momentum
d.parameters_DYNAMICS.iadv_thl = 5 # 5th order advection for temperature
d.parameters_RUN.krand = 0  # initial state randomization off

d.parameters_RUN.ladaptive = True
d.parameters_RUN.courant  = 0.5
d.parameters_RUN.peclet   = 0.1

d.parameters_PHYSICS.lcoriol = False 
d.parameters_PHYSICS.igrw_damp = 3

# initialize all velocities to 0 and a low spec. humidity
d.fields[:,:,:].U = 0 | units.m / units.s
d.fields[:,:,:].V = 0 | units.m / units.s
d.fields[:,:,:].W = 0 | units.m / units.s
d.fields[:,:,:].QT = 0.001 | units.kg / units.kg

# add perturbation in temperature - Gaussian bubble at (cx,cy,cz), radius r 
cx,cy,cz,r = 3200|units.m, 3200|units.m, 500|units.m, 500|units.m
d.fields[:,:,:].THL += (0.5 | units.K) * numpy.exp(
    -((d.fields.x-cx)**2 + (d.fields.y-cy)**2 + (d.fields.z-cz)**2)/(2*r**2))

times = numpy.linspace(0, 44, 12) | units.minute # times for snapshots
fig, axes = plt.subplots(3, 4, sharex=True, sharey=True)
extent = (0, d.fields.y[0,-1,0].value_in(units.m),
          0, d.fields.z[0,0,-1].value_in(units.m))
for t,ax in zip(times, axes.flatten()):
    print('Evolving to', t)
    d.evolve_model(t)
    im = ax.imshow(d.fields[16,:,:].THL.value_in(units.K).transpose(), 
              extent=extent, origin='bottom', vmin=292.5, vmax=292.75)
    ax.text(.1, .1, str(t.in_(units.minute)),
             color='w', transform=ax.transAxes)
plt.show()
