"""
A warm air bubble experiment using the DALES Python interface

This script performs the following steps:
* create a Dales object
* set parameters
* create an initial state with a warm air bubble on a constant liquid-water-potential temperature (thl) background
* perform time evolution, periodically storing the model variables
* plot a time series of the model using matplotlib.

"""

from __future__ import division
from __future__ import print_function

from omuse.community.dales.interface import Dales
from omuse.units import units
import matplotlib.pyplot as plt
import numpy

# create Dales object
d = Dales(workdir='bubble', channel_type='sockets', number_of_workers=1)
# add parameter redirection='none' to see DALES diagnostics output


# Set parameters

# Domain size and resolution
d.parameters_DOMAIN.itot = 32  # number of grid cells in x
d.parameters_DOMAIN.jtot = 32  # number of grid cells in y
d.parameters_DOMAIN.xsize = 6400 | units.m
d.parameters_DOMAIN.ysize = 6400 | units.m

# Select advection schemes
d.parameters_DYNAMICS.iadv_mom = 6  # 6th order advection for momentum
d.parameters_DYNAMICS.iadv_thl = 5  # 5th order advection for scalars, less overshoots than 6th order
d.parameters_DYNAMICS.iadv_qt = 5
d.parameters_DYNAMICS.iadv_tke = 5

# turn off randomization of the initial state
d.parameters_RUN.randqt = 0 | units.shu
d.parameters_RUN.randthl = 0 | units.K
d.parameters_RUN.randu = 0 | units.m / units.s

# turn on adaptive time stepping and set more conservative time step limits
d.parameters_RUN.ladaptive = True
d.parameters_RUN.courant = 0.5
d.parameters_RUN.peclet = 0.1

d.parameters_PHYSICS.lcoriol = False
d.parameters_PHYSICS.igrw_damp = 3  # Wind in the sponge layer dampened towards average wind (for symmetric evolution)

# set all velocities to 0
d.fields[:, :, :].U = 0 | units.m / units.s
d.fields[:, :, :].V = 0 | units.m / units.s
d.fields[:, :, :].W = 0 | units.m / units.s

# set a low specific humidity -> no cloud formation
d.fields[:, :, :].QT = 0.001 | units.kg / units.kg


# create a bubble perturbation, given a DALES grid which is used for grid size and coordinates
# if gaussian=True, a gaussian perturbation is generated, with standard deviation r, otherwise a
# constant perturbation is generated inside a sphere of radius r.
#
# r, center are quantities, i.e. numbers with units.

def make_bubble(grid, r, center=None, gaussian=False):
    if center is None:
        ci = ((numpy.array(grid.THL.shape) - 1) * .5)
        ci = (int(ci[0]), int(ci[1]), int(ci[2]))
        print('ci', ci)
        center = (grid[ci].x.value_in(units.m), grid[ci].y.value_in(units.m), grid[ci].z.value_in(units.m))
        print('center', center)
    else:
        center = [c.value_in(units.m) for c in center]

    x = grid[:, 0, 0].x.value_in(units.m)  # fetch coordinate grids once, for speed
    y = grid[0, :, 0].y.value_in(units.m)  # drop the units here for faster calculation below
    z = grid[0, 0, :].z.value_in(units.m)
    r = r.value_in(units.m)

    print(x,y,z,center)


    thl_bubble = numpy.zeros(grid.THL.shape)
    for index, v in numpy.ndenumerate(thl_bubble):
        i, j, k = index
        rx = x[i] - center[0]
        ry = y[j] - center[1]
        rz = z[k] - center[2]
        rr = rx ** 2 + ry ** 2 + rz ** 2
        if gaussian:
            thl_bubble[index] = numpy.exp(-rr / (2 * r ** 2))
        else:
            thl_bubble[index] = 1 if (rr <= r * r) else 0
    return thl_bubble


# create a perturbation: bubble of warm air.
bubble = make_bubble(d.fields, r=500 | units.m, center=(3200 | units.m, 3200 | units.m, 500 | units.m), gaussian=True)
d.fields[:, :, :].THL += 0.5 * bubble | units.K

# evolve model, save a sequence of 3D snapshots at times specified below:
times = numpy.linspace(0, 40, 11) | units.minute

states = []
for t in times:
    state = {}
    print("Evolving to", str(t))
    d.evolve_model(t)

    # save model variables
    state['thl'] = d.fields[:, :, :].THL
    state['qt'] = d.fields[:, :, :].QT
    state['ql'] = d.fields[:, :, :].QL
    state['T'] = d.fields[:, :, :].T
    state['time'] = t
    states.append(state)

C, R = 3, 4  # number of columns and rows in plot

# select field to plot
field, unit = 'thl', units.K  # liquid water potential temperature
# field,unit = 'T', units.K     # temperature
# field,unit = 'qt', units.shu  # total specific humidity
# field,unit = 'ql', units.shu  # specific cloud liquid water


# find range of the variable over all the saved snapshots
vmin = 1e10
vmax = -1e10
for ind in range(len(states)):
    vmin = min(vmin, numpy.amin(states[ind][field].value_in(unit)))
    vmax = max(vmax, numpy.amax(states[ind][field].value_in(unit)))
delta = vmax - vmin
vmax = vmin + delta / 4  # adjust the range for better view of the later stages

# set up the grid extents for proper y and z axes on the plot
e = (0, d.fields.y[0, -1, 0].value_in(units.m), 0, d.fields.z[0, 0, -1].value_in(units.m))  # (left, right, bottom, top)

# plot yz slices at a given x index xi
xi = int(d.fields.THL.shape[0] / 2)  # middle of the system in x

# make a grid of plots of all the snapshots
fig, axes = plt.subplots(R, C, sharex=True, sharey=True)
ind = 0
for j in range(R):
    for i in range(C):
        if ind < len(states):
            f = states[ind][field].value_in(unit)
            time = states[ind]['time']

            im = axes[j, i].imshow(f[xi, :, :].transpose(), origin='bottom', extent=e, vmin=vmin, vmax=vmax)
            axes[j, i].text(.1, .1, str(time.in_(units.minute)), color='w', transform=axes[j, i].transAxes)
            ind += 1
        else:
            # remove un-used axes
            fig.delaxes(axes[j, i])

# set up a color bar
bar_axes = plt.axes((.7, .05, .03, .18))
cbar = plt.colorbar(im, cax=bar_axes)
cbar.set_label('%s (%s)' % (field, unit))

plt.savefig('bubble.png', dpi=200)
plt.savefig('bubble.svg')

plt.show()
