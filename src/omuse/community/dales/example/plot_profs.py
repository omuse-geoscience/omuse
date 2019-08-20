#!/usr/bin/env python

# Run DALES using the Omuse interface, retrieve fields and plot them using Matplotlib
# 2016
#
# Draw 2D slices of selected fields.
# For each field, a horizontal and a vertical slice is drawn. 
# The color maps are independent for each slice

from omuse.units import units
from omuse.community.dales.interface import Dales

import numpy
import matplotlib.pyplot as plot
import matplotlib.animation as movie
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.widgets import CheckButtons
import os
import sys

dales_procs = 1
case = "bomex"
workdir = "plot1d_" + case
forcingFlag = False
dt = 30 | units.s
fields = [('THL', units.K),
          ('T', units.K),
          ('QT', units.mfu),
          ('QL', units.mfu),
          ('E12', units.m / units.s)]


def cleanup():
    if os.path.isdir(workdir):
        for f in os.listdir(workdir):
            os.remove(os.path.join(workdir, f))
        os.rmdir(workdir)
    elif os.path.isfile(workdir):
        os.remove(workdir)


def checkbtn_clicked(label):
    global forcingFlag
    forcingFlag = not forcingFlag
    print('Click', label, forcingFlag)


def plot_movie(dales, num_frames):
    n = len(fields)  # number of fields to draw
    fig, axes = plot.subplots(1, n, sharey=True)
    handles = []
    time_label = fig.text(0.1, 0.1, 'time')
    zf = dales.get_zf().value_in(units.m)
    ind = 0
    colors = "bgrcmyk"
    for ax in axes:
        f = fields[ind]
        var_name, var_unit = f[0], f[1]
        ax = axes[ind]
        ax.set_xlabel(var_name)
        ax.set_ylabel("z")
        var_prof = getattr(dales.profiles, var_name).value_in(var_unit)
        line = ax.plot(var_prof[:], zf[:], '-', color=colors[ind], linewidth=2)
        handles.append((ax, line))
        ind += 1

    # plot.tight_layout()
    # draw a check box for forcing
    rax = plot.axes([0.25, 0.05, 0.15, 0.08])
    check = CheckButtons(rax, ('Forcing',), (False,))
    check.on_clicked(checkbtn_clicked)

    def update_plot(lines, prof):
        lines[0].set_xdata(prof[:])

    def update(i):
        if forcingFlag:
            tendency = numpy.zeros(dales.get_ktot())
            print 'Applying forcing'
            for i in range(-15, 15):
                tendency[30 + i] = numpy.exp(-(i / 5.0) ** 2)
            dales.set_tendency_V(tendency * .2 | units.m / units.s**2)
            dales.set_tendency_QT(tendency * .000002 | units.mfu / units.s)

        t = dales.get_model_time()
        dales.evolve_model(t + dt)

        t = dales.get_model_time()
        print "updating with arg", i, 'time = ', t
        time_label.set_text('time %.2f' % t.value_in(units.s))

        for f, h in zip(fields, handles):
            name, unit = f[0], f[1]
            prof = getattr(dales.profiles, name).value_in(unit)
            update_plot(h[1], prof)

    a = movie.FuncAnimation(fig, update, frames=num_frames, repeat=False)
    plot.show()


def main(args):
    dales = Dales(case=case, workdir=workdir, number_of_workers=dales_procs, redirection="none")
    dales.commit_parameters()
    plot_movie(dales, 1000)
    cleanup()


if __name__ == "__main__":
    main(sys.argv[1:])
