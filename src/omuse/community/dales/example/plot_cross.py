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
case = "rico"
workdir = "plot2d_" + case
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
    k = 20  # z plane to show in the horizontal plot
    j = 32  # y plane to show in the vertical plot
    fig = plot.figure()
    grid = ImageGrid(fig, 111,  # as in plt.subplot(111)
                     nrows_ncols=(2, n),
                     axes_pad=0.25,
                     share_all=False,  # True,
                     cbar_location="right",
                     cbar_mode="each",  # "single",
                     cbar_size="7%",
                     cbar_pad=0.01,
                     )
    handles = []
    ind = 0
    time_label = fig.text(0.1, 0.1, 'time')
    for f in fields:
        var_name, var_unit = f[0], f[1]
        ax = grid[ind]
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        var_slab = getattr(dales.fields[:, :, k], var_name).value_in(var_unit)
        im = ax.imshow(var_slab.T, origin='lower')
        bar = None  # plot.colorbar()
        grid.cbar_axes[ind].colorbar(im)
        axv = grid[ind + n]
        var_slice = getattr(dales.fields[:, j, :], var_name).value_in(var_unit)
        axv.set_xlabel("x")
        axv.set_ylabel("z")
        imv = axv.imshow(var_slice.T, origin='lower')
        barv = None  # plot.colorbar()
        grid.cbar_axes[ind + n].colorbar(imv)
        handles.append((ax, im, bar, axv, imv, barv))
        ind += 1

    # plot.tight_layout()
    # draw a check box for forcing
    rax = plot.axes([0.25, 0.05, 0.15, 0.08])
    check = CheckButtons(rax, ('Forcing',), (False,))
    check.on_clicked(checkbtn_clicked)

    def set_min_max(img, fld):
        vmax = numpy.max(fld)
        vmin = numpy.min(fld)
        img.set_clim(vmin, vmax)

    def update_img(img, fld):
        img.set_data(fld.T)
        set_min_max(img, fld)

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
            hslab = getattr(dales.fields[:, :, k], name).value_in(unit)
            update_img(h[1], hslab)
            vslice = getattr(dales.fields[:, j, :], name).value_in(unit)
            update_img(h[4], vslice)

    a = movie.FuncAnimation(fig, update, frames=num_frames, repeat=False)
    plot.show()


def main(args):
    dales = Dales(case=case, workdir=workdir, number_of_workers=dales_procs, redirection="none")
    dales.commit_parameters()
    plot_movie(dales, 1000)
    cleanup()


if __name__ == "__main__":
    main(sys.argv[1:])
