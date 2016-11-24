#!/usr/bin/env python

from amuse.community import *
from omuse.community.oifs.interface import OpenIFS

import numpy
import matplotlib.pyplot as plot
import matplotlib.animation as movie
from mpl_toolkits.basemap import Basemap
import os
import sys
import shutil

oifs_procs = 4
oifs = None

def init():
    global oifs
    oifs = OpenIFS(number_of_workers=oifs_procs,redirection="none")
    oifs.initialize_code()


def plotmovie(frames,steps):
    size = oifs.itot
    x = oifs.longitudes.value_in(units.deg)
    y = oifs.latitudes.value_in(units.deg)
    fig = plot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    m = Basemap(lat_0=-90,lat_1=90,lon_0=180,resolution='c')
    m.drawcoastlines()
    def update(i):
        t = oifs.get_model_time()
        oifs.evolve_model(t + (steps/frames)*oifs.get_timestep())
        z = oifs.get_layer_field("SH",oifs.ktot-1)
        plot.scatter(x,y,c = z,s = 150,cmap = plot.cm.coolwarm,edgecolors = "none")
    a = movie.FuncAnimation(fig,update,frames=frames,repeat=False)
    plot.show()


def main(args):
    init()
    oifs.commit_parameters()
    tim=oifs.get_model_time()
    oifs.commit_grid()
    plotmovie(frames=10,steps=10.)
    oifs.cleanup_code()

if __name__=="__main__":
    main(sys.argv[1:])
