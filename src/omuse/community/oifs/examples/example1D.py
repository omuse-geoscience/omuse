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

oifs_procs = 1
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
    sp = fig.add_subplot(1,1,1)
    sp.set_xlabel("T(K)")
    sp.set_ylabel("z")
    height = numpy.arange(0,oifs.ktot)
    def update(i):
        t = oifs.get_model_time()
        oifs.evolve_model(t + (steps/frames)*oifs.get_timestep())
        z = oifs.get_profile_field("T",oifs.itot/2).value_in(units.K)
        print z
        sp.plot(z,height)
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
