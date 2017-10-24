#!/usr/bin/env python

from amuse.community import *
from omuse.community.dales.interface import Dales

import numpy
import matplotlib.pyplot as plot
import matplotlib.animation as movie
import os
import sys
import shutil

dales_procs=1
case="example"
dales=None
copied_files=[]

def init():
    global dales
    srcdir=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","dales-repo","cases",case))
    destdir=os.path.dirname(os.path.abspath(__file__))
    files=os.listdir(srcdir)
    for f in files:
        src=os.path.join(srcdir,f)
        if(os.path.isfile(src)):
            shutil.copy(src,destdir)
            copied_files.append(os.path.join(destdir,f))
    dales = Dales(number_of_workers=dales_procs,redirection="none")

def cleanup():
    dales.cleanup_code()
    dales.stop()
    for f in copied_files:
        os.remove(f)

def plotmovie(n,tstart):
    depth=numpy.arange(1,96)
    fig=plot.figure()
    ax=fig.add_subplot(1,1,1)
    ax.set_ylabel("layer")
    ax.set_xlabel("T(K)")
    def update(i):
        dales.evolve_model(tstart + (i | units.s))
        profile=dales.get_profile_T_(depth).value_in(units.K)
        ax.clear()
        ax.plot(profile[::-1],depth)
    a=movie.FuncAnimation(fig,update,frames=n,repeat=False)
    plot.show()


def main(args):
    init()
    dales.commit_parameters()
    tim=dales.get_model_time()
    dales.commit_grid()
    plotmovie(30,tim)
    cleanup()

if __name__=="__main__":
    main(sys.argv[1:])
