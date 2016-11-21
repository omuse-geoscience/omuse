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

def plotmovie(n):
    size=64
    x=numpy.zeros(size*size)
    y=numpy.zeros(size*size)
    k=0
    for i in range(1,size+1):
        for j in range(1,size+1):
            x[k]=i
            y[k]=j
            k+=1
    fig=plot.figure()
    ax=fig.add_subplot(1,1,1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    def update(i):
        print "updating with arg",i
        t=dales.get_model_time()
        dales.evolve_model(t + (1 | units.s))
        field=dales.get_layer_field(60,x,y).value_in(units.K).reshape(size,size)
        plot.imshow(field)
    a=movie.FuncAnimation(fig,update,frames=n,repeat=False)
    plot.show()


def main(args):
    init()
    dales.commit_parameters()
    tim=dales.get_model_time()
    dales.commit_grid()
    plotmovie(30)
    cleanup()

if __name__=="__main__":
    main(sys.argv[1:])
