#!/usr/bin/env python

# Run DALES using the Omuse interface, retrieve fields and plot them using Matplotlib
# 2016
#
# Draw 2D slices of selected fields.
# For each field, a horizontal and a vertical slice is drawn. 
# The color maps are independent for each slice

from amuse.community import *
from omuse.community.dales.interface import Dales

import numpy
import matplotlib.pyplot as plot
import matplotlib.animation as movie
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.widgets import CheckButtons 
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
    print 'Copying input files'
    for f in files:
        src=os.path.join(srcdir,f)
        if(os.path.isfile(src)):
            print src, '->', destdir
            shutil.copy(src,destdir)
            copied_files.append(os.path.join(destdir,f))
    dales = Dales(number_of_workers=dales_procs,redirection="none")

def cleanup():
    dales.cleanup_code()
    dales.stop()
    print 'Removing input files:', copied_files
    for f in copied_files:
        os.remove(f)

forcingFlag = False
        
def checkButtonClicked(label):
    global forcingFlag
    forcingFlag = not forcingFlag
    print('Click', label, forcingFlag)

        
# a list of the fields to draw
plots = [('THL', units.K),
         ('T', units.K),
         ('QT', None),
         ('QL', None),
         ('E12', None),
         
    ]

def plotmovie(n):
    N = len(plots) # number of fields to draw
    k = 20 # z plane to show in the horizontal plot
    j = 32 # y plane to show in the vertical plot
    
    fig=plot.figure()
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(2,N),
                 axes_pad=0.25,
                 share_all=False, #True,
                 cbar_location="right",
                 cbar_mode= "each", # "single",
                 cbar_size="7%",
                 cbar_pad=0.01,
                 )
    
    handles = []        
    ind = 0
    timelabel = fig.text(0.1, 0.1, 'time')
    for p in plots:
        #ax=fig.add_subplot(2,N,ind) #rows, columns, current
        ax = grid[ind]
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        ax.set_title(p[0])
        print p
        field = dales.get_field(p[0], kmin=k, kmax=k+1)[:,:,0]
        if p[1] != None: # drop units
            field = field.value_in(p[1])
            
        im = ax.imshow(field.T, origin='lower')
        bar = None #plot.colorbar()
        grid.cbar_axes[ind].colorbar(im)
        
        #axv=fig.add_subplot(2,N,N+ind)
        print ind
        axv = grid[ind+N]
        field  = dales.get_field(p[0], jmin=j, jmax=j+1)[:,0,:]
        if p[1] != None: # drop units
            field = field.value_in(p[1])
            
        axv.set_xlabel("x")
        axv.set_ylabel("z")
        imv = axv.imshow(field.T, origin='lower')
        barv = None #plot.colorbar()
        grid.cbar_axes[ind+N].colorbar(imv)
        handles.append((ax, im, bar, axv, imv, barv))
        ind += 1

    #plot.tight_layout()

    #draw a check box for forcing
    rax = plot.axes([0.25, 0.05, 0.15, 0.08])
    check = CheckButtons(rax, ('Forcing',), (False,))
    check.on_clicked(checkButtonClicked)

    
    def setMinMax(im, field):
        vmax     = numpy.max(field)
        vmin     = numpy.min(field)
        im.set_clim(vmin, vmax)

    def updateIm(im, field):
        im.set_data(field.T)
        setMinMax(im, field)
        
    def update(i):
        t=dales.get_model_time()

        tendency = numpy.zeros(dales.k)
        if forcingFlag:
            print 'Applying forcing'
            for i in range(-15,15):
                tendency[30+i] =  numpy.exp(-(i/5.0)**2)
        #dales.set_tendency_U(tendency)
        dales.set_tendency_V(tendency * .2)
        dales.set_tendency_QT(tendency* .000002)
        
        dales.evolve_model(t + (30 | units.s))

        t=dales.get_model_time()
        print "updating with arg",i, 'time = ', t
        timelabel.set_text('time %.2f'%t.value_in(units.s))
        
        for p,h in zip(plots,handles):
            field = dales.get_field(p[0], kmin=k, kmax=k+1)[:,:,0]
            if p[1] != None: # drop units
                field = field.value_in(p[1])
            updateIm(h[1], field)
            
            field  = dales.get_field(p[0], jmin=j, jmax=j+1)[:,0,:]
            if p[1] != None: # drop units
                field = field.value_in(p[1])
            updateIm(h[4], field)
                



        
    a=movie.FuncAnimation(fig,update,frames=n,repeat=False)
    plot.show()


def main(args):
    init()
    dales.commit_parameters()
    tim=dales.get_model_time()
    dales.commit_grid()
    plotmovie(1000)
    cleanup()

if __name__=="__main__":
    main(sys.argv[1:])
