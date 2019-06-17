#!/usr/bin/env python

# Run DALES using the Omuse interface.
# Minimal test case
#  * initialize the model
#  * evolve it 30 s
#  * get the Qt profile
#  * quit

from amuse.community import *
from omuse.community.dales.interface import Dales

#~ import logging
#~ logging.basicConfig(level=logging.DEBUG)
#~ logging.getLogger("code").setLevel(logging.DEBUG)

import numpy
import os
import sys
import shutil
import datetime

dales_procs=2

inputdir = os.path.abspath(os.path.join(os.path.dirname(__file__),"dales-repo/cases/bomex"))
workdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dales-work')

copied_files=[]

def init(inputdir, workdir):

    try:
        os.makedirs(workdir) # Create workdir
    except:                  # in case the directory already exists
        print('Working directory %s already exists.'%workdir)
        #sys.exit(1)
        
    files=os.listdir(inputdir)

    print ('Copying input files')
    for f in files:
        src=os.path.join(inputdir,f)
        if(os.path.isfile(src)):
            print (src, '->', workdir)
            shutil.copy(src,workdir)
            copied_files.append(os.path.join(workdir,f))

    os.chdir(workdir)
    dales = Dales(number_of_workers=dales_procs, channel_type='sockets')
    #~ dales = Dales(number_of_workers=dales_procs,redirection="none", channel_type='sockets')
    #~ dales = Dales(number_of_workers=dales_procs, channel_type='sockets', debugger='gdb')
    return dales


def main(args):

    dales=init(inputdir, workdir)
    dales.parameters.starttime=datetime.datetime(2013,12,30,23,55,10)
    dales.parameters.evolve_to_exact_time=True

    print('Initialization done.')
        
    t=dales.model_time

    print("starting time", t)

    dales.evolve_model(t + (30 | units.s))
    t=dales.model_time
    print('Dales evolved to ', t)

    qt = dales.get_profile_QT()
    print('Qt = ', qt)
    
    dales.stop()
    print('Done.')

    
if __name__=="__main__":
    main(sys.argv[1:])
