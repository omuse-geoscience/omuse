#!/usr/bin/env python

from omuse.community.pop.interface import POP
from amuse.units import units
import numpy

p=POP(channel_type="sockets",redirection="none",number_of_workers=8)

#set the grid we want to use
p.set_horiz_grid_file('data/input/grid/horiz_grid_20010402.ieeer8')
p.set_vert_grid_file('data/input/grid/in_depths.dat')
p.set_topography_file('data/input/grid/topography_20010702.ieeei4')

#set the restart file
p.set_ts_file('data/input/restart/r.x1_SAMOC_control.00750101')

#setup the forcing
p.set_shf_monthly_file('data/input/shf_monthly/shf.normal_year+flux.mon')
p.set_sfwf_monthly_file('data/input/sfwf/sfwf_phc0-50_ncarp_r46_flux.mon')
p.set_ws_monthly_file('data/input/ws_monthly/ws.1958-2000.mon')



#retrieve sst values
sst = p.elements.temp.value_in(units.C)

#create a plot using sst values
#from pop_plot import *
#pop_plot(p, sst)

#from matplotlib import pyplot
#pyplot.imshow(sst)
#pyplot.show()


#go interactive
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
import code
code.interact(local=dict(globals(), **locals()) )

