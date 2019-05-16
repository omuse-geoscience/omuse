#!/usr/bin/env python

from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
from amuse.units import units

import numpy
from matplotlib import pyplot

instance = DistributedAmuse(redirection="none")
instance.parameters.debug = False
instance.parameters.worker_queue_timeout=1 | units.hour
instance.commit_parameters()

print instance.parameters.webinterface_port

resource = Resource()
resource.name = "DAS-5"
resource.location = "bwn200@fs0.das5.cs.vu.nl"
resource.scheduler_type = "slurm"
resource.amuse_dir = "/home/bwn200/amuse/amuse"
resource.tmp_dir = "/home/bwn200/tmp"

instance.resources.add_resource(resource)

pilot = Pilot()
pilot.resource_name="DAS-5"
pilot.queue_name="defq" 
pilot.node_count=56
pilot.time= 1|units.hour
pilot.slots_per_node=16
pilot.label="DAS-5-Pilot"

instance.pilots.add_pilot(pilot)
instance.use_for_all_workers()


from omuse.community.pop.interface import POP
p=POP(channel_type="distributed", redirection="none", mode='3600x2400x42', number_of_workers=896, max_message_length=1000000)

#set grid info
p.set_horiz_grid_file('/var/scratch/bwn200/pop/input/grid/grid.3600x2400.fob.da')
p.set_vert_grid_file('/var/scratch/bwn200/pop/input/grid/in_depths.42.dat')
p.set_topography_file('/var/scratch/bwn200/pop/input/grid/kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black')
p.set_bottom_cell_file('/var/scratch/bwn200/pop/input/grid/dzbc_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black')
p.set_ts_file('/var/scratch/bwn200/pop/input/r.t0.1_42l_greenland.01150501')

#setup forcing files
p.set_shf_monthly_file('/var/scratch/bwn200/pop/input/forcing/shf.NY+H+f.mon')
p.set_sfwf_monthly_file('/var/scratch/bwn200/pop/input/forcing/sfwf.C+r+g8+f.mon')
p.set_ws_monthly_file('/var/scratch/bwn200/pop/input/forcing/ws.o_n_avg.mon')

#setup output files
#p.set_tavg_option('nday')
#p.set_tavg_freq_option(1)
#p.set_tavg_file('/var/scratch/bwn200/pop/output/tavg/t')



sst = p.elements.temp.value_in(units.C)

pyplot.imshow(numpy.swapaxes(sst,0,1), origin="lower", cmap=pyplot.cm.jet)
pyplot.show()

raw_input()




#go interactive
import readline
import rlcompleter
readline.parse_and_bind("tab: complete")
import code
code.interact(local=dict(globals(), **locals()) )



