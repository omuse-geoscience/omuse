#!/home/ben/amuse/prerequisites/bin/python

#for plotting
import numpy
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot
from matplotlib.mlab import griddata

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

numpy.random.seed(123451)

def get_lats_lons(p):
    lats= p.nodes.lat.value_in(units.rad).flatten() / numpy.pi*180.0
    lons= p.nodes.lon.value_in(units.rad).flatten() / numpy.pi*180.0
    lats+=numpy.random.random(len(lats))*1.e-4
    lons+=numpy.random.random(len(lons))*1.e-4    
    a=numpy.where(lons>180.0)[0]
    lons[a]=lons[a]-360.0
    return lats,lons

def get_xy():
    N = 500j
    extent = (-180,180.,-90.,90.)
    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    return xs,ys

def myplot(lats, lons, xs, ys, sst):
    sst=sst.flatten()
    sst_=griddata(lons,lats,sst,xs,ys, interp='linear')

    pyplot.clf()
    m = Basemap(projection='kav7',lon_0=0.,llcrnrlat=-80,urcrnrlat=80,\
                llcrnrlon=-180,urcrnrlon=180)
    x, y = m(lons,lats)
    xs,ys= m(xs,ys)

    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    im1 = m.pcolormesh(xs,ys,sst_,cmap=pyplot.cm.jet,latlon=False)

    pyplot.ylim(ys.min(),ys.max())
    pyplot.show()
    pyplot.draw()




from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
from amuse.units import units

#initialize code, print output of code to console
instance = DistributedAmuse(redirection='none')
instance.commit_parameters()

print instance.parameters.webinterface_port

resource = Resource()
resource.name = "Cartesius"
resource.location = "ben@int2.cartesius.surfsara.nl"
resource.scheduler_type = "slurm"
resource.amuse_dir = "/home/ben/amuse/amuse-svn"

instance.resources.add_resource(resource)

pilot = Pilot()
pilot.resource_name="Cartesius"
pilot.queue_name="short" 
pilot.node_count=1
pilot.time= 1|units.hour
pilot.slots_per_node=24
pilot.label="CartesiusNode" 

instance.pilots.add_pilot(pilot)

instance.use_for_all_workers()


from omuse.community.pop.interface import POP
p=POP(channel_type="distributed", redirection="none", number_of_workers=24)
p.change_directory('/home/ben/amuse/amuse-svn/src/omuse/community/pop/')


p.set_horiz_grid_file('data/input/grid/horiz_grid_20010402.ieeer8')
p.set_vert_grid_file('data/input/grid/in_depths.dat')
p.set_topography_file('data/input/grid/topography_20010702.ieeei4')
p.set_ts_file('data/input/restart/r.x1_SAMOC_control.00750101')

#raw_input()


#prepare the plot stuff
pyplot.ion()
lats,lons = get_lats_lons(p)
xs,ys = get_xy()
sst = p.elements.temp.value_in(units.C)
sst = p.elements.shf.value_in(units.W/units.m**2)

print sst[100,120]
print sst.min(),sst.max()
myplot(lats, lons, xs, ys, sst)

#raw_input()


#run the model for 30 days
time = p.get_model_time()
index = 0
while index<30:
    index = index + 1
    time = p.get_model_time()
    tend = time + (1 | units.day)
    p.evolve_model(tend)
    print 'model day ' + str(index) +  ' completed'
    sst = p.elements.temp.value_in(units.C)
    myplot(lats, lons, xs, ys, sst)
    pyplot.savefig("pop-" + str(index) + ".png")







