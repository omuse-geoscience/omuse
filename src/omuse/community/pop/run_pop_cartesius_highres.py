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
    nodes = p.nodes
    lats= nodes.lat.value_in(units.rad).flatten() / numpy.pi*180.0
    lons= nodes.lon.value_in(units.rad).flatten() / numpy.pi*180.0
   # lats+=numpy.random.random(len(lats))*1.e-4
   # lons+=numpy.random.random(len(lons))*1.e-4    
    a=numpy.where(lons>180.0)[0]
    lons[a]=lons[a]-360.0
    return lats,lons

def get_xy(p):
    size = p.get_domain_size()
    N = max(size)*1j
   # N = 500j
    extent = (-180,180.,-90.,90.)
    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    return xs,ys


#true_plot is used for visualizing how the grid sits in memory
#without doing any grid transformations
def true_plot(p, lats, lons, value):
    extent = (-180,180.,-90.,90.)
    size = p.get_domain_size()
    xs,ys = numpy.mgrid[extent[0]:extent[1]:size[0]*1j, extent[2]:extent[3]:size[1]*1j] 

    pyplot.clf()
    m = Basemap(projection='kav7',lon_0=0.,llcrnrlat=-80,urcrnrlat=80,\
                llcrnrlon=-180,urcrnrlon=180)
    x, y = m(lons,lats)
    xs,ys= m(xs,ys)
    
    m.drawmapboundary(fill_color='#99ffff')
    im1 = m.pcolormesh(xs,ys,value,cmap=pyplot.cm.jet,latlon=False)

    pyplot.ylim(ys.min(),ys.max())
    pyplot.show()
    pyplot.draw()


#plot is nice for visualizing the data, it uses griddata() to perform the
#transformation and interpolation from any grid to something that can be
#visualized. It is unclear how accurate/conserving/sphere-aware the trans-
#formation is at the moment.
def myplot(lats, lons, xs, ys, sst):
    sst=sst.ravel()

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

"""
from scipy.interpolate import SmoothSphereBivariateSpline
def remap_to_latlon_grid(p, value):
    lats= p.nodes.lat.value_in(units.rad)
    lons= p.nodes.lon.value_in(units.rad)

    lut = SmoothSphereBivariateSpline(lats.ravel(), lons.ravel(), data.ravel())
    size = p.get_domain_size()
    nx_global = size[0]
    ny_global = size[1]
    rect_lats = numpy.linspace(0., numpy.pi, ny_global)
    rect_lons = numpy.linspace(0., 2 * numpy.pi, nx_global)
    rect_data = lut(fine_lats, fine_lons)

    return rect_data
"""



from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
from amuse.units import units

#initialize code, print output of code to console
instance = DistributedAmuse(redirection='none')
#instance.parameters.debug = True
instance.parameters.worker_queue_timeout=1 | units.hour
instance.commit_parameters()

print instance.parameters.webinterface_port

resource = Resource()
resource.name = "Cartesius"
resource.location = "ben@int1.cartesius.surfsara.nl"
resource.scheduler_type = "slurm"
resource.amuse_dir = "/home/ben/amuse/amuse-svn"
resource.tmp_dir = "/scratch-local/ben/distributed-amuse/"

instance.resources.add_resource(resource)

pilot = Pilot()
pilot.resource_name="Cartesius"
pilot.queue_name="short" 
pilot.node_count=50
pilot.time= 1|units.hour
pilot.slots_per_node=12
pilot.label="CartesiusNode" 

instance.pilots.add_pilot(pilot)

instance.use_for_all_workers()


from omuse.community.pop.interface import POP
p=POP(channel_type="distributed", redirection="none", number_of_workers=600)
p.change_directory('/home/ben/amuse/amuse-svn/src/omuse/community/pop/')
p.set_namelist_filename('pop_in_highres')

p.set_horiz_grid_file('/home/ben/pop/input/grid/grid.3600x2400.fob.da')
p.set_vert_grid_file('/home/ben/pop/input/grid/in_depths.42.dat')
p.set_topography_file('/home/ben/pop/input/grid/kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black')
p.set_bottom_cell_file('/home/ben/pop/input/grid/dzbc_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black')
p.set_ts_file('/home/ben/pop/input/r.t0.1_42l_greenland.01150501')

p.set_monthly_shf_file('/home/ben/pop/input/forcing/shf.NY+H+f.mon')
p.set_monthly_sfwf_file('/home/ben/pop/input/forcing/sfwf.C+r+g8+f.mon')
p.set_monthly_ws_file('/home/ben/pop/input/forcing/ws.o_n_avg.mon')



#raw_input()


#prepare the plot stuff
pyplot.ion()
sst = p.elements.temp.value_in(units.C)
#sst = p.elements.shf.value_in(units.W/units.m**2)

print sst[100,120]
print sst.min(),sst.max()

#lats,lons = get_lats_lons(p)
#xs,ys = get_xy(p)
#myplot(lats, lons, xs, ys, sst)
#true_plot(p, lats, lons, sst)

pyplot.imshow(numpy.swapaxes(sst,0,1)[::-1,:], cmap=pyplot.cm.jet)

raw_input()

#see if it crashes
temp3d = p.elements3d.temp[:,:,1]

raw_input()



"""

#run the model for 30 days
time = p.get_model_time()
index = 0
while True:
    index = index + 1
    time = p.get_model_time()
    tend = time + p.get_timestep_next()
    p.evolve_model(tend)
    print 'model step ' + str(index) +  ' completed'
    sst = p.elements.temp.value_in(units.C)
    print sst.min(),sst.max()
    pyplot.imshow(numpy.swapaxes(sst,0,1)[::-1,:], cmap=pyplot.cm.jet)
    #myplot(lats, lons, xs, ys, sst)

    raw_input()

"""




