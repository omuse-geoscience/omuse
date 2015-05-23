#!/home/ben/amuse/prerequisites/bin/python

#for plotting
import numpy
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot
from matplotlib.mlab import griddata



#enable logging
import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)



def lons_and_lats():
    raw=numpy.fromfile("/home/ben/x1_Files/grid/horiz_grid_20010402.ieeer8").newbyteorder('S')
    grid=raw.reshape((7,384,320))
    lats=grid[0,:,:]/numpy.pi*180
    lons=grid[1,:,:]/numpy.pi*180
    return lats,lons


#plot it
def plot(lats, lons, sst):

    lats=lats.flatten()
    lats+=numpy.random.random(len(lats))*1.e-4
    lons=lons.flatten()
    lons+=numpy.random.random(len(lons))*1.e-4
    sst=sst.flatten()

    N = 500j
    extent = (-180,180.,-90.,90.)

    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

    a=numpy.where(lons>180)[0]
    lons[a]=lons[a]-360

    sst_=griddata(lons,lats,sst,xs,ys,interp='linear')

    pyplot.clf()
    m = Basemap(projection='kav7',lon_0=0.,llcrnrlat=-80,urcrnrlat=80,\
                llcrnrlon=-180,urcrnrlon=180)
    x, y = m(lons,lats)
    print xs.min(),xs.max()
    xs,ys= m(xs,ys)
    print xs.min(),xs.max()

    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')

    im1 = m.pcolormesh(xs,ys,sst_,cmap=pyplot.cm.jet,latlon=False)

    pyplot.ylim(ys.min(),ys.max())
  #  pyplot.savefig("pop.png")
    pyplot.show()

    pyplot.draw()






#from amuse.community.pop.interface import POPInterface
from omuse.community.pop.interface import POPInterface

p=POPInterface(channel_type="sockets",redirection="none")
p.initialize_code()

#prepare index arrays to fetch SST
size=p.get_domain_size()
i_array = []
j_array = []
for j in range(1,size['ny']+1):
    for i in range(1,size['nx']+1):
         i_array.append(i)
         j_array.append(j)

i=i_array
j=j_array



#prepare the plot stuff
pyplot.ion()

#print i

#print j


#actually fetch and plot SST

import time

start=time.clock()
sst = p.get_element_surface_state(i,j)['temp']
end=time.clock()
print 'gather state took:' + str(end-start) + 'secs.'


#print sst

start=time.clock()
lat_lon = p.get_node_position(i,j)
end=time.clock()
print 'reduce position took:' + str(end-start) + 'secs.'


lats=lat_lon['lat']/numpy.pi*180
lons=lat_lon['lon']/numpy.pi*180


#lats,lons=lons_and_lats()

plot(lats, lons, sst)
pyplot.savefig("pop.png")




time = p.get_model_time()['time']
dt = p.get_timestep()['dt']

index = 0
tend = time + dt

while True:
    time = p.get_model_time()['time']
    tend = time + dt
    p.evolve_model(tend)
    print 'timestep completed, time=', tend
    sst = p.get_element_surface_state(i,j)['temp']
    index = index + 1
    plot(lats, lons, sst)
    pyplot.savefig("pop-" + str(index) + ".png")







