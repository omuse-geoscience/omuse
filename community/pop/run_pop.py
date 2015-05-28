#!/home/ben/amuse/prerequisites/bin/python

#for plotting
import numpy
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot
from matplotlib.mlab import griddata


def get_lats_lons(p):
    lats= p.nodes.lat.value_in(units.deg).flatten() / numpy.pi*180.0
    lons= p.nodes.lon.value_in(units.deg).flatten() / numpy.pi*180.0
    a=numpy.where(lons>180)[0]
    lons[a]=lons[a]-360
    return lats,lons

def get_xy():
    N = 500j
    extent = (-180,180.,-90.,90.)
    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    return xs,ys

def myplot(lats, lons, xs, ys, sst):
    sst=sst.flatten()
    sst_=griddata(lons,lats,sst,xs,ys,interp='linear')

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


from omuse.community.pop.interface import POP
from amuse.units import units

p=POP(channel_type="sockets",redirection="none",number_of_workers=8)

p.set_horiz_grid_file('/home/ben/x1_Files/grid/horiz_grid_20010402.ieeer8')
p.set_vert_grid_file('/home/ben/x1_Files/grid/in_depths.dat')
p.set_topography_file('/home/ben/x1_Files/grid/topography_20010702.ieeei4')
p.set_ts_file('/home/ben/x1_Files/restart/r.x1_SAMOC_control.00750101')



#prepare the plot stuff
pyplot.ion()
lats,lons = get_lats_lons(p)

xs,ys = get_xy()
sst = p.elements.temp.value_in(units.C)

myplot(lats, lons, xs, ys, sst)


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







