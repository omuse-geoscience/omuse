#!/usr/bin/env python

from omuse.community.pop.read_grid import pop_grid_reader
import numpy

class pop_grid_converter:

    def __init__(self, filename, nx, ny):

        r=pop_grid_reader(filename, nx, ny)

        #build the t-grid the same way SCRIP does
        #ensure all values are in [0,2Pi]
        u_lat = r.lat
        u_lon = r.lon
        for j in range(ny):
            for i in range(nx):
                while u_lon[j,i] > 2.0*numpy.pi:
                    u_lon[j,i] -= 2.0*numpy.pi
                while u_lon[j,i] < 0.0:
                    u_lon[j,i] += 2.0*numpy.pi

        #the u points are considered the ne corners of t-cells
        t_ne_lat = u_lat
        t_ne_lon = u_lon

        #nw corners are obtained from rotating the grid
        t_nw_lat = numpy.take(t_ne_lat, [nx-1] + range(0,nx-1), axis=1)
        t_nw_lon = numpy.take(t_ne_lon, [nx-1] + range(0,nx-1), axis=1)

        tiny = 1.0e-14
        t_se_lat = numpy.zeros([ny,nx], dtype=numpy.double)
        t_sw_lat = numpy.zeros([ny,nx], dtype=numpy.double)
        t_se_lon = numpy.zeros([ny,nx], dtype=numpy.double)
        t_sw_lon = numpy.zeros([ny,nx], dtype=numpy.double)

        #copy ne and nw corners to se and sw corners
        t_se_lat[1:,:] = t_ne_lat[:ny-1,:]
        t_se_lon[1:,:] = t_ne_lon[:ny-1,:]
        t_sw_lat[1:,:] = t_nw_lat[:ny-1,:]
        t_sw_lon[1:,:] = t_nw_lon[:ny-1,:]

        #mock up the lower row boundaries
        t_se_lat[0,:] = (-numpy.pi*0.5) + tiny
        t_se_lon[0,:] = t_ne_lon[0,:]
        t_sw_lat[0,:] = (-numpy.pi*0.5) + tiny
        t_sw_lon[0,:] = t_nw_lon[0,:]

        #correct for 0,2pi longitude crossings the same way as SCRIP does
        for j in range(ny):
            for i in range(nx):
                if t_sw_lon[j,i] > 2.0*numpy.pi:
                    t_sw_lon[j,i] -= 2.0*numpy.pi
                if t_sw_lon[j,i] < 0.0:
                    t_sw_lon[j,i] += 2.0*numpy.pi

                tmplon = t_se_lon[j,i] - t_sw_lon[j,i]
                if tmplon < -1.5*numpy.pi:
                    t_se_lon[j,i] += 2.0*numpy.pi
                if tmplon > 1.5*numpy.pi:
                    t_se_lon[j,i] -= 2.0*numpy.pi

                tmplon = t_ne_lon[j,i] - t_se_lon[j,i]
                if tmplon < -1.5*numpy.pi:
                    t_ne_lon[j,i] += 2.0*numpy.pi
                if tmplon > 1.5*numpy.pi:
                    t_ne_lon[j,i] -= 2.0*numpy.pi

                tmplon = t_nw_lon[j,i] - t_ne_lon[j,i]
                if tmplon < -1.5*numpy.pi:
                    t_nw_lon[j,i] += 2.0*numpy.pi
                if tmplon > 1.5*numpy.pi:
                    t_nw_lon[j,i] -= 2.0*numpy.pi

        #compute ocean cell centers by averaging corner values
        t_lat = (t_sw_lat + t_se_lat + t_ne_lat + t_nw_lat) / 4.0
        t_lon = (t_sw_lon + t_se_lon + t_ne_lon + t_nw_lon) / 4.0

        for j in range(ny):
            for i in range(nx):
                if t_lon[j,i] > 2.0*numpy.pi:
                    t_lon[j,i] -= 2.0*numpy.pi
                if t_lon[j,i] < 0.0:
                    t_lon[j,i] += 2.0*numpy.pi

        #store info for reading by the remapper
        self.size = ny*nx
        self.dims = [nx, ny]
        self.corners = 4

        self.t_cell_center_lat = t_lat.flatten()
        self.t_cell_center_lon = t_lon.flatten()

        print self.t_cell_center_lat
        print self.t_cell_center_lon

        self.t_cell_corner_lat = numpy.dstack((t_sw_lat, t_se_lat, t_ne_lat, t_nw_lat)).flatten()
        self.t_cell_corner_lon = numpy.dstack((t_sw_lon, t_se_lon, t_ne_lon, t_nw_lon)).flatten()

        print self.t_cell_corner_lat
        print self.t_cell_corner_lon





if __name__ == "__main__":

    #nx=320
    #ny=384
    #c = pop_grid_converter("../pop/data/input/grid/horiz_grid_20010402.ieeer8", nx, ny)

    nx=3600
    ny=2400
    c = pop_grid_converter("/home/ben/Downloads/SCRIP/trunk/SCRIP/grids/input/grid.3600x2400.dat", nx, ny)


    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import axes3d

    lon = c.t_cell_corner_lon[3::4]
    lat = c.t_cell_corner_lat[3::4]

    #compute convenient strides, designed to go up for high resolution grids,
    #but not so fast that low resolution grids look poorly
    rstride = int(numpy.rint((2*nx/48 + numpy.sqrt(nx))/3))
    cstride = int(numpy.rint((2*ny/60 + numpy.sqrt(ny))/3))

    #"""
    f=pyplot.figure()
    ax = f.add_subplot(111, projection='3d')
    X = numpy.cos(lat)*numpy.cos(lon)
    Y = numpy.cos(lat)*numpy.sin(lon)
    Z = numpy.sin(lat)

    ax.plot_surface(X.reshape(ny,nx), Y.reshape(ny,nx), Z.reshape(ny,nx), color="white", shade=True, rstride=rstride, cstride=cstride)
    #"""

    #pyplot.imshow(lon.reshape(ny, nx), cmap=pyplot.cm.jet)

    pyplot.show()
    raw_input()
