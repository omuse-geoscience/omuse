#!/usr/bin/python

import numpy
import sys
import os.path

total = len(sys.argv)
if len(sys.argv) != 2:
    sys.exit("Usage: view_grid.py filename")

filename = sys.argv[1]
if not os.path.isfile(filename):
    sys.exit("Error: No such file " + filename)


from netCDF4 import Dataset
fh = Dataset(filename, mode='r')

ny = fh.variables['grid_dims'][0]
nx = fh.variables['grid_dims'][1]

latitudes = fh.variables['grid_center_lat'][:]
longitudes = fh.variables['grid_center_lon'][:]
imask = fh.variables['grid_imask'][:]

lats = latitudes.__array__().reshape(nx,ny)
lons = longitudes.__array__().reshape(nx,ny)
mask = imask.__array__().reshape(nx,ny)

import matplotlib.pyplot as pyplot

f, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
ax1.set_adjustable('box-forced')
ax2.set_adjustable('box-forced')
ax3.set_adjustable('box-forced')
ax4.set_adjustable('box-forced')

#print lats[::-1,:]
ax1.imshow(lats[::-1,:], cmap=pyplot.cm.jet)

#print lons[::-1,:]
ax2.imshow(lons[::-1,:], cmap=pyplot.cm.jet)

#print mask[::-1,:]
ax3.imshow(mask[::-1,:], cmap=pyplot.cm.bone)

if 'grid_area' in fh.variables:
    ax4.imshow(fh.variables['grid_area'][:].reshape(nx,ny)[::-1,:], cmap=pyplot.cm.jet)
else:
    ax4.imshow((numpy.zeros(nx*ny)).reshape(nx,ny), cmap=pyplot.cm.bone)



from mpl_toolkits.mplot3d import axes3d
import numpy

lon = lons
lat = lats

#print min(lon.flatten()), max(lon.flatten())

if (fh.variables['grid_center_lon'].units == "degrees"):
    lon = lon/180*numpy.pi
if (fh.variables['grid_center_lat'].units == "degrees"):
    lat = lat/180*numpy.pi

#compute convenient strides, designed to go up for high resolution grids,
#but not so fast that low resolution grids look poorly
rstride = int(numpy.rint((2*nx/48 + numpy.sqrt(nx))/3))
cstride = int(numpy.rint((2*ny/60 + numpy.sqrt(ny))/3))

f=pyplot.figure(2)
ax = f.add_subplot(111, projection='3d')
X = numpy.cos(lat)*numpy.cos(lon)
Y = numpy.cos(lat)*numpy.sin(lon) 
Z = numpy.sin(lat)
ax.plot_surface(X, Y, Z, color='white', rstride=rstride, cstride=cstride, shade=True)




pyplot.ion()
pyplot.show()

raw_input()
