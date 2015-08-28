#!/usr/bin/python

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

src_nx = fh.variables['src_grid_dims'][0]
src_ny = fh.variables['src_grid_dims'][1]
dst_nx = fh.variables['dst_grid_dims'][0]
dst_ny = fh.variables['dst_grid_dims'][1]
num_links = fh.dimensions['num_links']
remap_matrix = fh.variables['remap_matrix'][:]

#print min(fh.variables['src_address']), max(fh.variables['src_address'])
#print min(fh.variables['dst_address']), max(fh.variables['dst_address'])
#print 'src_ny=', src_ny, ' src_nx=', src_nx

import numpy
src_weights = numpy.zeros(src_nx*src_ny).reshape(src_ny, src_nx)
dst_weights = numpy.zeros(dst_nx*dst_ny).reshape(dst_ny, dst_nx)
for (i,j), weight in numpy.ndenumerate(remap_matrix):
    src_i = fh.variables['src_address'][i] -1 #correct for Fortran indexing starting at 1
    src_y = src_i/src_nx
    src_x = src_i-(src_y*src_nx)
#    print src_y,src_x,src_i,j,i
    src_weights[src_y,src_x] += remap_matrix[i,j]
    dst_i = fh.variables['dst_address'][i] -1 #correct for Fortran indexing starting at 1
    dst_y = dst_i/dst_nx
    dst_x = dst_i-(dst_y*dst_nx)
#    print dst_y,dst_x,dst_i,j,i
    dst_weights[dst_y,dst_x] += remap_matrix[i,j]

src_imask = fh.variables['src_grid_imask'][:]
dst_imask = fh.variables['dst_grid_imask'][:]

src_mask = src_imask.__array__().reshape(src_ny, src_nx)
dst_mask = dst_imask.__array__().reshape(dst_ny, dst_nx)


import matplotlib.pyplot as pyplot

f, (ax1, ax2) = pyplot.subplots(2, sharex=True, sharey=True)
ax1.set_adjustable('box-forced')
ax2.set_adjustable('box-forced')
f, (ax3, ax4) = pyplot.subplots(2, sharex=True, sharey=True)
ax3.set_adjustable('box-forced')
ax4.set_adjustable('box-forced')
ax1.imshow(src_mask[::-1,:], cmap=pyplot.cm.bone)
ax2.imshow(src_weights[::-1,:], cmap=pyplot.cm.jet)

ax3.imshow(dst_mask[::-1,:], cmap=pyplot.cm.bone)
ax4.imshow(dst_weights[::-1,:], cmap=pyplot.cm.jet)


pyplot.ion()
pyplot.show()


raw_input()
