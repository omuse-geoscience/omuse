#!/usr/bin/python

import numpy

import matplotlib.pyplot as pyplot
from netCDF4 import Dataset

class weights_file_reader():

    def __init__(self,filename):

        fh = Dataset(filename, mode='r')

        self.src_dims = src_dims = fh.variables['src_grid_dims'][::-1]
        self.dst_dims = dst_dims = fh.variables['dst_grid_dims'][::-1]

        num_links = len(fh.dimensions['num_links'])
        remap_matrix = fh.variables['remap_matrix'][:]

        #print min(fh.variables['src_address']), max(fh.variables['src_address'])
        #print min(fh.variables['dst_address']), max(fh.variables['dst_address'])
        #print 'src_ny=', src_ny, ' src_nx=', src_nx

        src_weights = numpy.zeros(len(fh.dimensions['src_grid_size']))
        dst_weights = numpy.zeros(len(fh.dimensions['dst_grid_size']))
        src_address = fh.variables['src_address'][:]
        dst_address = fh.variables['dst_address'][:]

        num_wts = len(fh.dimensions['num_wgts'])
        remap_links = remap_matrix.flatten()
        for i in range(num_links):
            src_weights[src_address[i] -1] += remap_links[i*num_wts] #correct for Fortran indexing starting at 1
            dst_weights[dst_address[i] -1] += remap_links[i*num_wts] #correct for Fortran indexing starting at 1

#        for (i,j), weight in numpy.ndenumerate(remap_matrix):
#            src_i = src_address[i] -1 #correct for Fortran indexing starting at 1
#            src_weights[src_i] += weight
#            dst_i = dst_address[i] -1 #correct for Fortran indexing starting at 1
#            dst_weights[dst_i] += weight

        self.src_weights = src_weights.reshape(src_dims)
        self.dst_weights = dst_weights.reshape(dst_dims)

        src_imask = fh.variables['src_grid_imask'][:]
        dst_imask = fh.variables['dst_grid_imask'][:]

        self.src_mask = src_mask = src_imask.__array__().reshape(src_dims)
        self.dst_mask = dst_mask = dst_imask.__array__().reshape(dst_dims)


def create_window():
    f, (ax1, ax2) = pyplot.subplots(2, sharex=True, sharey=True)
    f.tight_layout()
    ax1.set_adjustable('box-forced')
    ax2.set_adjustable('box-forced')
    return ax1, ax2


def view_weights_scrip(gridfile, weights):
    fh = Dataset(gridfile, mode='r')

    dims = fh.variables['grid_dims'][::-1]
    imask = fh.variables['grid_imask'][:]
    mask = imask.__array__().reshape(dims)

    ax1, ax2 = create_window()
    ax1.imshow(mask[::-1,:], cmap=pyplot.cm.bone)
    ax2.imshow(weights.reshape(dims)[::-1,:], cmap=pyplot.cm.jet)


def view_weights_adcirc(gridfile, weights):

    from omuse.community.cdo.view_adcirc_grid import adcirc_grid_viewer
    v = adcirc_grid_viewer(filename=gridfile, coordinates='spherical')

    ax3, ax4 = create_window()
    ax3.triplot(v.x, v.y, v.triangles)
    ax4.tripcolor(v.x, v.y, v.triangles, weights)



if __name__ == "__main__":

    import sys
    import os.path

    total = len(sys.argv)
    if not len(sys.argv) in [2, 6]:
        sys.exit("Usage: view_weights.py weightsfile [src_grid_file [scrip|adcirc] dst_grid_file [scrip|adcirc]]")

    use_mask_from_weightsfile = False
    if len(sys.argv) == 2:
        use_mask_from_weightsfile = True

    weightsfile = sys.argv[1]
    files = []
    files.append(weightsfile)

    if not use_mask_from_weightsfile:
        src_grid_file = sys.argv[2]
        files.append(src_grid_file)
        dst_grid_file = sys.argv[4]
        files.append(dst_grid_file)

    for filename in files:
        if not os.path.isfile(filename):
            sys.exit("Error: No such file " + filename)

    pyplot.ion()

    r=weights_file_reader(weightsfile)

    if use_mask_from_weightsfile:
        ax1, ax2 = create_window()
        ax1.imshow(r.src_mask[::-1,:], cmap=pyplot.cm.bone)
        ax2.imshow(r.src_weights.reshape(r.src_dims)[::-1,:], cmap=pyplot.cm.jet)

        ax3, ax4 = create_window()
        ax3.imshow(r.dst_mask[::-1,:], cmap=pyplot.cm.bone)
        ax4.imshow(r.dst_weights.reshape(r.dst_dims)[::-1,:], cmap=pyplot.cm.jet)

    
    else:
        if sys.argv[3] == "scrip":
            view_weights_scrip(src_grid_file, r.src_weights)
        elif sys.argv[3] == "adcirc":
            view_weights_adcirc(src_grid_file, r.src_weights)

        if sys.argv[5] == "scrip":
            view_weights_scrip(dst_grid_file, r.dst_weights)
        elif sys.argv[5] == "adcirc":
            view_weights_adcirc(dst_grid_file, r.dst_weights)
    

    pyplot.show()

    raw_input()

