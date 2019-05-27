#!/usr/bin/env python

import numpy

from amuse.units import units

from amuse.datamodel import Grid

class pop_grid_reader(object):
  
  def __init__(self,filename,nx,ny):
    self.filename=filename

    import os.path
    import sys
    if not os.path.isfile(filename):
        sys.exit("Error: No such file " + filename)

    raw=numpy.fromfile(filename).newbyteorder('S')

    grid=raw.reshape((7,ny,nx))

#    lats=grid[0,:,:]/numpy.pi*180
#    lons=grid[1,:,:]/numpy.pi*180
    self.lat = lats = grid[0,:,:]
    self.lon = lons = grid[1,:,:]

    self.htn = htn = grid[2,:,:]
    self.hte = hte = grid[3,:,:]
    self.hus = hus = grid[4,:,:]
    self.huw = huw = grid[5,:,:]
    self.angle = angle = grid[6,:,:]

    #nodes is the U-grid
    self.nodes = nodes=Grid(ny*nx)
    nodes.lats = lats | units.rad
    nodes.lons = lons | units.rad

  def get_sets(self):
    return self.nodes



if __name__=="__main__":

    import sys
    total = len(sys.argv)
    if len(sys.argv) != 4:
        sys.exit("Usage: read_grid.py filename nx ny")
    filename = sys.argv[1]
    nx = int(sys.argv[2])
    ny = int(sys.argv[3])

    r=pop_grid_reader(filename, nx, ny)
    nodes=r.get_sets()
    
    print nodes

    print nodes.lats


