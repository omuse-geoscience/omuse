import numpy

from amuse.units import units
from amuse.datamodel import Grid
from interface import QGmodel

from matplotlib import pyplot

#import logging
#logging.basicConfig(level=logging.DEBUG)
#logging.getLogger("code").setLevel(logging.DEBUG)


def reference(tend=1. | units.hour,dt=3600. | units.s):

  q=QGmodel(redirection="none")
  
  q.parameters.dt=dt
    
  q.evolve_model(tend)
  
  psi=q.grid[:,:,0].psi.number

  return psi

def interface_bc(tend=1. | units.hour,dt=3600. | units.s,correct=True):

  q=QGmodel(redirection="none")
  
  q.parameters.dt=dt

  q.parameters.xbound1="interface"
    
  while q.model_time<tend:
    q.evolve_model(q.model_time+dt)
    if correct:
      xb1=q.boundaries(1).copy()
      xb1.psi[0,1:-1,0]=-q.grid.psi[1,:,0]
      channel=xb1.new_channel_to(q.boundaries(1))
      channel.copy()
      
  psi=q.grid[:,:,0].psi.number

  return psi

def test1():
  
  tend=192. | units.hour
  psi1=reference(tend)
  psi2=interface_bc(tend)
    
  d=abs(psi2-psi1)
  print d.max(),d.mean(),abs(psi1).max()
  
  pyplot.ion()
  f=pyplot.figure(figsize=(12,4))
  pyplot.show()

  f.clf()
  f1=pyplot.subplot(131)
  f1.imshow(psi1.transpose()/psi1.max(),vmin=0,vmax=1,origin="lower")

  f2=pyplot.subplot(132)
  f2.imshow(psi2.transpose()/psi1.max(),vmin=0,vmax=1,origin="lower")

  f3=pyplot.subplot(133)
  f3.imshow(abs(psi2-psi1).transpose(),vmin=0,origin="lower")

  pyplot.draw()
  raw_input()

def semi_domain_test(tend=1. | units.hour,dt=3600. | units.s):
  
  q1=QGmodel(redirection="none")  
  q1.parameters.dt=dt/2    
  Lx=q1.parameters.Lx.value_in(1000*units.km)
  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2
  q2.parameters.boundary_east="interface"
  
  Nx,Ny,Nm=q1.grid.shape
  
  pyplot.ion()
  f=pyplot.figure(figsize=(12,5))
  pyplot.show()
  
  i=0
  while q1.model_time<tend:
    i=i+1
    tnow=q1.model_time
    q1.evolve_model(tnow+dt/2)
    psi=q1.grid[(Nx-1)/2:(Nx-1)/2+2,:,0].psi
    dpsi_dt=q1.grid[(Nx-1)/2:(Nx-1)/2+2,:,0].dpsi_dt
    west=q2.boundaries("west").copy()
    west[:,1:-1,0].psi=psi
    west[:,1:-1,0].dpsi_dt=dpsi_dt
    channel=west.new_channel_to(q2.boundaries("west"))
    channel.copy()
    q1.evolve_model(tnow+dt)
    q2.evolve_model(tnow+dt)

#  print q1.grid[(Nx+1)/2,1,0].dpsi_dt
#  print q2.boundaries("west")[0:2,0,0].x
    if i%5==0:  
      psi1_complete=q1.grid.psi.number[:,:,0]
      psi1=q1.grid[:(Nx+1)/2,:,0].psi.number
      psi2=q2.grid[:,:,0].psi.number
    
      d=abs(psi2-psi1)
      print d.max(),d.mean(),abs(psi1).max()
        
      f.clf()
      f1=pyplot.subplot(131)
      f1.imshow(psi1_complete.transpose()/psi1.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
      f1.set_xlabel("x (x1000 km)")
    
      f2=pyplot.subplot(132)
      f2.imshow(psi2.transpose()/psi1.max(),vmin=0,vmax=1,extent=[0,Lx/2,0,Lx],origin="lower")
      f2.set_xlabel("x (x1000 km)")
    
      f3=pyplot.subplot(133)
      f3.imshow(d.transpose()/psi1.max(),vmin=0,vmax=0.001,extent=[0,Lx/2,0,Lx],origin="lower")
      f3.set_xlabel("x (x1000 km)")
    
      pyplot.draw()
      pyplot.savefig("snapshots/half_domain_test-%6.6i.png"%i)

  raw_input()

def semi_domain_test_interpolation(tend=1. | units.hour,dt=3600. | units.s):
  
  q1=QGmodel(redirection="none")  
  q1.parameters.dt=dt/2    
  Lx=q1.parameters.Lx.value_in(1000*units.km)
  dx=q1.parameters.dx
  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2
  q2.parameters.boundary_east="interface"
  
  Nx,Ny,Nm=q1.grid.shape
  
  pyplot.ion()
  f=pyplot.figure(figsize=(12,5))
  pyplot.show()
  
  i=0
  while q1.model_time<tend:
    i=i+1
    tnow=q1.model_time
    q1.evolve_model(tnow+dt/2)
    west=q2.boundaries("west").copy()    
    x=west[:,1:-1,0].x.flatten()
    y=west[:,1:-1,0].y.flatten()
    psi,dpsi_dt=q1.get_psi_state_at_point(0.*x+dx,x,y)
    west[:,1:-1,0].psi=psi.reshape(west[:,1:-1,0].shape)
    west[:,1:-1,0].dpsi_dt=dpsi_dt.reshape(west[:,1:-1,0].shape)
    
    channel=west.new_channel_to(q2.boundaries("west"))
    channel.copy()
    q1.evolve_model(tnow+dt)
    q2.evolve_model(tnow+dt)

#  print q1.grid[(Nx+1)/2,1,0].dpsi_dt
#  print q2.boundaries("west")[0:2,0,0].x
    if i%5==0:  
      psi1_complete=q1.grid.psi.number[:,:,0]
      psi1=q1.grid[:(Nx+1)/2,:,0].psi.number
      psi2=q2.grid[:,:,0].psi.number
    
      d=abs(psi2-psi1)
      print d.max(),d.mean(),abs(psi1).max()
        
      f.clf()
      f1=pyplot.subplot(131)
      f1.imshow(psi1_complete.transpose()/psi1.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
      f1.set_xlabel("x (x1000 km)")
    
      f2=pyplot.subplot(132)
      f2.imshow(psi2.transpose()/psi1.max(),vmin=0,vmax=1,extent=[0,Lx/2,0,Lx],origin="lower")
      f2.set_xlabel("x (x1000 km)")
    
      f3=pyplot.subplot(133)
      f3.imshow(d.transpose()/psi1.max(),vmin=0,vmax=0.001,extent=[0,Lx/2,0,Lx],origin="lower")
      f3.set_xlabel("x (x1000 km)")
    
      pyplot.draw()
      pyplot.savefig("snapshots/half_domain_test-%6.6i.png"%i)

  raw_input()

def semi_domain_test_multires(tend=1. | units.hour,dt=3600. | units.s):
  
  q1=QGmodel(redirection="none")  
  q1.parameters.dt=dt/2    
  Lx=q1.parameters.Lx.value_in(1000*units.km)
  dx=q1.parameters.dx*8
  q1.parameters.dx=dx
  q1.parameters.dy=dx
  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2
  q2.parameters.boundary_east="interface"
  
  Nx,Ny,Nm=q1.grid.shape
  
  pyplot.ion()
  f=pyplot.figure(figsize=(12,5))
  pyplot.show()
  
  i=0
  while q1.model_time<tend:
    i=i+1
    tnow=q1.model_time
    q1.evolve_model(tnow+dt/2)
    west=q2.boundaries("west").copy()    
    x=west[:,1:-1,0].x.flatten()
    y=west[:,1:-1,0].y.flatten()
    psi,dpsi_dt=q1.get_psi_state_at_point(0.*x+dx,x,y)
    west[:,1:-1,0].psi=psi.reshape(west[:,1:-1,0].shape)
    west[:,1:-1,0].dpsi_dt=dpsi_dt.reshape(west[:,1:-1,0].shape)
    
    channel=west.new_channel_to(q2.boundaries("west"))
    channel.copy()
    q1.evolve_model(tnow+dt)
    q2.evolve_model(tnow+dt)

#  print q1.grid[(Nx+1)/2,1,0].dpsi_dt
#  print q2.boundaries("west")[0:2,0,0].x
    if i%5==0:  
      psi1_complete=q1.grid.psi.number[:,:,0]
      psi1=q1.grid[:(Nx+1)/2,:,0].psi.number
      psi2=q2.grid[:,:,0].psi.number
    
#      d=abs(psi2-psi1)
#      print d.max(),d.mean(),abs(psi1).max()
        
      f.clf()
      f1=pyplot.subplot(131)
      f1.imshow(psi1_complete.transpose()/psi1.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
      f1.set_xlabel("x (x1000 km)")
    
      f2=pyplot.subplot(132)
      f2.imshow(psi2.transpose()/psi1.max(),vmin=0,vmax=1,extent=[0,Lx/2,0,Lx],origin="lower")
      f2.set_xlabel("x (x1000 km)")
    
#      f3=pyplot.subplot(133)
#      f3.imshow(d.transpose()/psi1.max(),vmin=0,vmax=0.001,extent=[0,Lx/2,0,Lx],origin="lower")
#      f3.set_xlabel("x (x1000 km)")
    
      pyplot.draw()
      pyplot.savefig("snapshots/half_domain_test-%6.6i.png"%i)

  raw_input()


class QGmodelWithRefinements(QGmodel):
  def __init__(self,*args,**kwargs):
    self.refinements=[]
    QGmodel.__init__(self,*args,**kwargs)
    
  def interpolate_grid(self,grid):
    copy=grid.empty_copy()    
    x=grid.x.flatten()
    y=grid.y.flatten()
    dx=max(self.grid.cellsize()[0],grid.cellsize()[0])
    psi,dpsi_dt=self.get_psi_state_at_point(dx+0.*x,x,y)
    copy.psi=psi.reshape(copy.shape)
    copy.dpsi_dt=dpsi_dt.reshape(copy.shape)    
    channel=copy.new_channel_to(grid)
    channel.copy_attributes(["psi","dpsi_dt"])

  def add_refinement(self,sys,offset=[0.,0.] | units.m):
    self.refinements.append(sys)
    sys.parameters.xoffset=offset[0]
    sys.parameters.yoffset=offset[1]
    
    minpos=sys.grid.get_minimum_position()
    maxpos=sys.grid.get_maximum_position()
    
    parent_minpos=self.grid.get_minimum_position()
    parent_maxpos=self.grid.get_maximum_position()
    dx,dy=self.grid.cellsize()
    
    if minpos[0]<parent_minpos[0]+dx/2:
      print "linking west boundary"
      sys.parameters.boundary_west=self.parameters.boundary_west
      sys.west_boundary_updater=self.west_boundary_updater
    else:
      print "west boundary interpolated"
      sys.parameters.boundary_west="interface"
      sys.west_boundary_updater=self.interpolate_grid
      
    if maxpos[0]>parent_maxpos[0]-dx/2:
      print "linking east boundary"
      sys.parameters.boundary_east=self.parameters.boundary_east
      sys.east_boundary_updater=self.east_boundary_updater
    else:
      print "east boundary interpolated"
      sys.parameters.boundary_east="interface"
      sys.east_boundary_updater=self.interpolate_grid

    if minpos[1]<parent_minpos[1]+dy/2:
      print "linking south boundary"
      sys.parameters.boundary_south=self.parameters.boundary_south
      sys.south_boundary_updater=self.south_boundary_updater
    else:
      print "south boundary interpolated"
      sys.parameters.boundary_south="interface"
      sys.south_boundary_updater=self.interpolate_grid
      
    if maxpos[1]>parent_maxpos[1]-dy/2:
      print "linking north boundary"
      sys.parameters.boundary_north=self.parameters.boundary_north
      sys.north_boundary_updater=self.north_boundary_updater
    else:
      print "north boundary interpolated"
      sys.parameters.boundary_north="interface"
      sys.north_boundary_updater=self.interpolate_grid
    if not hasattr(sys,"get_psi_dpsidt"):
      sys.get_psi_dpsidt=sys.get_psi_state_at_point

  def evolve_model(self,tend,dt=None):
    if dt is None:
      dt=2*self.parameters.dt
    tnow=self.model_time
    while tnow<tend-dt/2:
      self.overridden().evolve_model(tnow+dt/2)
      for sys in self.refinements:
        print "update boundaries...",
        sys.update_boundaries()
        print "done"
        sys.evolve_model(tnow+dt)
      self.overridden().evolve_model(tnow+dt)
#      self.update_refined_regions()
      tnow=self.model_time

  def get_psi_dpsidt(self,dx,x,y,k=None):
    minx=x-dx/2
    maxx=x+dx/2
    miny=y-dx/2
    maxy=y+dx/2
    done=numpy.zeros(minx.shape,dtype=bool)
    psi=numpy.zeros(minx.shape) | units.m**2/units.s
    dpsi=numpy.zeros(minx.shape) | units.m**2/units.s**2
    for sys in self.refinements+[self]:
      gridminx,gridminy=sys.grid.get_minimum_position()
      gridmaxx,gridmaxy=sys.grid.get_maximum_position()
      select=- ( (gridminx>minx)+(gridmaxx<maxx)+(gridminy>miny)+(gridmaxy<maxy)+done )
      if select.sum() and sys is not self: 
        p,d=sys.get_psi_dpsidt(dx[select],x[select],y[select])
        psi[select]=p
        dpsi[select]=d
      done=done+select
    if select.sum():
      p,d=self.get_psi_state_at_point(dx[select],x[select],y[select])
      psi[select]=p
      dpsi[select]=d
    print "points done:", done.sum()
    print "points skipped:", (1-done).sum()
    return psi,dpsi

def test_refinement_east(tend=1. | units.hour,dt=3600. | units.s,dtplot=None):
  if dtplot is None:
    dtplot=dt
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2

  q1.add_refinement(q2,offset=[0,0]| units.m)
  
  pyplot.ion()
  f=pyplot.figure(figsize=(10,10))
  pyplot.show()
  
  i=0
  
  Lx=q1.parameters.Lx
  grid=Grid.create((400,400), (Lx,Lx))
  dx,dy=grid.cellsize()
  Lx=q1.parameters.Lx.value_in(1000*units.km)  
  Nx,Ny,Nm=q1.grid.shape
  
  while q1.model_time<tend-dtplot/2:
    i=i+1
    q1.evolve_model(q1.model_time+dtplot,dt=dt)
    
    x=grid.x.flatten()
    y=grid.y.flatten()
    psi,dpsi=q1.get_psi_dpsidt(dx+0.*x,x,y)
    psi=psi.reshape(grid.shape)

    f.clf()
    f1=pyplot.subplot(111)
    f1.imshow(psi.transpose()/psi.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
    f1.set_xlabel("x (x1000 km)")
        
    pyplot.draw()
  raw_input()

def test_refinement_north(tend=1. | units.hour,dt=3600. | units.s,dtplot=None):
  if dtplot is None:
    dtplot=dt
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Ly/=2

  q1.add_refinement(q2,offset=[0,0]| units.m)
  
  pyplot.ion()
  f=pyplot.figure(figsize=(10,10))
  pyplot.show()
  
  i=0
  
  Lx=q1.parameters.Lx
  grid=Grid.create((400,400), (Lx,Lx))
  dx,dy=grid.cellsize()
  Lx=q1.parameters.Lx.value_in(1000*units.km)  
  Nx,Ny,Nm=q1.grid.shape
  
  while q1.model_time<tend-dtplot/2:
    i=i+1
    q1.evolve_model(q1.model_time+dtplot,dt=dt)
    
    x=grid.x.flatten()
    y=grid.y.flatten()
    psi,dpsi=q1.get_psi_dpsidt(dx+0.*x,x,y)
    psi=psi.reshape(grid.shape)

    f.clf()
    f1=pyplot.subplot(111)
    f1.imshow(psi.transpose()/psi.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
    f1.set_xlabel("x (x1000 km)")
        
    pyplot.draw()
  raw_input()


def test_2_refinements(tend=1. | units.hour,dt=3600. | units.s,dtplot=None):
  if dtplot is None:
    dtplot=dt
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=16
  q1.parameters.dy*=16

  Lx=q1.parameters.Lx

  q2=QGmodelWithRefinements(redirection="none")  
  q2.parameters.dt=dt/4    
  q2.parameters.Lx/=2
  q2.parameters.dx*=8
  q2.parameters.dy*=8

  q3=QGmodel(redirection="none")  
  q3.parameters.dt=dt/4    
  q3.parameters.Lx/=2
  q3.parameters.Ly/=2

  q1.add_refinement(q2,offset=[0,0]| units.m)
  q2.add_refinement(q3,offset=[0,0]| units.m)

  
  pyplot.ion()
  f=pyplot.figure(figsize=(12,5))
  pyplot.show()

#  Lx=q1.parameters.Lx.value_in(1000*units.km)  
#  dx=q1.parameters.dx

  raise
  i=0
  
  grid=Grid.create((400,400), (Lx,Lx))
  dx,dy=grid.cellsize()
  Lx=q1.parameters.Lx.value_in(1000*units.km)  
  Nx,Ny,Nm=q1.grid.shape
  
  while q1.model_time<tend-dtplot/2:
    i=i+1
    q1.evolve_model(q1.model_time+dtplot,dt=dt)
    
    x=grid.x.flatten()
    y=grid.y.flatten()
    psi,dpsi=q1.get_psi_dpsidt(dx+0.*x,x,y)
    psi=psi.reshape(grid.shape)

    f.clf()
    f1=pyplot.subplot(121)
    f1.imshow(psi.transpose()/psi.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
    f1.set_xlabel("x (x1000 km)")
        
    pyplot.draw()
  raw_input()


if __name__=="__main__":
  dt=1800 | units.s
#  test_refinement_east(tend=100*dt,dt=dt,dtplot=4*dt)
  test_refinement_north(tend=100*dt,dt=dt,dtplot=4*dt)
#  test_2_refinements(tend=100*dt,dt=dt,dtplot=4*dt)