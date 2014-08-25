import numpy

from amuse.units import units
from amuse.datamodel import Grid
from interface import QGmodel,QGmodelWithRefinements,jans_wind_model

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

def test_evolve_w_plot(sysfac,tend=1. | units.hour,dt=3600. | units.s,dtplot=None):
  sys=sysfac()
  
  if dtplot is None:
    dtplot=dt

  pyplot.ion()
  f=pyplot.figure(figsize=(10,10))
  pyplot.show()
  
  i=0
  
  Lx=sys.parameters.Lx
  grid=Grid.create((400,400), (Lx,Lx))
  dx,dy=grid.cellsize()
  Lx=Lx.value_in(1000*units.km)  
  x=grid.x.flatten()
  y=grid.y.flatten()
  
  while sys.model_time<tend-dtplot/2:
    i=i+1
    sys.evolve_model(sys.model_time+dtplot,dt=dt)
    
    psi,dpsi=sys.get_psi_dpsidt(dx+0.*x,x,y)
    psi=psi.reshape(grid.shape)
#    psi=sys.grid[:,:,0].psi

    f.clf()
    f1=pyplot.subplot(111)
    f1.imshow(psi.transpose()/psi.max(),vmin=0,vmax=1,extent=[0,Lx,0,Lx],origin="lower")
    f1.set_xlabel("x (x1000 km)")
        
    pyplot.draw()
  raw_input()

def refinement_east(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2

  q1.add_refinement(q2,offset=[0,0]| units.m)

  return q1  

def refinement_north(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8
  q1.parameters.interface_wind=True
  jans_wind_model(q1.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Ly/=2
  q2.parameters.interface_wind=True
  jans_wind_model(q2.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q1.add_refinement(q2,offset=[0,0]| units.m)

  return q1

def refinement_south(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8
  q1.parameters.interface_wind=True
  jans_wind_model(q1.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Ly/=2
  q2.parameters.interface_wind=True

  q1.add_refinement(q2,offset=[0| units.m, q2.parameters.Ly])
  jans_wind_model(q2.wind_field,q1.parameters.Ly,q1.parameters.tau)

  return q1

def refinement_south_west(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8
  q1.parameters.interface_wind=True
  jans_wind_model(q1.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2
  q2.parameters.Ly/=2
  q2.parameters.interface_wind=True

  q1.add_refinement(q2,offset=[q2.parameters.Lx, q2.parameters.Ly])
  jans_wind_model(q2.wind_field,q1.parameters.Ly,q1.parameters.tau)

  return q1

def refinement_north_east_south(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8
  q1.parameters.interface_wind=True
  jans_wind_model(q1.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2
  q2.parameters.Ly/=2
  q2.parameters.interface_wind=True

  q1.add_refinement(q2,offset=[q2.parameters.Lx, q2.parameters.Ly/2])
  jans_wind_model(q2.wind_field,q1.parameters.Ly,q1.parameters.tau)

  return q1

def refinement_central(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8
  q1.parameters.interface_wind=True
  jans_wind_model(q1.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q2=QGmodel(redirection="none")  
  q2.parameters.dt=dt/2    
  q2.parameters.Lx/=2
  q2.parameters.Ly/=2
  q2.parameters.interface_wind=True

  q1.add_refinement(q2,offset=[q2.parameters.Lx/2, q2.parameters.Ly/2])
  jans_wind_model(q2.wind_field,q1.parameters.Ly,q1.parameters.tau)

  return q1

def nested_refinement(dt=3600. | units.s):
  q1=QGmodelWithRefinements(redirection="none")
  q1.parameters.dt=dt/2    
  q1.parameters.dx*=8
  q1.parameters.dy*=8
  q1.parameters.interface_wind=True
  jans_wind_model(q1.wind_field,q1.parameters.Ly,q1.parameters.tau)

  Lx=q1.parameters.Lx

  q2=QGmodelWithRefinements(redirection="none")  
  q2.parameters.dt=dt/4    
  q2.parameters.Lx/=2
  q2.parameters.dx*=4
  q2.parameters.dy*=4

  q3=QGmodel(redirection="none")  
  q3.parameters.dt=dt/4    
  q3.parameters.Lx/=2
  q3.parameters.Ly/=2
  q3.parameters.interface_wind=True
  jans_wind_model(q3.wind_field,q1.parameters.Ly,q1.parameters.tau)

  q1.add_refinement(q2,offset=[0,0]| units.m)
  q2.add_refinement(q3,offset=[0,0]| units.m)

  return q1
  
if __name__=="__main__":
  dt=1800 | units.s
  
  def sysfac():
    return refinement_east(dt)
  
  test_evolve_w_plot(sysfac,tend=100*dt,dt=dt,dtplot=4*dt)
