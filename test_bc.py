import numpy

from amuse.units import units
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
  q2.parameters.xbound2="interface"
  
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


if __name__=="__main__":
  dt=1800 | units.s
  semi_domain_test(tend=500*dt,dt=dt)
