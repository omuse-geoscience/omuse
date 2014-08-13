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

if __name__=="__main__":
