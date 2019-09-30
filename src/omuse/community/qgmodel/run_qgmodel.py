import numpy

from os import path
 
#try:
#from amuse.community.qgmodel.interface import QGmodel,QGmodelInterface
#except ImportError as ex:
from interface import QGmodel,QGmodelInterface

from amuse.units import units

from matplotlib import pyplot


def low_level():

  q=QGmodelInterface(redirection="none")
  
  print 1
  q.initialize_code()
  print 2
  
  q.set_Lx(4.e6)
  q.set_Ly(4.e6)
  q.set_dx(1.e4)
  q.set_dy(1.e4)
  q.set_dt(1800)
  
  q.commit_parameters()
  print 3
  
  q.initialize_grid()
  print 4
  
  q.evolve_model(86400.)
  print 5
  print q.get_time()
  
  x,y=numpy.mgrid[0:400,0:400]
  
  x=x.flatten()+1
  y=y.flatten()+1
  
  psi,err= q.get_psi1_state(x,y,1)
  
  psi=psi.reshape((400,400))
  
  print psi.shape
  pyplot.imshow(psi)
  
  pyplot.show()

def high_level():

  qg=QGmodel(redirection="none")
  qg.initialize_code()

  Nx=128
  Ny=128

  qg.parameters.Lx=1.e6 | units.m
  qg.parameters.Ly=1.e6 | units.m
  qg.parameters.dx=qg.parameters.Lx/Nx
  qg.parameters.dy=qg.parameters.Ly/Ny
  qg.parameters.ocean_depth=600 | units.m
  qg.parameters.tau=0.15 | units.Pa
  qg.parameters.beta0=1.6e-11 | (units.m*units.s)**-1

  qg.parameters.dt=1800/10 | units.s
  
  rho0=qg.parameters.rho
  beta0=qg.parameters.beta0
  H=qg.parameters.ocean_depth
  L=qg.parameters.Lx
  T=qg.parameters.tau

  U_dijkstra=1.6e-2 | units.m/units.s
  U=(T/(beta0*rho0*H*L)) 
  print "actual, target U:", U.in_(units.m/units.s), U_dijkstra
  
  
  Reynolds_number=1.
  
  A_H=U*L/Reynolds_number
  
  qg.parameters.A_H=A_H
  
  timescale=1/(beta0*L)
  
  print "timescale:", timescale.in_(units.s)
  print qg.parameters 

  qg.commit_parameters()
  
#  qg.grid.psi=(qg.grid.x/L ) | units.m**2/units.s
  
  qg.initialize_grid()

  pyplot.ion()
  f=pyplot.figure()
  pyplot.show()

  for i in range(101):
    tend=i*12.| units.hour

    qg.evolve_model(tend)
    psi=qg.grid.psi[:,:,0]

    print qg.model_time.in_(units.day),psi.max(),psi.min()
  
    pyplot.imshow(psi.transpose()/psi.max(),vmin=0,vmax=1,origin="lower")

    pyplot.draw()
  raw_input()
    



if __name__=="__main__":
  high_level()
  
