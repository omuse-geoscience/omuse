import numpy

from os import path
 
#try:
#from amuse.community.qgmodel.interface import QGmodel,QGmodelInterface
#except ImportError as ex:
from interface import QGmodel,QGmodelInterface

from amuse.units import units

from matplotlib import pyplot

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

  qg.parameters.dt=1800 | units.s
  
  qg.parameters.free_slip=False
  
  rho0=qg.parameters.rho
  beta0=qg.parameters.beta0
  H=qg.parameters.ocean_depth
  L=qg.parameters.Lx
  T=qg.parameters.tau

  U_dijkstra=1.6e-2 | units.m/units.s
  U=(T/(beta0*rho0*H*L)) 
  print "actual, target U:", U.in_(units.m/units.s), U_dijkstra
  
  
  Reynolds_number=10.
  
  A_H=U*L/Reynolds_number
  
  qg.parameters.A_H=A_H
  
  timescale=1/(beta0*L)
  
  print "timescale:", timescale.in_(units.s)
  print qg.parameters 

  qg.commit_parameters()
  
#  qg.grid.psi=(qg.grid.x/L ) | units.m**2/units.s
  
  qg.initialize_grid()

  pyplot.ion()
  f=pyplot.figure(figsize=(12,6))
  pyplot.show()

  dtplot=12.| units.hour

  for i in range(101):
    tend=i*dtplot

    qg.evolve_model(tend)
    psi=qg.grid.psi[:,:,0]

    print qg.model_time.in_(units.day),psi.max(),psi.min()
  
    f.clf()
    f1=pyplot.subplot(121)
    f1.imshow(psi.transpose()/psi.max(),vmin=0,vmax=1,origin="lower")

    f2=pyplot.subplot(122)
    x=qg.grid.x[:,64,0].value_in(units.km)
    psi=qg.grid.psi[:,64,0]
    psi=psi/psi.max()
    
    f2.plot(x,psi)

    pyplot.draw()
    c=raw_input()
    if c=="q": break



if __name__=="__main__":
  high_level()
  
