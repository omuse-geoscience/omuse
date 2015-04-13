import numpy

from amuse.units import units

from matplotlib import pyplot
 
from interface import QGmodel,QGmodelInterface

# run with dijkstra 2005 parameters
def sample_run(Nx=128,Ny=128,Reynolds_number=50.,wind_sigma=1.,dtplot=24.| units.hour,nstep=10):

  qg=QGmodel(redirection="none")
  qg.initialize_code()

  L=1.e6 | units.m

  qg.parameters.Lx=L
  qg.parameters.Ly=L
  qg.parameters.dx=qg.parameters.Lx/Nx
  qg.parameters.dy=qg.parameters.Ly/Ny
  qg.parameters.ocean_depth=600 | units.m
  qg.parameters.tau=0.15 | units.Pa
  qg.parameters.beta0=1.6e-11 | (units.m*units.s)**-1
  qg.parameters.dt=900 | units.s
  
  rho0=qg.parameters.rho
  beta0=qg.parameters.beta0
  H=qg.parameters.ocean_depth
  L=qg.parameters.Lx
  T=qg.parameters.tau

  U_dijkstra=1.6e-2 | units.m/units.s   # mentioned in paper
  U=(T/(beta0*rho0*H*L)) # actual derived from parameters
  print "actual, target U:", U.in_(units.m/units.s), U_dijkstra
    
  A_H=U*L/Reynolds_number
  
  qg.parameters.A_H=A_H
  
  qg.parameters.wind_sigma=wind_sigma
  
  timescale=1/(beta0*L)
  
  print "timescale:", timescale.in_(units.s)
  print
  print qg.parameters 
  raw_input()

  qg.commit_parameters()
  
#  qg.grid.psi=(qg.grid.x/L ) | units.m**2/units.s
  
  qg.initialize_grid()

  pyplot.ion()
  f=pyplot.figure()
  pyplot.show()


  for i in range(nstep+1):
    tend=i*dtplot

    qg.evolve_model(tend)
    psi=qg.grid.psi[:,:,0]

    print qg.model_time.in_(units.day),(psi.max()/L).in_(units.km/units.hour), \
      (psi.min()/L).in_(units.km/units.hour)
  
    f.clf()
    
    extent=[0,L.value_in(units.km),0,L.value_in(units.km)]
    pyplot.imshow(psi.transpose()/max(psi.max(),-psi.min()),vmin=-1,vmax=1,origin="lower",extent=extent)
    pyplot.xlabel("km")
    pyplot.title(str(qg.model_time.in_(units.day)))

    pyplot.draw()
  raw_input()
    



if __name__=="__main__":
  sample_run(Nx=256,Ny=256,Reynolds_number=60.,wind_sigma=0.)
  
