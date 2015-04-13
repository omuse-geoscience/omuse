import numpy

from os import path
 
from interface import QGmodel,single_gyre_wind_model

from amuse.units import units

from matplotlib import pyplot

from amuse.io import write_set_to_file

def viebahn2014(N=50,reynolds_number=1,dm=0.04):

  beta0=1.8616e-11 |(units.m * units.s)**-1
  L=1.e6 | units.m
  H=4000.| units.m
  rho=1000. | units.kg/units.m**3
  dx=L/N

  A_H=beta0*(dm*L)**3
  tau=reynolds_number*A_H*rho*beta0*H
  U=tau/beta0/rho/H/L
  delta_m=(A_H/beta0)**(1./3)/L
  delta_i=(U/beta0)**0.5/L
  timescale=1./(beta0*L)

  print "Viebahn 2014 setup"
  print "N=%i, Reynolds_number=%f"%(N,reynolds_number)
  print "dm (derived):", (A_H/beta0)**(1./3)/L
  print "tau:", tau.value_in(units.Pa)
  print "A:", A_H
  print "timescale:", timescale.in_(units.s)
  print "delta_m:", delta_m
  print "delta_i:", delta_i

  qg=QGmodel(redirection="none")

  qg.parameters.Lx=L
  qg.parameters.Ly=L
  qg.parameters.dx=dx
  qg.parameters.dy=dx
  qg.parameters.dt=900 | units.s
  qg.parameters.A_H=A_H
  qg.parameters.interface_wind=True
  qg.parameters.rho=rho
  qg.parameters.beta0=beta0
  qg.parameters.ocean_depth=H  
  qg.parameters.tau=tau
    
  def wind_function(x,y):
    return single_gyre_wind_model(x,y,L,tau)
  
  qg.set_wind(wind_function)

  return qg


def evolve_to_eq(qg,f=0.01,label=""):

  dtplot=10.| units.day

  psi=qg.grid[:,:,0].psi
  for i in range(101):
    tend=i*dtplot

    qg.evolve_model(tend)
    prev=psi
    psi=qg.grid[:,:,0].psi

    d=abs(psi-prev).sum()/psi.sum()
    print i,d
    if d<0.01: break
  
  write_set_to_file(qg.grid,"viebahn2014_grid"+label,"amuse")

if __name__=="__main__":
  sys=viebahn2014(400,20.)  
  evolve_to_eq(sys,label="_reference")
