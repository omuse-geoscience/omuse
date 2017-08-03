from amuse.units.constants import *
from amuse.units.trigo import sin,cos
from amuse.units import units

Rearth=6371.0088 | units.km

def coriolis_frequency(lat,omega=(1.|units.rev)/(sidereal_day)):
  return 2*omega*sin(lat)
def coriolis_beta(lat,omega=(1.|units.rev)/(sidereal_day)):
  return 2*omega*cos(lat)/Rearth
