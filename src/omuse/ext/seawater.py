from amuse.units.trigo import *
from omuse.units import units
from . import _seawater

def dens(S,T,P=0. | units.dbar):
    """  Density of sea water """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.dens(_S, _T, _P) | units.kg/units.m**3

def svan(S,T,P=0. | units.dbar):
    """ Specific volume anomaly """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.svan(_S,_T,_P) | units.m**3/units.kg
delta=svan

def sigma(S,T,P=0. | units.dbar):
    """ Density anomaly """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.sigma(_S,_T,_P) | units.kg/units.m**3

def drhodt(S,T,P=0. | units.dbar):
    """Temperature derivative of density"""
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.drhodt(_S,_T,_P) | units.kg/(units.K*units.m**3)

def alpha(S,T,P=0. | units.dbar):
    """ Thermal expansion coefficient """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.alpha(_S,_T,_P) | units.K**-1
    
def drhods(S,T,P=0. | units.dbar):
    """ Salinity derivative of density """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.drhods(_S,_T,_P) | units.kg/units.m**3

def beta(S,T,P=0. | units.dbar):
    """ Salinity expansion coefficient """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.beta(_S,_T,_P)

def salt(R_cond, T, P=0. | units.dbar):
    """ Salinity """
    _R=R_cond
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.salt( _R,_T,_P)

def R_cond(S,T,P=0. | units.dbar):
    """ Conductivity ratio """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.cond(_S,_T,_P)
cond=R_cond

def heatcap(S,T,P=0. | units.dbar):
    """ Heat capacity """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.heatcap( _S,_T,_P) | units.J/(units.kg*units.K)
    
def adtgrad(S,T,P=0. | units.dbar):
    """ Adiabatic lapse rate """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)   
    return _seawater.adtgrad( _S,_T,_P) | units.K/units.dbar
    
def temppot(S,T,P=0. | units.dbar ,Pref=0. | units.dbar):
    """ Potential temperature """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)
    _Pref=Pref.value_in(units.dbar)
    return _seawater.temppot( _S,_T,_P,_Pref) | units.Celsius

def temppot0(S,T,P=0. | units.dbar):
    """ Potential temperature """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)
    return _seawater.temppot0( _S,_T,_P) | units.Celsius

def freezept(S,P=0. | units.dbar):
    """ Freezing point """
    _S=S.value_in(units.psu)
    _P=P.value_in(units.dbar)
    return _seawater.freezept(_S,_P) | units.Celsius

def soundvel(S,T,P=0. | units.dbar):
    """ Sound velocity  """
    _S=S.value_in(units.psu)
    _T=T.value_in(units.Celsius)
    _P=P.value_in(units.dbar)
    return _seawater.soundvel(_S,_T,_P) | units.m/units.s
    
def depth(P,lat):
    """ depth from pressure and lattitude """
    _P=P.value_in(units.dbar)
    _lat=to_deg(lat)
    return _seawater.depth(_P,_lat) | units.m
   
def EOS_linear(S, T, rho0=1028. | units.kg/units.m**3, T0=10. | units.Celsius, S0=35. | units.psu, alpha=0.00017 | units.Celsius**-1, beta=0.00076 | units.psu**-1):
    """
    the equation of state given in Cushman-Roisin, B., Introduction to Geophysical Fluid Dynamics, Prentice-Hall, 1994, 320 pp. (1994)
    (as reported by ADCIRC)
    """
    return rho0 * (1-alpha*(T-T0)+beta*(S-S0))

def EOS_UNESCO(S,T, P=0 | units.Pa):
    """
    Unesco seawater density as imported from Bjorn Adlandsvik's seawater package
    """
    return dens(S,T,P)

if __name__=="__main__":
  S=35. | units.psu
  T=10 | units.Celsius
  P=0. | units.dbar
  print(EOS_linear( S,T))
  print(EOS_linear.__doc__)
  print(EOS_UNESCO( S,T))
  print(EOS_UNESCO.__doc__)

  for x in [dens,svan,sigma,drhodt,alpha,drhods,beta, R_cond,cond, heatcap,
       adtgrad,temppot,temppot0,soundvel]:
      print(x(S,T,P))
      print(x.__doc__)
      
  R=cond(S,T,P)    
  print(salt(R,T,P))
  print(salt.__doc__)
  print(freezept(S,P))
  print(freezept.__doc__)
  print(depth(P,0.))
  print(depth.__doc__)
  
  


