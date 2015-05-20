from omuse.units import units

def EOS_linear(T, S, rho0=1028. | units.kg/units.m**3, T0=10. | units.C, S0=35. | units.psu, alpha=0.00017 | units.C**-1, beta=0.00076 | units.psu**-1):
    """
    the equation of state given in Cushman-Roisin, B., Introduction to Geophysical Fluid Dynamics, Prentice-Hall, 1994, 320 pp. (1994)
    (as reported by ADCIRC)
    """
    return rho0 * (1-alpha*(T-T0)+beta*(S-S0))
  
  
if __name__=="__main__":
  print EOS_linear(10 | units.C, 35. | units.psu)
  print EOS_linear.__doc__
