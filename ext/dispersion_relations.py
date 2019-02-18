import numpy

from omuse.units import units

# todo:
# k_from_omega with root finding
# finite rho_air effects: interface waves
# capillary only

class DispersionRelation(object):
    def omega2(self, k, h=None):
        raise Exception("implement")

    def domega2_dk(self, k, h=None):
        raise Exception("implement")

    def domega_dk(self, k, h=None):
        return 1./(2*self.omega(k,h))*self.domega2_dk(k,h)

    def omega(self, k,h=None):
        return self.omega2(k,h)**0.5

    def group_velocity(self,k,h=None):
        return self.domega_dk(k,h)

    def phase_velocity(self,k,h=None):
        return self.omega(k,h)/k

    def k_from_omega(self,omega,h=None):
        pass

    def omega_from_phase_velocity(self,v,h=None):
        pass

    def k_from_lambda(self, l):
        return 2*numpy.pi/l

class SurfaceGravityWaves(DispersionRelation):
    def __init__(self, g=9.81 | units.m/units.s**2,depth=None):
        self.g=g
        self.depth=depth
    def omega2(self,k,h=None):
        h=h or self.depth 
        return self.g*k*numpy.tanh(k*h)
    def domega2_dk(self,k,h=None):
        h=h or self.depth
        th=numpy.tanh(k*h)
        return self.g*(th+k*h*(1-th**2))
  
class DeepWaterWaves(DispersionRelation):
    def __init__(self, g=9.81 | units.m/units.s**2):
        self.g=g
    def omega2(self,k,h=None):
        return self.g*k
    def domega2_dk(self,k,h=None):
        return self.g
    def k_from_omega(self,omega,h=None):
        return omega**2/self.g
    def omega_from_phase_velocity(self,v,h=None):
        return self.g/v

class ShallowWaterWaves(DispersionRelation):
    def __init__(self, g=9.81 | units.m/units.s**2,depth=None):
        self.g=g
        self.depth=depth
    def omega2(self,k,h=None):
        h=h or self.depth 
        return self.g*k**2*h
    def domega2_dk(self,k,h=None):
        h=h or self.depth
        return 2*self.g*k*h  
    def k_from_omega(self,omega,h=None):
        return omega/((self.g*h)**0.5)
 
class GravityCapillaryWaves(DispersionRelation):
    def __init__(self, g=9.81 | units.m/units.s**2,depth=None, 
                       fluid_density=1025. | units.kg/units.m**3,
                       surface_tension=0.074 | units.N/units.m):
        self.g=g
        self.depth=depth
        self.surface_tension=surface_tension
        self.fluid_density=fluid_density
        self.sigma_over_rho=surface_tension/fluid_density
    def omega2(self,k,h=None):
        h=h or self.depth
        return (self.g*k+self.sigma_over_rho*k**3)*numpy.tanh(k*h)
          
    def domega2_dk(self,k,h=None):
        h=h or self.depth
        th=numpy.tanh(k*h)
        return (self.g+3*self.sigma_over_rho*k**2)*th + \
               (self.g*k+self.sigma_over_rho*k**3)*h*(1-th**2)

class DeepGravityCapillaryWaves(DispersionRelation):
    def __init__(self, g=9.81 | units.m/units.s**2, 
                       fluid_density=1025. | units.kg/units.m**3,
                       surface_tension=0.074 | units.N/units.m):
        self.g=g
        self.surface_tension=surface_tension
        self.fluid_density=fluid_density
        self.sigma_over_rho=surface_tension/fluid_density
    def omega2(self,k,h=None):
        return (self.g*k+self.sigma_over_rho*k**3)          
    def domega2_dk(self,k,h=None):
        return (self.g+3*self.sigma_over_rho*k**2)

