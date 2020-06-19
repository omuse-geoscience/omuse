from math import isnan, sqrt
from amuse.test.amusetest import TestWithMPI

from omuse.community.iemic.interface import iemicInterface
from omuse.community.iemic.interface import iemic

from amuse.datamodel import Particle, Particles

kwargs=dict()

def newton(interface, x0, tol=1.e-7, maxit=1000):
    x = x0
    for k in range(maxit):
      fval = interface.rhs(x)
      interface.jacobian(x)
      dx = -interface.solve(fval)

      x = x + dx
 
      dxnorm = dx.norm()
      print("norm:", dxnorm)
      if dxnorm < tol:
        break

    return x

def setup_2dmoc(interface):
    interface.set_parameter('Ocean->Belos Solver->FGMRES tolerance', 1e-3)

    interface.set_parameter('Ocean->THCM->Global Bound xmin', 286.0)
    interface.set_parameter('Ocean->THCM->Global Bound xmax', 350.0)
    interface.set_parameter('Ocean->THCM->Global Bound ymin', -60.0)
    interface.set_parameter('Ocean->THCM->Global Bound ymax', 60.0)
    interface.set_parameter('Ocean->THCM->Depth hdim', 4000.0)
    interface.set_parameter('Ocean->THCM->Periodic', True)

    interface.set_parameter('Ocean->THCM->Global Grid-Size n', 3)
    interface.set_parameter('Ocean->THCM->Global Grid-Size m', 6)
    interface.set_parameter('Ocean->THCM->Global Grid-Size l', 6)

    interface.set_parameter('Ocean->THCM->Topography', 1)
    interface.set_parameter('Ocean->THCM->Flat Bottom', True)

    interface.set_parameter('Ocean->THCM->Starting Parameters->Combined Forcing', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Solar Forcing', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Salinity Forcing', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Wind Forcing', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Temperature Forcing', 10.0)

    interface.set_parameter('Ocean->THCM->Starting Parameters->P_VC', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Rossby-Number', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->CMPR', 0.0)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Horizontal Ekman-Number', 371.764)
    interface.set_parameter('Ocean->THCM->Starting Parameters->Rayleigh-Number', 15.6869)

    interface.set_parameter('Ocean->THCM->Rho Mixing', False)
    interface.set_parameter('Ocean->THCM->Coriolis Force', 0)
    interface.set_parameter('Ocean->THCM->Restoring Temperature Profile', 1)
    interface.set_parameter('Ocean->THCM->Restoring Salinity Profile', 0)
    interface.set_parameter('Ocean->THCM->Forcing Type', 1)
    interface.set_parameter('Ocean->THCM->Wind Forcing Type', 2)

def newtoncorrector(interface, par, ds, x, x0, l, l0, tol):
    # Set some parameters
    maxit = 20
    zeta = 1 / x.length()
    delta = 1e-6

    # Do the main iteration
    for k in range(maxit):
        # Set the parameter value and compute F (RHS of 2.2.9)
        interface.set_parameter(par, l)
        fval = interface.rhs(x)

        # Compute F_mu (bottom part of the RHS of 2.2.9)
        interface.set_parameter(par, l + delta)
        dflval = (interface.rhs(x) - fval) / delta
        interface.set_parameter(par, l)

        # Compute the jacobian at x
        interface.jacobian(x)

        # Solve twice with F_x (2.2.9)
        z1 = -interface.solve(fval)
        z2 = interface.solve(dflval)

        # Compute r (2.2.8)
        diff = x - x0
        rnp1 = zeta*diff.dot(diff) + (1-zeta)*(l-l0)**2 - ds**2

        # Compute dl (2.2.13)
        dl = (-rnp1 - 2*zeta*diff.dot(z1)) / (2*(1-zeta)*(l-l0) - 2*zeta*diff.dot(z2))

        # Compute dx (2.2.12)
        dx = z1 - dl*z2

        # Compute a new x and l (2.2.10 - 2.2.11)
        x = x + dx
        l = l + dl

        if dx.norm() < tol:
            print('Newton corrector converged in %d steps' % k)
            return (x, l)

    print('No convergence achieved by Newton corrector')

def continuation(interface, x0, par_name, target, ds, maxit):
    x = x0

    # Get the initial tangent (2.2.5 - 2.2.7). 'l' is called mu in Erik's thesis.
    delta = 1e-6
    l = interface.get_parameter(par_name)
    fval = interface.rhs(x)
    interface.set_parameter(par_name, l + delta)
    dl = (interface.rhs(x) - fval) / delta
    interface.set_parameter(par_name, l)

    # Compute the jacobian at x and solve with it (2.2.5)
    interface.jacobian(x)
    dx = -interface.solve(dl)

    # Scaling of the initial tangent (2.2.7)
    dl = 1
    zeta = 1 / x.length()
    nrm = sqrt(zeta * dx.dot(dx) + dl**2)
    dl = dl / nrm
    dx = dx / nrm

    dl0 = dl
    dx0 = dx

    # Perform the continuation
    for j in range(maxit):
        l0 = l
        x0 = x

        # Predictor (2.2.3)
        l = l0 + ds * dl0
        x = x0 + ds * dx0

        # Corrector (2.2.9 and onward)
        x2, l2 = newtoncorrector(interface, par_name, ds, x, x0, l, l0, 1e-4)

        print("%s:" % par_name, l2)

        if (l2 >= target and l0 < target) or (l2 <= target and l0 > target):
            # Converge onto the end point (we usually go past it, so we
            # use Newton to converge)
            l = target;
            interface.set_parameter(par_name, l);
            x = newton(interface, x, 1e-4);

            return x

        # Set the new values computed by the corrector
        dl = l2 - l0
        l = l2
        dx = x2 - x0
        x = x2

        if abs(dl) < 1e-10:
            return

        # Compute the tangent (2.2.4)
        dx0 = dx / ds
        dl0 = dl / ds

    return x

class iemicStateTests(TestWithMPI):
    def test1(self):
        instance = iemic(**kwargs)

        state0=instance.new_state()

        self.assertEqual(state0._id, 0)
        
        state1=instance.new_state()

        self.assertEqual(state1._id, 1)

        #~ del state0, state1

        instance.cleanup_code()
        instance.stop()


    def test2(self):
        instance = iemic(**kwargs)
                
        state=instance.new_state()
        
        del state
        
        self.assertEqual(instance.new_state()._id, 0)
        
        instance.cleanup_code()
        instance.stop()

    def test3(self):
        instance = iemic(**kwargs)
                
        x0=instance.new_state()
        
        sol=newton(instance, x0)
        
        instance.cleanup_code()
        instance.stop()

    def test4(self):
        instance = iemic(**kwargs)
                
        x=instance.new_state()
        
        x=x/2.
        
        print(x[:,:,:].u_velocity)
        
        instance.cleanup_code()
        instance.stop()

    def test5(self):
        instance = iemic(**kwargs)

        setup_2dmoc(instance)

        x = instance.get_state()
        x = newton(instance, x, 1e-10)

        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Combined Forcing', 1, 0.1, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->CMPR', -0.2, -0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.02, 0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->CMPR', 0, 0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.04, 0.5, 500)
        psi_min1, psi_max1 = instance.get_psi_m(x)
        self.assertAlmostEqual(psi_max1, 14.7, 1);
        self.assertAlmostEqual(psi_min1, 0);

        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.03, -0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.04, -0.5, 500)
        psi_min2, psi_max2 = instance.get_psi_m(x)
        self.assertAlmostEqual(psi_min1, -psi_max2, 4);
        self.assertAlmostEqual(psi_max1, -psi_min2, 4);

if __name__=="__main__":
    iemicStateTests().test5()
