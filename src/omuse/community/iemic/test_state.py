from amuse.test.amusetest import TestWithMPI

from omuse.community.iemic.interface import iemicInterface
from omuse.community.iemic.interface import iemic

from omuse.community.iemic.implicit_utils import newton, continuation, time_stepper

from amuse.datamodel import Particle, Particles

kwargs=dict()


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

class iemicStateTests(TestWithMPI):
    _multiprocess_can_split_ = True

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

        # Setup a simplified 2D AMOC problem
        setup_2dmoc(instance)

        # Converge to an initial steady state
        x = instance.get_state()
        x = newton(instance, x, 1e-10)

        # Check that we can perform a continuation in multiple parameters
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Combined Forcing', 1, 0.1, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->CMPR', -0.2, -0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.02, 0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->CMPR', 0, 0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.04, 0.5, 500)
        psi_min1, psi_max1 = instance.get_psi_m(x)
        self.assertAlmostEqual(psi_max1, 14.7, 1)
        self.assertAlmostEqual(psi_min1, 0)

        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.03, -0.5, 500)
        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.035, -0.5, 500)
        psi_min3, psi_max3 = instance.get_psi_m(x)
        x3 = x

        x = continuation(instance, x, 'Ocean->THCM->Starting Parameters->Salinity Forcing', 0.04, 0.5, 500)
        psi_min2, psi_max2 = instance.get_psi_m(x)
        self.assertAlmostEqual(psi_min1, -psi_max2, 4)
        self.assertAlmostEqual(psi_max1, -psi_min2, 4)

        b = instance.rhs(x)
        self.assertAlmostEqual(b.norm(), 0, 3)

        instance.set_parameter('Ocean->THCM->Starting Parameters->Salinity Forcing', 0.035)

        b = instance.rhs(x)
        self.assertGreaterEqual(b.norm(), 0.01)
        self.assertGreaterEqual(x.norm(), 0.1)

        # Check that the time stepper converges to the same steady state,
        # but not in a single time step
        x = time_stepper(instance, x, 0.5, 1, 1)

        b = instance.rhs(x)
        self.assertGreaterEqual(b.norm(), 0.01)
        self.assertGreaterEqual(x.norm(), 0.1)

        x = time_stepper(instance, x, 0.5, 30, 3000)

        b = instance.rhs(x)
        self.assertAlmostEqual(b.norm(), 0, 3)
        self.assertGreaterEqual(x.norm(), 0.1)

        psi_min4, psi_max4 = instance.get_psi_m(x)
        self.assertAlmostEqual(psi_min3, psi_min4, 3)
        self.assertAlmostEqual(psi_max3, psi_max4, 3)

        x3 -= x
        self.assertAlmostEqual(x3.norm(), 0, 1)

        instance.cleanup_code()
        instance.stop()

    def test6(self):
        instance = iemic(**kwargs)

        # Setup a simplified 2D AMOC problem
        setup_2dmoc(instance)

        # Converge to an initial steady state
        x = instance.get_state()
        x = newton(instance, x, 1e-10)
        
        instance.set_parameter('Ocean->THCM->Starting Parameters->Combined Forcing', 1.)
        instance.set_parameter('Ocean->THCM->Starting Parameters->Salinity Forcing', .035)

        self.assertEqual(instance.parameters.Ocean__THCM__Starting_Parameters__Salinity_Forcing,0.035)

        instance.evolve_model(100)

        instance.cleanup_code()
        instance.stop()

    def test7(self):
        instance = iemic(redirection="none")
        print(instance.parameters)

        x = instance.get_state()
        
        lm=instance.get_surface_mask()
        print(lm)
        
        instance.cleanup_code()
        instance.stop()
        


if __name__=="__main__":
    iemicStateTests().test7()
