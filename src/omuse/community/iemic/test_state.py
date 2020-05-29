from math import isnan
from amuse.test.amusetest import TestWithMPI

from omuse.community.iemic.interface import iemicInterface
from omuse.community.iemic.interface import iemic

from amuse.datamodel import Particle, Particles

kwargs=dict()

def newton(interface, x0, tol=1.e-7, maxit=1000):

    x = x0

    k=0
    while k<maxit:
      fval = interface.rhs(x)
      interface.jacobian(x)
      dx = -1*interface.solve(fval)

      x = x + dx
 
      dxnorm=dx.norm()
      print("norm:", dxnorm)
      if dxnorm < tol:
        break
        
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
        
        self.assertEqual(instance.new_state()._id, 1)
        
        instance.cleanup_code()
        instance.stop()

    def test3(self):
        instance = iemic(**kwargs)
                
        x0=instance.new_state()
        
        sol=newton(instance, x0)
        
        instance.cleanup_code()
        instance.stop()

                
if __name__=="__main__":
    iemicStateTests().test3()
