import os
import sys
import numpy

from amuse.units import units

from amuse.test.amusetest import TestWithMPI

from interface import QGmodelInterface,QGmodel

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)


class TestQGmodelInterface(TestWithMPI):
    
    def test0(self):
        instance=QGmodelInterface()
        instance.initialize_code()
        instance.stop()

    def test1(self):
        instance=QGmodelInterface()
        instance.initialize_code()
        
        for key,val in [("Lx",4000000.),("Ly",4000000.),
         ("dx",10000.),("dy",10000.), ("dt",3600.),("H",4000.),
         ("rho",1000.), ("beta0",1.8616e-11),("tau",0.05),
         ("R_H",0),("A_H",100), ("lambda0",0.),("lambda1",2.e-5),
         ("Nm",1),("free_slip",1),("begin_time",0.),("wind_sigma",-99.),
         ("e111",0),("phi1z0",1.4142135623731),("ra_alpha",0.1)]:
            result,err=getattr(instance, 'get_'+key)()
            self.assertEquals( result, val)
            newvalue=type(val)(123.)
            err=getattr(instance, 'set_'+key)(newvalue)
            result,err=getattr(instance, 'get_'+key)()
            self.assertEquals( result,newvalue)

        instance.stop()

    def test2(self):
        instance=QGmodelInterface()
        instance.initialize_code()
        
        err=instance.set_Lx(100.)
        err=instance.set_dx(1)
        err=instance.set_Ly(1000.)
        err=instance.set_dy(20)

        instance.commit_parameters()

        Nx,err=instance.get_Nx()
        self.assertEquals(Nx,101)
        Ny,err=instance.get_Ny()
        self.assertEquals(Ny,51)


        instance.stop()

    def test3(self):
        instance=QGmodelInterface()
        instance.initialize_code()
        
        instance.commit_parameters()

        psi,err=instance.get_psi1_state(1,1,1)
        self.assertEqual(psi,0.)

        Nx,err=instance.get_Nx()
        Ny,err=instance.get_Ny()

        psi,err=instance.get_psi1_state(Nx,Ny,1)
        self.assertEqual(psi,0.)
         
        instance.stop()

    def test4(self):
        instance=QGmodelInterface(redirection="none")
        instance.initialize_code()
        
        instance.commit_parameters()

        instance.evolve_model(3600.)
        time,err=instance.get_time()
        self.assertEqual(time,3600.) 
        instance.stop()


class TestQGmodel(TestWithMPI):
    
    def test0(self):
        instance=QGmodel()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.evolve_model(1 | units.hour)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')        
        psi=instance.grid.psi
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()

    def test1(self):
        instance=QGmodel()
        instance.evolve_model(1 | units.hour)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')        
        psi=instance.grid.psi
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()

    def test2(self):
        instance=QGmodel(redirection="none")
        instance.grid[0,0,0].psi       
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.initialize_grid()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        
    def test3(self):
        instance=QGmodel(redirection="none")
        dt=0.123 | units.hour
        instance.parameters.dt=dt
        instance.evolve_model(2*dt)
        self.assertEqual(instance.model_time,2*dt)
        instance.stop()

    def test4(self):
        instance=QGmodel(redirection="none")
        instance.grid[0,0,0].dpsi_dt       
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()

        
    def test5(self):
        instance=QGmodel(redirection="none")
        dt=instance.parameters.dt
        instance.evolve_model(dt)

        dpsi_dt1=instance.grid.dpsi_dt
        psi1=instance.grid.psi
        instance.evolve_model(2*dt)

        psi2=instance.grid.psi
        dpsi_dt2=instance.grid.dpsi_dt
        
        dpsi_estimate=((psi2-psi1)/dt).mean()
        dpsi=0.5*(dpsi_dt1+dpsi_dt2).mean()
        self.assertTrue( abs(dpsi-dpsi_estimate)/dpsi < 0.1)
        
        instance.stop()
        
        

        
