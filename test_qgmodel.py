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
         ("Nm",1),("begin_time",0.),("wind_sigma",-99.),
         ("e111",0),("phi1z0",1.4142135623731),("ra_alpha",0.1)]:
            result,err=getattr(instance, 'get_'+key)()
            self.assertEquals( result, val)
            newvalue=type(val)(123.)
            err=getattr(instance, 'set_'+key)(newvalue)
            result,err=getattr(instance, 'get_'+key)()
            self.assertEquals( result,newvalue)

        instance.stop()

    def test1b(self):
        instance=QGmodelInterface()
        instance.initialize_code()
        
        lowx,highx,lowy,highy,err= instance.get_boundary_conditions()
        for x in [lowx,highx,lowy,highy]:
          self.assertEqual(x,"free_slip")

        instance.set_boundary_conditions("no_slip","interface","interface","no_slip")
        lowx,highx,lowy,highy,err= instance.get_boundary_conditions()
        self.assertEqual(lowx,"no_slip")
        self.assertEqual(highy,"no_slip")
        self.assertEqual(lowy,"interface")
        self.assertEqual(highx,"interface")

        # check for failure
        err=instance.set_boundary_conditions("no_slip","interface","interface","wrong")
        self.assertEqual(err,8)

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

    def test5(self):
        instance=QGmodelInterface()
        instance.initialize_code()
        instance.commit_parameters()
        Nx,err=instance.get_Nx()
        Ny,err=instance.get_Ny()
        Nm,err=instance.get_Nm()
        
        minx,maxx,miny,maxy,minz,maxz,err=instance.get_index_range_inclusive()
        self.assertEqual(minx,1)
        self.assertEqual(miny,1)
        self.assertEqual(minz,1)
        self.assertEqual(maxx,Nx)
        self.assertEqual(maxy,Ny)
        self.assertEqual(maxz,Nm)
        
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
        
    def test6(self):
        instance=QGmodel()
        Nx=instance.get_Nx()
        Ny=instance.get_Ny()
        Nm=instance.get_Nm()

        self.assertEqual(401,Nx)
        
        minx,maxx,miny,maxy,minz,maxz=instance.get_index_range_inclusive()
        self.assertEqual(minx,1)
        self.assertEqual(miny,1)
        self.assertEqual(minz,1)
        self.assertEqual(maxx,Nx)
        self.assertEqual(maxy,Ny)
        self.assertEqual(maxz,Nm)
        
        instance.stop()
        
    def test7(self):
        instance=QGmodel()
        Nx=instance.get_Nx()
        Ny=instance.get_Ny()
        Nm=instance.get_Nm()
        self.assertEqual(401,Nx)
        
        minx,maxx,miny,maxy,minz,maxz=instance.get_boundary_index_range_inclusive(1)
        self.assertEqual(minx,0)
        self.assertEqual(miny,0)
        self.assertEqual(minz,1)
        self.assertEqual(maxx,1)
        self.assertEqual(maxy,Ny+1)
        self.assertEqual(maxz,Nm)

        minx,maxx,miny,maxy,minz,maxz=instance.get_boundary_index_range_inclusive(2)
        self.assertEqual(minx,Nx)
        self.assertEqual(miny,0)
        self.assertEqual(minz,1)
        self.assertEqual(maxx,Nx+1)
        self.assertEqual(maxy,Ny+1)
        self.assertEqual(maxz,Nm)

        minx,maxx,miny,maxy,minz,maxz=instance.get_boundary_index_range_inclusive(3)
        self.assertEqual(minx,0)
        self.assertEqual(miny,0)
        self.assertEqual(minz,1)
        self.assertEqual(maxx,Nx+1)
        self.assertEqual(maxy,1)
        self.assertEqual(maxz,Nm)        

        minx,maxx,miny,maxy,minz,maxz=instance.get_boundary_index_range_inclusive(4)
        self.assertEqual(minx,0)
        self.assertEqual(miny,Ny)
        self.assertEqual(minz,1)
        self.assertEqual(maxx,Nx+1)
        self.assertEqual(maxy,Ny+1)
        self.assertEqual(maxz,Nm)        
        
        instance.stop()

    def test8(self):
        instance=QGmodel()
        instance.parameters.Ly=instance.parameters.Lx/2
        Nx=instance.get_Nx()
        Ny=instance.get_Ny()
        Nm=instance.get_Nm()
        s=instance.boundaries(1).shape
        self.assertEqual(s,(2,Ny+2,1))
        s=instance.boundaries(2).shape
        self.assertEqual(s,(2,Ny+2,1))
        s=instance.boundaries(3).shape
        self.assertEqual(s,(Nx+2,2,1))
        s=instance.boundaries(4).shape
        self.assertEqual(s,(Nx+2,2,1))

        instance.stop()        
        
    def test9(self):
        instance=QGmodel()
        org=instance.grid
        copy=org.copy()
        instance.stop()        

    def test10(self):
        instance=QGmodel()
        org=instance.boundaries(1)
        self.assertEquals( org.psi.number,0.)
        self.assertEquals( org.dpsi_dt.number,0.)
        shape=org.shape
        copy=org.copy()
        instance.stop()

    def test11(self):
        instance=QGmodel(redirection="none")
        instance.parameters.Ly=instance.parameters.Ly/8
        instance.parameters.Lx=instance.parameters.Lx/8
  
        for i in [1,2,3,4]:
          org=instance.boundaries(i)          
          shape=org.shape
          copy=org.copy()
          channel=copy.new_channel_to(org)
          copy.psi=numpy.random.random(shape) | units.m**2/units.s
          copy.dpsi_dt=numpy.random.random(shape) | units.m**2/units.s**2
  
          channel.copy_attributes(["psi","dpsi_dt"])
          self.assertEquals(copy.psi, instance.boundaries(i).psi)
          self.assertEquals(copy.dpsi_dt, instance.boundaries(i).dpsi_dt)
        
      
        instance.stop()

    def test12(self):
        instance=QGmodel()
        instance.parameters.Ly=instance.parameters.Ly/16
        instance.parameters.Lx=instance.parameters.Lx/16

        Nx=instance.parameters.Nx
        Ny=instance.parameters.Ny

        ix=numpy.arange(Nx)+1
        iy=numpy.arange(Ny)+1
        onesx=numpy.ones(Nx)
        onesy=numpy.ones(Ny)

        i1= numpy.outer(ix,onesx).flatten()
        j1= numpy.outer(onesy,iy).flatten()

        x=instance.grid.x.flatten()
        y=instance.grid.y.flatten()
        i,j=instance.get_index_of_position(x,y)

        self.assertEquals(i-i1,0.)
        self.assertEquals(j-j1,0.)

        instance.stop()        

    def test13(self):
        instance=QGmodel()
        instance.parameters.xbound1="interface"
        instance.commit_parameters()
        b1,b2,b3,b4=instance.get_boundary_conditions()
        self.assertEqual(b1,"interface")
        instance.stop()    

    def test14(self):
        instance=QGmodel()
        instance.parameters.ybound1="interface"
        instance.commit_parameters()
        b1,b2,b3,b4=instance.get_boundary_conditions()
        self.assertEqual(b3,"interface")
        instance.stop()    
        
        
