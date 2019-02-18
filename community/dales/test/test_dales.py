from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.dales.interface import DalesInterface
from omuse.community.dales.interface import Dales

from omuse.units import units

default_options={}
#default_options=dict(number_of_workers=4)
#default_options=dict(channel_type="sockets",redirection="none")
#default_options=dict(redirection="none",debugger="gdb")

import logging
import pickle
import numpy
import numpy.random

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from nose.tools import nottest

class TestDalesInterface(TestWithMPI):

    def test1(self):

        print "Test 1: instantiate and clean up"

        instance = Dales(**default_options)
        instance.cleanup_code()
        instance.stop()

    def test2(self):

        print "Test 2: instantiate with 4 workers and clean up"

        instance = Dales(number_of_workers=4,**default_options)
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):

        print "Test 3: instantiate, run one minute and clean up"

        instance = Dales(**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (1 | units.minute))
        newtim=instance.get_model_time()
        self.assertTrue(newtim>=tim)
        instance.cleanup_code()
        instance.stop()

    def test4(self):

        print "Test 4: instantiate with 4 workers, run one minute and clean up"

        instance = Dales(number_of_workers=4,**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (1 | units.minute))
        newtim=instance.get_model_time()
        self.assertTrue(newtim>=tim)
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):

        print "Test 5: instantiate, run 10 seconds and retrieve temperature profile"

        instance = Dales(number_of_workers=4,**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (10 | units.s))
        profile=instance.get_profile_T_(numpy.arange(1,96))
        print "The retrieved profile is:",profile
        instance.cleanup_code()
        instance.stop()

    def test6(self):
        print "Test 6: instantiate, run 2 seconds and retrieve profiles"

        instance = Dales(number_of_workers=4,**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()
        instance.evolve_model(tim + (2 | units.s))

        
        print "The retrieved U profile is:", instance.get_profile_U()
        print "The retrieved V profile is:", instance.get_profile_V()
        print "The retrieved W profile is:", instance.get_profile_W()
        print 

        print "The retrieved THL profile is:", instance.get_profile_THL()
        print "The retrieved QT profile is:", instance.get_profile_QT()
        print "The retrieved zf levels:", instance.get_zf()
        print "The retrieved zh levels:", instance.get_zh()

        i,j,k,xsize,ysize = instance.get_params_grid()
        print "Grid size", (i, j, k)
        print "Horizontal extent of model", (xsize,ysize)

        instance.cleanup_code()
        instance.stop()

    def test7(self):
        print "Test 7: instantiate, retrieve profiles without time stepping for different number of threads"
        
        profile = numpy.loadtxt("prof.inp.001")
        zf_in  = profile[:,0]
        thl_in = profile[:,1]
        qt_in  = profile[:,2]

        # root mean square difference between vectors a and b
        def rms(a, b):
            return numpy.sqrt(((a-b)**2).mean())
        
        for n in (1,2,4):
            print "---------------------------------------------------- %d threads -------------------------------------------------"%n
            instance = Dales(number_of_workers=n,**default_options)
            tim=instance.get_model_time()
            instance.commit_grid()

            zf  = instance.get_zf()
            thl = instance.get_profile_THL()
            qt  = instance.get_profile_QT()
            for i in range(0,instance.k):
                print "%3d %6.1f - %7.3f %7.3f - %8.4f %8.4f"%(i, zf[i].value_in(units.m), thl_in[i], thl[i].value_in(units.K), qt_in[i], qt[i].value_in(units.shu))
            # get rid of units
            thl = [t.value_in(units.K) for t in thl]
            rmsT  = rms(thl, thl_in)
            rmsQT = rms(qt, qt_in)
            print n, ' threads:'
            print 'RMS error in thl', rmsT
            print 'RMS error in qt', rmsQT

            # demand that the rms error btw input profiles and extracted average profiles is small
            # DALES adds some random fluctuations, so the average != input
            assert rmsT < 1e-2 
            assert rmsQT < 1e-6
            
            #print "The retrieved THL profile is:", thl
            #print "The retrieved QT profile is:", qt
            #print "The retrieved zf levels:", instance.get_zf()
            #print "The retrieved zh levels:", instance.get_zh()
            
            i,j,k,xsize,ysize = instance.get_params_grid()
            print "Grid size", (i, j, k)
            print "Horizontal extent of model", (xsize,ysize)

            instance.cleanup_code()
            instance.stop()
            print
            
    def test8(self):
        print "Test 8: instantiate, set a forcing, run 30 seconds and retrieve profiles"
        # the effect of the forcing is seen in the U profile at the same height as the forcing was places
        # no automatic test for this

        instance = Dales(number_of_workers=4,**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()

        #print "Number of height levels:", instance.k
        #tend = numpy.arange(instance.k)*.01
        #print ('U_tend in python', tend)
        #instance.set_tendency_U(tend)

        tend = numpy.zeros(instance.k)
        tend[5] = -.1
        tend[6] = -.2   # layer 9 sees the most effect of this
        tend[7] = -.1

        instance.set_tendency_U(tend | units.m/units.s**2)
        
        instance.evolve_model(tim + (30 | units.s))

        U = instance.get_profile_U()
        V = instance.get_profile_V()
        W = instance.get_profile_W()
        THL = instance.get_profile_THL()
        QT = instance.get_profile_QT()

        print "The retrieved U profile is:", U
        print "The retrieved V profile is:", V
        #print "The retrieved W profile is:", W
        print 

        #print "The retrieved THL profile is:", THL
        print "The retrieved QT profile is:", QT
        #print "The retrieved zf levels:", instance.get_zf()
        #print "The retrieved zh levels:", instance.get_zh()

        

        
        i,j,k,xsize,ysize = instance.get_params_grid()
        print "Grid size", (i, j, k)
        print "Horizontal extent of model", (xsize,ysize)

        instance.cleanup_code()
        instance.stop()

        

    def test9(self):
        print "Test 9: get 3D data and compare with slab averages"

        # root mean square difference between vectors a and b
        def rms(a, b):
            return (((a-b)**2).mean())**0.5
        
        instance = Dales(number_of_workers=2,**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()

        tend = numpy.zeros(instance.k)
        tend[5] = -1
        instance.set_tendency_U(tend | units.m/units.s**2)
        
        instance.evolve_model(tim + (180 | units.s))
        
        def getU(imin=0, imax=instance.itot, jmin=0, jmax=instance.jtot, kmin=0, kmax=instance.k):
            points = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
            points = points.reshape(3, -1)
            i,j,k = points
            print 'i', i
            print 'j', j
            print 'k', k
        
            field = instance.get_field_U(i,j,k)
            field = field.reshape ((imax-imin, jmax-jmin, kmax-kmin))
            return field
            
        
        #field = field.reshape ((isize, jsize, ksize))

        #print 'vertical profile U'
        #print getU(imin=3,imax=4,jmin=8,jmax=9) # a vertical profile

        #print 'horizontal profile U at k=5'
        #print getU(imin=1,imax=10, jmin=8,jmax=9, kmin=5, kmax=6) # a vertical profile


        #print 'horizontal profile U at k=6'
        #print getU(imin=1,imax=10, jmin=8,jmax=9, kmin=6, kmax=7) # a vertical profile

        # get slab averages
        U = instance.get_profile_U()
        V = instance.get_profile_V()
        W = instance.get_profile_W()
        THL = instance.get_profile_THL()
        QT = instance.get_profile_QT()

        # get 3D profiles
        U3   = instance.get_field('U')
        V3   = instance.get_field('V')
        W3   = instance.get_field('W')
        THL3 = instance.get_field('THL')
        QT3  = instance.get_field('QT')

        
        # calculate own slab averages of 3D profiles
        U3a   = U3.mean(axis=0).mean(axis=0)
        V3a   = V3.mean(axis=0).mean(axis=0)
        W3a   = W3.mean(axis=0).mean(axis=0)
        THL3a = THL3.mean(axis=0).mean(axis=0)
        QT3a  = QT3.mean(axis=0).mean(axis=0)

        # rms difference btw retrieved and calculated slab averages
        rmsU    = rms(U,   U3a)
        rmsV    = rms(V,   V3a)
        rmsW    = rms(W,   W3a)
        rmsTHL  = rms(THL, THL3a)
        rmsQT   = rms(QT,  QT3a)
        
        print 'RMS errors (U, V, W, THL, QT)', rmsU, rmsV, rmsW, rmsTHL, rmsQT
        print THL3 [0, 0, 0:4]
#        print instance.get_field('THL', imin)
        
        assert rmsU   < 1e-12 | units.m / units.s
        assert rmsV   < 1e-12 | units.m / units.s
        assert rmsW   < 1e-12 | units.m / units.s
        assert rmsTHL < 1e-12 | units.K
        assert rmsQT  < 1e-12
        
        instance.cleanup_code()
        instance.stop()


    def test10(self):
        fileName = 'reference-data/test10.pkl'
        save = False #True
        
        print "Test 10: get 3D data, compare with stored reference"

        # summed root mean square difference between vectors a and b
        def rmsSum(a, b):
            return (((a-b)**2).sum())**0.5
        
        instance = Dales(number_of_workers=1,**default_options)
        tim=instance.get_model_time()
        instance.commit_grid()

       # tend = numpy.zeros(instance.k)
       # tend[5] = -1
       # instance.set_tendency_U(tend)
        
        instance.evolve_model(tim + (300 | units.s), exactEnd=0)

        # get 3D profiles
        U3   = instance.get_field('U')
        V3   = instance.get_field('V')
        W3   = instance.get_field('W')
        THL3 = instance.get_field('THL')
        QT3  = instance.get_field('QT')

        if save:
            outFile = open(fileName, 'wb')
            pkl = (U3, V3, W3, THL3, QT3)
            pickle.dump(pkl, outFile)
            print('Saved reference profiles in', fileName)
        else:
            inFile = open(fileName, 'rb')
            pkl = pickle.load(inFile)
            U3p, V3p, W3p, THL3p, QT3p = pkl
            rmsU   = rmsSum(U3, U3p)
            rmsV   = rmsSum(V3, V3p)
            rmsW   = rmsSum(W3, W3p)
            rmsTHL = rmsSum(THL3, THL3p)
            rmsQT  = rmsSum(QT3, QT3p)
            
            print 'RMS errors (U, V, W, THL, QT)', rmsU, rmsV, rmsW, rmsTHL, rmsQT
        
            assert rmsU   < 1e-9 | units.m / units.s
            assert rmsV   < 1e-9 | units.m / units.s
            assert rmsW   < 1e-9 | units.m / units.s
            assert rmsTHL < 1e-9 | units.K
            assert rmsQT  < 1e-9
        
        instance.cleanup_code()
        instance.stop()



    def test11(self):
        print("Test 11: set 3D data, get it back, compare")

        fieldunits = {
                'U'   : units.m / units.s,
                'V'   : units.m / units.s,
                'W'   : units.m / units.s,
                'THL' : units.K,
                'QT'  : units.shu,
            }

        # write a single grid point at i,j,k, read back
        def testSingle(i,j,k):
            #a = numpy.random.random((1,1,1)) # 3D matrix with one element
            a = numpy.random.random()  | fieldunits[field]# scalar
            b = dales.get_field(field, i, i+1, j, j+1, k, k+1)
#            b = drop_units(field, b)
            dales.set_field(field, a, i, j, k)
            c = dales.get_field(field, i, i+1, j, j+1, k, k+1)
#            c = drop_units(field, c)
            #~ print field, (i,j,k), 'was:', b, 'set: ', a, 'get:', c
            self.assertEquals(a,c)

        # write to a block [i:i+si, j:j+sj, k:k+sk], read back
        def testBlock(i,j,k, si, sj, sk):
            a = numpy.random.random((si,sj,sk)) | fieldunits[field]
            b = dales.get_field(field, i, i+si, j, j+sj, k, k+sk)
#            b = drop_units(field, b)
            dales.set_field(field, a, i, j, k)
            c = dales.get_field(field, i, i+si, j, j+sj, k, k+sk)
#            c = drop_units(field, c)
            print field, (i,j,k),  (si,sj,sk) #'was:', b, 'set: ', a, 'get:', c
            self.assertEquals(a,c)
            
        for n in (1,2,4):
            print ' --- ', n, 'dales workers --- '
            dales = Dales(number_of_workers=n,**default_options)
            dales.commit_grid()

            # set random grid elements for all fields, read back
            for s in range(0,50):
                for field in ('U', 'V', 'W', 'THL', 'QT'):
                    i = numpy.random.randint(0,dales.itot)
                    j = numpy.random.randint(0,dales.jtot)
                    k = numpy.random.randint(0,dales.k)
                    testSingle(i, j, k)

            for s in range(0,50):
                for field in ('U', 'V', 'W', 'THL', 'QT'):
                    i = numpy.random.randint(0,dales.itot)
                    j = numpy.random.randint(0,dales.jtot)
                    k = numpy.random.randint(0,dales.k)

                    si = numpy.random.randint(1,dales.itot-i+1)
                    sj = numpy.random.randint(1,dales.jtot-j+1)
                    sk = numpy.random.randint(1,dales.k-k+1)
                    testBlock(i, j, k, si, sj, sk)
            


            
            dales.cleanup_code()
            dales.stop()

class TestDales(TestWithMPI):

    def test1(self):

        print "Test 1: instantiate and stop"

        instance = Dales(**default_options)
        instance.stop()

    def test2(self):

        print "Test 2: instantiate, run one minute and clean up"

        instance = Dales(**default_options)
        tim=instance.model_time
        instance.evolve_model(tim + (1 | units.minute))
        newtim=instance.model_time
        self.assertTrue(newtim>=tim)
        instance.stop()

    def test3(self):

        print "Test 3: instantiate, run one minute and clean up"

        instance = Dales(**default_options)
        instance.parameters.exactEndFlag=True
        tim=instance.model_time
        instance.evolve_model(tim + (1 | units.minute))
        newtim=instance.model_time
        self.assertTrue(newtim-tim==(1 | units.minute))
        instance.stop()

    def xtest4(self):

        print "Test 3: instantiate, run one minute and clean up"

        instance = Dales(**default_options)
        print instance.get_name_of_current_state()
        print instance.parameters
        print instance.get_name_of_current_state()
        instance.parameters.input_file="dummy"
        print instance.get_name_of_current_state()
        print instance.parameters
        print instance.get_name_of_current_state()
