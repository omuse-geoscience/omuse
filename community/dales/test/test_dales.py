from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.dales.interface import Dales

from omuse.units import units

default_options = {}
# default_options=dict(number_of_workers=4)
# default_options=dict(channel_type="sockets",redirection="none")
# default_options=dict(redirection="none",debugger="gdb")

import logging
import pickle
import numpy
import numpy.random

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)


def cleanup_data(rundir):
    if os.path.isdir(rundir):
        for f in os.listdir(rundir):
            os.remove(os.path.join(rundir, f))
        os.rmdir(rundir)
    elif os.path.isfile(rundir):
        os.remove(rundir)


class TestDalesInterface(TestWithMPI):

    def test_namopt_file_is_written(self):
        rundir = "work0"
        os.mkdir(rundir)
        os.chdir(rundir)
        try:
            instance = Dales(redirection="none")
            instance.commit_parameters()
            assert os.path.isfile("namoptions.001")
            instance.cleanup_code()
            instance.stop()
        finally:
            for f in os.listdir("."):
                os.remove(f)
            os.chdir("..")
            os.rmdir(rundir)

    def test_prof_file_is_written(self):
        rundir = "work1"
        os.mkdir(rundir)
        os.chdir(rundir)
        try:
            instance = Dales(redirection="none")
            instance.commit_parameters()
            assert os.path.isfile("prof.inp.001")
            instance.cleanup_code()
            instance.stop()
        finally:
            for f in os.listdir("."):
                os.remove(f)
            os.chdir("..")
            os.rmdir(rundir)

    def test_lscale_file_is_written(self):
        rundir = "work2"
        os.mkdir(rundir)
        os.chdir(rundir)
        try:
            instance = Dales(redirection="none")
            instance.commit_parameters()
            assert os.path.isfile("lscale.inp.001")
            instance.cleanup_code()
            instance.stop()
        finally:
            for f in os.listdir("."):
                os.remove(f)
            os.chdir("..")
            os.rmdir(rundir)

    def test_set_workdir(self):
        rundir = "work3"
        instance = Dales(workdir=rundir, redirection="none")
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_namopt_file_written_in_workdir(self):
        rundir = "work4"
        instance = Dales(workdir=rundir, redirection="none")
        instance.commit_parameters()
        assert os.path.isfile(os.path.join(rundir, "namoptions.001"))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_load_rico(self):
        rundir = "work-rico"
        instance = Dales(case="rico", workdir=rundir, redirection="none")
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_rico(self):
        rundir = "work-rico2"
        instance = Dales(case="rico", workdir=rundir, redirection="none")
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_load_bomex(self):
        rundir = "work-bomex"
        instance = Dales(case="bomex", workdir=rundir, redirection="none")
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_bomex(self):
        rundir = "work-bomex2"
        instance = Dales(case="bomex", workdir=rundir, redirection="none")
        tim = instance.get_model_time()
        instance.evolve_model(tim + (1 | units.minute))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_load_sp_case(self):
        rundir = "work-sp"
        instance = Dales(case="sp-testcase", workdir=rundir, redirection="none")
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_get_grid_dimensions(self):
        rundir = "work-sp2"
        instance = Dales(case="sp-testcase", workdir=rundir, redirection="none")
        instance.commit_parameters()
        assert (instance.get_itot(), instance.get_jtot(), instance.get_ktot()) == (200, 200, 160)
        assert (instance.get_dx().value_in(units.m), instance.get_dy().value_in(units.m)) == (200, 200)
        assert numpy.array_equal(instance.get_zf().value_in(units.m), numpy.arange(12.5, 4000., 25.))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_sp_case(self):
        rundir = "work-sp3"
        instance = Dales(case="sp-testcase", workdir=rundir, redirection="none")
        tim = instance.get_model_time()
        instance.evolve_model(tim + (1 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_upscale_arm_brown_case(self):
        rundir = "work-arm-brown"
        instance = Dales(case="arm_brown", workdir=rundir, redirection="none")
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        instance.commit_parameters()
        assert os.path.isfile(os.path.join(rundir, "ls_flux.inp.001"))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_grid_shapes(self):
        rundir = "work-sp4"
        instance = Dales(case="sp-testcase", workdir=rundir, redirection="none")
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        tim = instance.get_model_time()
        instance.evolve_model(tim + (1 | units.s))
        temp_profile = instance.profiles.T
        assert temp_profile.shape == (instance.get_ktot())
        v_block = instance.grid.V
        assert v_block.shape == (instance.get_itot(), instance.get_jtot(), instance.get_ktot())
        twp_slab = instance.surface_field.TWP
        assert twp_slab.shape == (instance.get_itot(), instance.get_jtot())
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_u_profile_bomex(self):
        rundir = "work-bomex3"
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (5 | units.s))
        u_profile = instance.profiles.U.value_in(units.m / units.s)
        u_field = instance.grid.U.value_in(units.m / units.s)
        print u_profile
        print numpy.mean(u_field, axis=(0, 1))
        assert numpy.allclose(u_profile, numpy.mean(u_field, axis=(0, 1)), rtol=1.e-9)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_qt_profile_sp(self):
        rundir = "work-sp5"
        instance = Dales(case="sp-testcase", workdir=rundir, number_of_workers=4, redirection="none")
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        tim = instance.get_model_time()
        instance.evolve_model(tim + (5 | units.s))
        q_profile = instance.profiles.QT.value_in(units.mfu)
        q_field = instance.grid.QT.value_in(units.mfu)
        assert numpy.allclose(q_profile, numpy.mean(q_field, axis=(0, 1)), rtol=1.e-9)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_set_bomex_t_field(self):
        rundir = "work-bomex4"
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        thlfld = numpy.copy(instance.grid.THL.value_in(units.K))
        thlfld[:, :, 1] = 304.
        instance.grid.THL = thlfld | units.K
        instance.evolve_model(tim + (60 | units.s))
        thl1 = numpy.mean(instance.grid.THL.value_in(units.K)[:, :, 1], axis=(0, 1))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (70 | units.s))
        thl2 = numpy.mean(instance.grid.THL.value_in(units.K)[:, :, 1], axis=(0, 1))
        assert thl2 < thl1 < 304.

    def test_bomex_3d_grid_interface(self):
        rundir = "work-bomex45"
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (2 | units.s))
        qblock = numpy.full((3, 3, 1), 1.234e-5) | units.mfu
        instance.set_field("QT", qblock, imin=8, jmin=4, kmin=10)
        qt = instance.get_field("QT", imin=8, imax=11, jmin=4, jmax=7, kmin=10, kmax=11)
        qtgrid = instance.grid[7:10, 3:6, 9:10].QT.value_in(units.mfu)
        assert numpy.array_equal(qblock, qt.value_in(units.mfu))
        assert numpy.array_equal(qtgrid, qt.value_in(units.mfu))

    def test_set_bomex_t_value(self):
        rundir = "work-bomex5"
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        instance.grid[1::2, 0::2, 1].THL = 302. | units.K
        instance.grid[0::2, 1::2, 1].THL = 302. | units.K
        instance.grid[0::2, 0::2, 1].THL = 298. | units.K
        instance.grid[1::2, 1::2, 1].THL = 298. | units.K
        thl1 = instance.grid.THL.value_in(units.K)[5, 6, 1]
        thl2 = instance.grid.THL.value_in(units.K)[8, 7, 1]
        thl3 = instance.grid.THL.value_in(units.K)[4, 2, 1]
        thl4 = instance.grid.THL.value_in(units.K)[3, 3, 1]
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)
        assert thl1 == thl2 == 302.
        assert thl3 == thl4 == 298.

    def test_set_bomex_q_forcing(self):
        rundir = "work-bomex6"
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        zf = instance.get_zf().value_in(units.km)
        instance.forcings[numpy.argwhere(zf < 2.2)].tendency_QT = 1.e-5 | (units.mfu / units.s)
        instance.forcings[numpy.argwhere(zf >= 2.2)].tendency_QT = 1.e-8 | (units.mfu / units.s)
        assert instance.forcings.tendency_QT[0].value_in(units.mfu / units.s) == 1.e-5
        assert instance.forcings.tendency_QT[-1].value_in(units.mfu / units.s) == 1.e-8
        tim = instance.get_model_time()
        instance.evolve_model(tim + (60 | units.s))
        assert instance.profiles.QT[0].value_in(units.mfu) > instance.profiles.QT[-1].value_in(units.mfu)
        instance.stop()
        cleanup_data(rundir)

    def test_set_bomex_q_nudging(self):
        rundir = "work-bomex7"
        instance = Dales(case="bomex", workdir=rundir)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        zf = instance.get_zf().value_in(units.km)
        instance.nudging[numpy.argwhere(zf < 2.2)].QT = 1.e-3 | units.mfu
        instance.nudging[numpy.argwhere(zf >= 2.2)].QT = 1.e-6 | units.mfu
        assert instance.nudging.QT[0].value_in(units.mfu) == 1.e-3
        assert instance.nudging.QT[-1].value_in(units.mfu) == 1.e-6
        tim = instance.get_model_time()
        instance.evolve_model(tim + (60 | units.s))
        assert instance.profiles.QT[0].value_in(units.mfu) > instance.profiles.QT[-1].value_in(units.mfu)
        instance.stop()
        cleanup_data(rundir)

    def test_atex_surface_fields(self):
        rundir = "work-atex"
        instance = Dales(case="atex", workdir=rundir)
        instance.parameters_DOMAIN.itot = 32
        instance.parameters_DOMAIN.jtot = 32
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        instance.set_wt_surf(0.1 | units.m * units.s ** -1 * units.K)
        wq = instance.get_wt_surf().value_in(units.m * units.s ** -1 * units.K)
        assert(wq == 0.1)
        instance.set_z0h_surf(0.2 | units.m)
        instance.evolve_model(tim + (1 | units.s))
        z0 = instance.get_z0h_surf().value_in(units.m)
        instance.stop()
        cleanup_data(rundir)

#     def test7(self):
#         print "Test 7: instantiate, retrieve profiles without time stepping for different number of threads"
#
#         profile = numpy.loadtxt("prof.inp.001")
#         zf_in  = profile[:,0]
#         thl_in = profile[:,1]
#         qt_in  = profile[:,2]
#
#         # root mean square difference between vectors a and b
#         def rms(a, b):
#             return numpy.sqrt(((a-b)**2).mean())
#
#         for n in (1,2,4):
#             print "---------------------------------------------------- %d threads -------------------------------------------------"%n
#             instance = Dales(number_of_workers=n,**default_options)
#             tim=instance.get_model_time()
#             instance.commit_grid()
#
#             zf  = instance.get_zf()
#             thl = instance.get_profile_THL()
#             qt  = instance.get_profile_QT()
#             for i in range(0,instance.k):
#                 print "%3d %6.1f - %7.3f %7.3f - %8.4f %8.4f"%(i, zf[i].value_in(units.m), thl_in[i], thl[i].value_in(units.K), qt_in[i], qt[i].value_in(units.shu))
#             # get rid of units
#             thl = [t.value_in(units.K) for t in thl]
#             rmsT  = rms(thl, thl_in)
#             rmsQT = rms(qt, qt_in)
#             print n, ' threads:'
#             print 'RMS error in thl', rmsT
#             print 'RMS error in qt', rmsQT
#
#             # demand that the rms error btw input profiles and extracted average profiles is small
#             # DALES adds some random fluctuations, so the average != input
#             assert rmsT < 1e-2
#             assert rmsQT < 1e-6
#
#             #print "The retrieved THL profile is:", thl
#             #print "The retrieved QT profile is:", qt
#             #print "The retrieved zf levels:", instance.get_zf()
#             #print "The retrieved zh levels:", instance.get_zh()
#
#             i,j,k,xsize,ysize = instance.get_params_grid()
#             print "Grid size", (i, j, k)
#             print "Horizontal extent of model", (xsize,ysize)
#
#             instance.cleanup_code()
#             instance.stop()
#             print
#
#     def test8(self):
#         print "Test 8: instantiate, set a forcing, run 30 seconds and retrieve profiles"
#         # the effect of the forcing is seen in the U profile at the same height as the forcing was places
#         # no automatic test for this
#
#         instance = Dales(number_of_workers=4,**default_options)
#         tim=instance.get_model_time()
#         instance.commit_grid()
#
#         #print "Number of height levels:", instance.k
#         #tend = numpy.arange(instance.k)*.01
#         #print ('U_tend in python', tend)
#         #instance.set_tendency_U(tend)
#
#         tend = numpy.zeros(instance.k)
#         tend[5] = -.1
#         tend[6] = -.2   # layer 9 sees the most effect of this
#         tend[7] = -.1
#
#         instance.set_tendency_U(tend | units.m/units.s**2)
#
#         instance.evolve_model(tim + (30 | units.s))
#
#         U = instance.get_profile_U()
#         V = instance.get_profile_V()
#         W = instance.get_profile_W()
#         THL = instance.get_profile_THL()
#         QT = instance.get_profile_QT()
#
#         print "The retrieved U profile is:", U
#         print "The retrieved V profile is:", V
#         #print "The retrieved W profile is:", W
#         print
#
#         #print "The retrieved THL profile is:", THL
#         print "The retrieved QT profile is:", QT
#         #print "The retrieved zf levels:", instance.get_zf()
#         #print "The retrieved zh levels:", instance.get_zh()
#
#
#
#
#         i,j,k,xsize,ysize = instance.get_params_grid()
#         print "Grid size", (i, j, k)
#         print "Horizontal extent of model", (xsize,ysize)
#
#         instance.cleanup_code()
#         instance.stop()
#
#
#
#     def test9(self):
#         print "Test 9: get 3D data and compare with slab averages"
#
#         # root mean square difference between vectors a and b
#         def rms(a, b):
#             return (((a-b)**2).mean())**0.5
#
#         instance = Dales(number_of_workers=2,**default_options)
#         tim=instance.get_model_time()
#         instance.commit_grid()
#
#         tend = numpy.zeros(instance.k)
#         tend[5] = -1
#         instance.set_tendency_U(tend | units.m/units.s**2)
#
#         instance.evolve_model(tim + (180 | units.s))
#
#         def getU(imin=0, imax=instance.itot, jmin=0, jmax=instance.jtot, kmin=0, kmax=instance.k):
#             points = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
#             points = points.reshape(3, -1)
#             i,j,k = points
#             print 'i', i
#             print 'j', j
#             print 'k', k
#
#             field = instance.get_field_U(i,j,k)
#             field = field.reshape ((imax-imin, jmax-jmin, kmax-kmin))
#             return field
#
#
#         #field = field.reshape ((isize, jsize, ksize))
#
#         #print 'vertical profile U'
#         #print getU(imin=3,imax=4,jmin=8,jmax=9) # a vertical profile
#
#         #print 'horizontal profile U at k=5'
#         #print getU(imin=1,imax=10, jmin=8,jmax=9, kmin=5, kmax=6) # a vertical profile
#
#
#         #print 'horizontal profile U at k=6'
#         #print getU(imin=1,imax=10, jmin=8,jmax=9, kmin=6, kmax=7) # a vertical profile
#
#         # get slab averages
#         U = instance.get_profile_U()
#         V = instance.get_profile_V()
#         W = instance.get_profile_W()
#         THL = instance.get_profile_THL()
#         QT = instance.get_profile_QT()
#
#         # get 3D profiles
#         U3   = instance.get_field('U')
#         V3   = instance.get_field('V')
#         W3   = instance.get_field('W')
#         THL3 = instance.get_field('THL')
#         QT3  = instance.get_field('QT')
#
#
#         # calculate own slab averages of 3D profiles
#         U3a   = U3.mean(axis=0).mean(axis=0)
#         V3a   = V3.mean(axis=0).mean(axis=0)
#         W3a   = W3.mean(axis=0).mean(axis=0)
#         THL3a = THL3.mean(axis=0).mean(axis=0)
#         QT3a  = QT3.mean(axis=0).mean(axis=0)
#
#         # rms difference btw retrieved and calculated slab averages
#         rmsU    = rms(U,   U3a)
#         rmsV    = rms(V,   V3a)
#         rmsW    = rms(W,   W3a)
#         rmsTHL  = rms(THL, THL3a)
#         rmsQT   = rms(QT,  QT3a)
#
#         print 'RMS errors (U, V, W, THL, QT)', rmsU, rmsV, rmsW, rmsTHL, rmsQT
#         print THL3 [0, 0, 0:4]
# #        print instance.get_field('THL', imin)
#
#         assert rmsU   < 1e-12 | units.m / units.s
#         assert rmsV   < 1e-12 | units.m / units.s
#         assert rmsW   < 1e-12 | units.m / units.s
#         assert rmsTHL < 1e-12 | units.K
#         assert rmsQT  < 1e-12
#
#         instance.cleanup_code()
#         instance.stop()
#
#
#     def test10(self):
#         fileName = 'reference-data/test10.pkl'
#         save = False #True
#
#         print "Test 10: get 3D data, compare with stored reference"
#
#         # summed root mean square difference between vectors a and b
#         def rmsSum(a, b):
#             return (((a-b)**2).sum())**0.5
#
#         instance = Dales(number_of_workers=1,**default_options)
#         tim=instance.get_model_time()
#         instance.commit_grid()
#
#        # tend = numpy.zeros(instance.k)
#        # tend[5] = -1
#        # instance.set_tendency_U(tend)
#
#         instance.evolve_model(tim + (300 | units.s), exactEnd=0)
#
#         # get 3D profiles
#         U3   = instance.get_field('U')
#         V3   = instance.get_field('V')
#         W3   = instance.get_field('W')
#         THL3 = instance.get_field('THL')
#         QT3  = instance.get_field('QT')
#
#         if save:
#             outFile = open(fileName, 'wb')
#             pkl = (U3, V3, W3, THL3, QT3)
#             pickle.dump(pkl, outFile)
#             print('Saved reference profiles in', fileName)
#         else:
#             inFile = open(fileName, 'rb')
#             pkl = pickle.load(inFile)
#             U3p, V3p, W3p, THL3p, QT3p = pkl
#             rmsU   = rmsSum(U3, U3p)
#             rmsV   = rmsSum(V3, V3p)
#             rmsW   = rmsSum(W3, W3p)
#             rmsTHL = rmsSum(THL3, THL3p)
#             rmsQT  = rmsSum(QT3, QT3p)
#
#             print 'RMS errors (U, V, W, THL, QT)', rmsU, rmsV, rmsW, rmsTHL, rmsQT
#
#             assert rmsU   < 1e-9 | units.m / units.s
#             assert rmsV   < 1e-9 | units.m / units.s
#             assert rmsW   < 1e-9 | units.m / units.s
#             assert rmsTHL < 1e-9 | units.K
#             assert rmsQT  < 1e-9
#
#         instance.cleanup_code()
#         instance.stop()
#
#
#
#     def test11(self):
#         print("Test 11: set 3D data, get it back, compare")
#
#         fieldunits = {
#                 'U'   : units.m / units.s,
#                 'V'   : units.m / units.s,
#                 'W'   : units.m / units.s,
#                 'THL' : units.K,
#                 'QT'  : units.shu,
#             }
#
#         # write a single grid point at i,j,k, read back
#         def testSingle(i,j,k):
#             #a = numpy.random.random((1,1,1)) # 3D matrix with one element
#             a = numpy.random.random()  | fieldunits[field]# scalar
#             b = dales.get_field(field, i, i+1, j, j+1, k, k+1)
# #            b = drop_units(field, b)
#             dales.set_field(field, a, i, j, k)
#             c = dales.get_field(field, i, i+1, j, j+1, k, k+1)
# #            c = drop_units(field, c)
#             #~ print field, (i,j,k), 'was:', b, 'set: ', a, 'get:', c
#             self.assertEquals(a,c)
#
#         # write to a block [i:i+si, j:j+sj, k:k+sk], read back
#         def testBlock(i,j,k, si, sj, sk):
#             a = numpy.random.random((si,sj,sk)) | fieldunits[field]
#             b = dales.get_field(field, i, i+si, j, j+sj, k, k+sk)
# #            b = drop_units(field, b)
#             dales.set_field(field, a, i, j, k)
#             c = dales.get_field(field, i, i+si, j, j+sj, k, k+sk)
# #            c = drop_units(field, c)
#             print field, (i,j,k),  (si,sj,sk) #'was:', b, 'set: ', a, 'get:', c
#             self.assertEquals(a,c)
#
#         for n in (1,2,4):
#             print ' --- ', n, 'dales workers --- '
#             dales = Dales(number_of_workers=n,**default_options)
#             dales.commit_grid()
#
#             # set random grid elements for all fields, read back
#             for s in range(0,50):
#                 for field in ('U', 'V', 'W', 'THL', 'QT'):
#                     i = numpy.random.randint(0,dales.itot)
#                     j = numpy.random.randint(0,dales.jtot)
#                     k = numpy.random.randint(0,dales.k)
#                     testSingle(i, j, k)
#
#             for s in range(0,50):
#                 for field in ('U', 'V', 'W', 'THL', 'QT'):
#                     i = numpy.random.randint(0,dales.itot)
#                     j = numpy.random.randint(0,dales.jtot)
#                     k = numpy.random.randint(0,dales.k)
#
#                     si = numpy.random.randint(1,dales.itot-i+1)
#                     sj = numpy.random.randint(1,dales.jtot-j+1)
#                     sk = numpy.random.randint(1,dales.k-k+1)
#                     testBlock(i, j, k, si, sj, sk)
#
#
#
#
#             dales.cleanup_code()
#             dales.stop()
