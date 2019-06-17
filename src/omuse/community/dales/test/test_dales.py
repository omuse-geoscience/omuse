import logging

from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from omuse.community.dales.interface import Dales
from omuse.units import units

kwargs = {}
# kwargs=dict(channel_type="sockets",redirection="none")
# kwargs=dict(redirection="none",debugger="gdb")

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)


def cleanup_data(rundir):
    if os.path.isdir(rundir):
        for f in os.listdir(rundir):
            os.remove(os.path.join(rundir, f))
        os.rmdir(rundir)
    elif os.path.isfile(rundir):
        os.remove(rundir)


# Unsupported test cases: aerosolrad, chem, example, fog, heterogen, hireslapse, neutral

class TestDalesInterface(TestWithMPI):

    def test_namopt_file_is_written(self):
        rundir = "work0"
        os.mkdir(rundir)
        os.chdir(rundir)
        try:
            instance = Dales(**kwargs)
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
            instance = Dales(**kwargs)
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
            instance = Dales(**kwargs)
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
        instance = Dales(workdir=rundir, **kwargs)
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_namopt_file_written_in_workdir(self):
        rundir = "work4"
        instance = Dales(workdir=rundir, **kwargs)
        instance.commit_parameters()
        assert os.path.isfile(os.path.join(rundir, "namoptions.001"))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_load_rico(self):
        rundir = "work-rico"
        instance = Dales(case="rico", workdir=rundir, **kwargs)
        instance.commit_parameters()
        assert instance.get_ktot() == 126
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_rico(self):
        rundir = "work-rico2"
        instance = Dales(case="rico", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_cblstrong(self):
        rundir = "work-cblstrong"
        instance = Dales(case="cblstrong", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (30 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_cblweak(self):
        rundir = "work-cblweak"
        instance = Dales(case="cblweak", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (30 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_load_bomex(self):
        rundir = "work-bomex"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_bomex(self):
        rundir = "work-bomex2"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (1 | units.minute))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_load_sp_case(self):
        rundir = "work-sp"
        instance = Dales(case="sp-testcase", workdir=rundir, **kwargs)
        instance.commit_parameters()
        assert os.path.exists(rundir)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_get_grid_dimensions(self):
        rundir = "work-sp2"
        instance = Dales(case="sp-testcase", workdir=rundir, **kwargs)
        instance.commit_parameters()
        assert (instance.get_itot(), instance.get_jtot(), instance.get_ktot()) == (200, 200, 160)
        assert (instance.get_dx().value_in(units.m), instance.get_dy().value_in(units.m)) == (200, 200)
        assert numpy.array_equal(instance.get_zf().value_in(units.m), numpy.arange(12.5, 4000., 25.))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_sp_case(self):
        rundir = "work-sp3"
        instance = Dales(case="sp-testcase", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (1 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_upscale_arm_brown_case(self):
        rundir = "work-arm-brown"
        instance = Dales(case="arm_brown", workdir=rundir, **kwargs)
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        instance.commit_parameters()
        assert os.path.isfile(os.path.join(rundir, "ls_flux.inp.001"))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_grid_shapes(self):
        rundir = "work-sp4"
        instance = Dales(case="sp-testcase", workdir=rundir, **kwargs)
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        tim = instance.get_model_time()
        instance.evolve_model(tim + (1 | units.s))
        temp_profile = instance.profiles.T
        assert temp_profile.shape == (instance.get_ktot())
        v_block = instance.grid.V
        assert v_block.shape == (instance.get_itot(), instance.get_jtot(), instance.get_ktot())
        u_block = instance.get_field("THL")
        assert u_block.shape == (instance.get_itot(), instance.get_jtot(), instance.get_ktot())
        twp_slab = instance.surface_fields.LWP
        assert twp_slab.shape == (instance.get_itot(), instance.get_jtot())
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_u_profile_bomex(self):
        rundir = "work-bomex3"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (5 | units.s))
        u_profile1 = instance.profiles.U.value_in(units.m / units.s)
        u_profile2 = instance.get_profile_U().value_in(units.m / units.s)
        u_profile3 = instance.get_profile('U').value_in(units.m / units.s)
        print "Lengths:", len(u_profile1), len(u_profile2), len(u_profile3)
        assert numpy.allclose(u_profile1, u_profile2, rtol=1.e-16)
        assert numpy.allclose(u_profile1, u_profile3, rtol=1.e-16)
        u_field = instance.grid.U.value_in(units.m / units.s)
        assert numpy.allclose(u_profile1, numpy.mean(u_field, axis=(0, 1)), rtol=1.e-9)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_qt_profile_sp(self):
        rundir = "work-sp5"
        instance = Dales(case="sp-testcase", workdir=rundir, number_of_workers=4, **kwargs)
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
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
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
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (70 | units.s))
        thl2 = numpy.mean(instance.grid.THL.value_in(units.K)[:, :, 1], axis=(0, 1))
        assert thl2 < thl1 < 304.
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_bomex_3d_grid_interface(self):
        rundir = "work-bomex45"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (2 | units.s))
        qblock = numpy.full((3, 3, 1), 1.234e-5) | units.mfu
        instance.set_field("QT", qblock, imin=8, jmin=4, kmin=10)
        qt = instance.get_field("QT", imin=8, imax=11, jmin=4, jmax=7, kmin=10, kmax=11)
        qtgrid = instance.grid[7:10, 3:6, 9:10].QT.value_in(units.mfu)
        assert numpy.array_equal(qblock, qt.value_in(units.mfu))
        assert numpy.array_equal(qtgrid, qt.value_in(units.mfu))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_set_bomex_t_value(self):
        rundir = "work-bomex5"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
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
        assert thl1 == thl2 == 302.
        assert thl3 == thl4 == 298.
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_set_bomex_q_forcing(self):
        rundir = "work-bomex6"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
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
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_set_bomex_q_nudging(self):
        rundir = "work-bomex7"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
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
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_atex_scalar_fields(self):
        rundir = "work-atex"
        instance = Dales(case="atex", workdir=rundir, **kwargs)
        instance.parameters_DOMAIN.itot = 32
        instance.parameters_DOMAIN.jtot = 32
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        instance.set_wt_surf(0.1 | units.m * units.s ** -1 * units.K)
        wt = instance.get_wt_surf().value_in(units.m * units.s ** -1 * units.K)
        wtgrid = instance.scalars.wt.value_in(units.m * units.s ** -1 * units.K)
        assert (wt == wtgrid == 0.1)
        instance.set_z0h_surf(0.2 | units.m)
        instance.evolve_model(tim + (1 | units.s))
        z0 = instance.get_z0h_surf().value_in(units.m)
        z0grid = instance.scalars.z0h[0].value_in(units.m)
        assert (z0 == z0grid == 0.2)
        z0h = instance.surface_fields.z0h.value_in(units.m)
        assert (numpy.all(z0h == 0.2))
        instance.set_z0m_surf(0.3 | units.m)
        z0 = instance.get_z0m_surf().value_in(units.m)
        z0grid = instance.scalars.z0m[0].value_in(units.m)
        assert (z0 == z0grid == 0.3)
        z0h = instance.surface_fields.z0m.value_in(units.m)
        assert (numpy.all(z0h == 0.3))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_bomex_2d_fields(self):
        rundir = "work-bomex8"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (60 | units.s))
        vars2d = {"TWP": units.kg / units.m ** 2,
                  "ustar": units.m / units.s,
                  "tskin": units.K,
                  "qskin": units.mfu,
                  "LE": units.W / units.m ** 2,
                  "H": units.W / units.m ** 2,
                  "obl": units.m}
        for variable, unit in vars2d.items():
            field = getattr(instance.surface_fields, variable).value_in(unit)
            assert (field.shape == (instance.get_itot(), instance.get_jtot()))
        twp = instance.surface_fields.TWP.value_in(vars2d["TWP"])
        qt = instance.grid.QT.value_in(units.mfu)
        dz = instance.get_zf().value_in(units.m)[-1] / instance.get_ktot()
        twp_check = (numpy.sum(qt, axis=2) * dz)
        assert (numpy.allclose(twp, twp_check, rtol=0.1))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_rico_restart(self):
        rundir = "work-rico3"
        instance1 = Dales(case="rico", workdir=rundir, **kwargs)
        tim = instance1.get_model_time()
        instance1.evolve_model(tim + (20 | units.s))
        instance1.write_restart()
        instance1.cleanup_code()
        instance1.stop()
        instance2 = Dales(workdir=rundir, **kwargs)
        instance2.evolve_model(tim + (40 | units.s))
        t1 = instance2.grid.T.value_in(units.K)
        instance2.cleanup_code()
        instance2.stop()
        cleanup_data(rundir)
        rundir = "work-rico4"
        instance3 = Dales(case="rico", workdir=rundir, **kwargs)
        tim = instance3.get_model_time()
        instance3.evolve_model(tim + (60 | units.s))
        t2 = instance3.grid.T.value_in(units.K)
        instance3.cleanup_code()
        instance3.stop()
        cleanup_data(rundir)
        assert (numpy.allclose(t1, t2, rtol=1.e-3))

    def test_bomex_nudge(self):
        rundir = "work-bomex9"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        instance.commit_parameters()
        zf = instance.nudging.z.value_in(units.m)
        qt1 = instance.profiles.QT.value_in(units.shu)
        zprof = 1. / ((zf / zf[-1] - 0.66) ** 2 + 1)
        instance.nudging.QT = numpy.sum(qt1) * zprof / numpy.sum(zprof) | units.shu
        instance.set_nudge_time_QT(100. | units.s)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (100 | units.s))
        qt2 = instance.profiles.QT.value_in(units.shu)
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)
        rundir = "work-bomex10"
        instance = Dales(case="bomex", workdir=rundir, **kwargs)
        instance.evolve_model(tim + (100 | units.s))
        qt3 = instance.profiles.QT.value_in(units.shu)
        assert (not numpy.array_equal(qt1, qt2) and not numpy.array_equal(qt2, qt3))
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_arm_brown(self):
        rundir = "arm_brown"
        instance = Dales(case="arm_brown", workdir=rundir, **kwargs)
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_arm_unstable(self):
        rundir = "arm_unstable"
        instance = Dales(case="arm_unstable", workdir=rundir, **kwargs)
        instance.parameters_DOMAIN.itot = 64
        instance.parameters_DOMAIN.jtot = 64
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_dycoms(self):
        rundir = "dycoms"
        instance = Dales(case="dycoms_rf02", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_gabls(self):
        rundir = "gabls1"
        instance = Dales(case="gabls1", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

    def test_run_smoke(self):
        rundir = "smoke"
        instance = Dales(case="smoke", workdir=rundir, **kwargs)
        tim = instance.get_model_time()
        instance.evolve_model(tim + (10 | units.s))
        newtim = instance.get_model_time()
        assert newtim > tim
        instance.cleanup_code()
        instance.stop()
        cleanup_data(rundir)

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
