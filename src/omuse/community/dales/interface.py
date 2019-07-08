import os.path
import shutil

import numpy
from amuse import datamodel
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community.interface.stopping_conditions import StoppingConditionInterface, StoppingConditions
from amuse.datamodel import new_rectilinear_grid
from amuse.rfi.core import CodeInterface
from amuse.rfi.core import remote_function
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.support.parameter_tools import CodeWithNamelistParameters
from .dalesreader import make_file_reader
from omuse.units import units
from parameters import namelist_parameters


# needs a bit of work for "non-local" runs
# (ie mechanism to copy input files...)

class DalesInterface(CodeInterface,
                     CommonCodeInterface,
                     StoppingConditionInterface,
                     LiteratureReferencesMixIn):
    """

    DALES - Dutch Atmospheric Large Eddy Simulation

    .. [#] Heus et al., Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications,  Geoscientific Model Development 3, 415, (2010)

    """
    use_modules = ["StoppingConditions", "dales_interface"]

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return "dales_worker"

    @remote_function
    def get_input_file():
        returns(input_file="s")

    @remote_function
    def set_input_file(input_file="namoptions"):
        pass

    @remote_function
    def set_qt_forcing(forcing_type=0):
        pass

    @remote_function
    def get_exact_end():
        returns(exactEndFlag=False)

    @remote_function
    def set_exact_end(exactEndFlag=False):
        pass

    @remote_function
    def change_dir(directory="."):
        pass

    @remote_function
    def set_start_date(date=0):
        pass

    @remote_function
    def set_start_time(time=0):
        pass

    @remote_function
    def get_model_time():
        returns(time=0. | units.s)

    @remote_function
    def get_timestep():
        returns(dt=0. | units.s)

    @remote_function
    def set_surface_pressure(p=1.0e5 | units.Pa):
        pass

    @remote_function
    def set_tendency_surface_pressure(p=1.0e5 | units.Pa / units.s):
        pass

    @remote_function
    def get_surface_pressure():
        returns(p=0. | units.Pa)

    @remote_function
    def commit_grid():
        pass

    # getter functions for vertical profiles - slab averages
    # these take a dummy array as input, and return output of the same length

    @remote_function(must_handle_array=True)
    def get_profile_U_(k=0):
        returns(out=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_V_(k=0):
        returns(out=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_W_(k=0):
        returns(out=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_THL_(k=0):
        returns(out=0. | units.K)

    @remote_function(must_handle_array=True)
    def get_profile_QT_(k=0):
        returns(out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_QL_(k=0):
        returns(out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_QL_ice_(k=0):
        returns(out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_QR_(k=0):
        returns(out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_E12_(k=0):
        returns(out=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_T_(k=0):
        returns(out=0. | units.K)

    # getter for cloud fraction. Uses the index array to define slabs.
    @remote_function(must_handle_array=True)
    def get_cloudfraction(k=0):
        returns(out=0. | units.m ** 2 / units.m ** 2)

    # getter for accumulated surface rain flux
#    @remote_function()
#    def get_rain():
#        returns(out=0. | units.kg / units.m ** 2)

    # getter functions for height levels
    # these take a dummy array as input, and return output of the same length
    @remote_function(must_handle_array=True)
    def get_zf_(k=0):
        returns(out=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_zh_(k=0):
        returns(out=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_presf_(k=0):
        returns(out=0. | units.Pa)

    @remote_function(must_handle_array=True)
    def get_presh_(k=0):
        returns(out=0. | units.Pa)

    @remote_function(must_handle_array=True)
    def get_rhof_(k=0):
        returns(out=0. | units.kg / units.m**3)

    @remote_function(must_handle_array=True)
    def get_rhobf_(k=0):
        returns(out=0. | units.kg / units.m**3)

    # setter functions for vertical tendencies / forcings
    @remote_function(must_handle_array=True)
    def set_tendency_U(a=0. | units.m / units.s ** 2):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_V(a=0. | units.m / units.s ** 2):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_THL(a=0. | units.K / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_QT(a=0. | units.mfu / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_QL(a=0. | units.mfu / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_ref_profile_QL(a=0. | units.mfu):
        returns()

    @remote_function(must_handle_array=True)
    def set_qt_variability_factor(a=0. | 1 / units.s):
        returns()

    # indexed getter/setter functions for vertical tendencies / forcings
    @remote_function(must_handle_array=True)
    def set_tendency_U_(g_i=0, a=0. | units.m / units.s ** 2):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_V_(g_i=0, a=0. | units.m / units.s ** 2):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_THL_(g_i=0, a=0. | units.K / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_tendency_QT_(g_i=0, a=0. | units.mfu / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def get_tendency_U_(g_i=0):
        returns(a=0. | units.m / units.s ** 2)

    @remote_function(must_handle_array=True)
    def get_tendency_V_(g_i=0):
        returns(a=0. | units.m / units.s ** 2)

    @remote_function(must_handle_array=True)
    def get_tendency_THL_(g_i=0):
        returns(a=0. | units.K / units.s)

    @remote_function(must_handle_array=True)
    def get_tendency_QT_(g_i=0):
        returns(a=0. | units.mfu / units.s)

    # indexed getter/setter functions for vertical nudging profiles
    @remote_function(must_handle_array=True)
    def set_nudge_U(g_i=0, a=0. | units.m / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_nudge_V(g_i=0, a=0. | units.m / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_nudge_THL(g_i=0, a=0. | units.K):
        returns()

    @remote_function(must_handle_array=True)
    def set_nudge_QT(g_i=0, a=0. | units.mfu):
        returns()

    @remote_function(must_handle_array=True)
    def get_nudge_U(g_i=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_nudge_V(g_i=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_nudge_THL(g_i=0):
        returns(a=0. | units.K)

    @remote_function(must_handle_array=True)
    def get_nudge_QT(g_i=0):
        returns(a=0. | units.mfu)

    @remote_function()
    def set_nudge_time_U(time=0. | units.s):
        returns()

    @remote_function()
    def set_nudge_time_V(time=0. | units.s):
        returns()

    @remote_function()
    def set_nudge_time_THL(time=0. | units.s):
        returns()

    @remote_function()
    def set_nudge_time_QT(time=0. | units.s):
        returns()

    @remote_function()
    def get_nudge_time_U():
        returns(time=0. | units.s)

    @remote_function()
    def get_nudge_time_V():
        returns(time=0. | units.s)

    @remote_function()
    def get_nudge_time_THL():
        returns(time=0. | units.s)

    @remote_function()
    def get_nudge_time_QT():
        returns(time=0. | units.s)

    # getter functions for 3D fields usning index arrays
    @remote_function(must_handle_array=True)
    def get_field_U(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_V(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_W(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_THL(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.K)

    @remote_function(must_handle_array=True)
    def get_field_QT(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_field_QL(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.mfu)

#    @remote_function(must_handle_array=True)
#    def get_field_Qsat(g_i=0, g_j=0, g_k=0):
#        returns(a=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_field_E12(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_T(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.K)

    @remote_function(must_handle_array=True)
    def get_field_rswd(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rswdir(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rswdif(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rswu(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rlwd(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rlwu(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rswdcs(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rswucs(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rlwdcs(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_rlwucs(g_i=0, g_j=0, g_k=0):
        returns(a=0. | units.W / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_LWP(g_i=0, g_j=0):
        returns(a=0. | units.kg / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_RWP(g_i=0, g_j=0):
        returns(a=0. | units.kg / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_TWP(g_i=0, g_j=0):
        returns(a=0. | units.kg / units.m ** 2)

    @remote_function(must_handle_array=True)
    def get_field_ustar(g_i=0, g_j=0):
        returns(a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_z0m(g_i=0, g_j=0):
        returns(a=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_field_z0h(g_i=0, g_j=0):
        returns(a=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_field_tskin(g_i=0, g_j=0):
        returns(a=0. | units.K)

    @remote_function(must_handle_array=True)
    def get_field_qskin(g_i=0, g_j=0):
        returns(a=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_field_LE(g_i=0, g_j=0):
        returns(a=0. | units.W / units.m**2)

    @remote_function(must_handle_array=True)
    def get_field_H(g_i=0, g_j=0):
        returns(a=0. | units.W / units.m**2)

    @remote_function(must_handle_array=True)
    def get_field_obl(g_i=0, g_j=0):
        returns(a=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_field_thlflux(g_i=0, g_j=0):
        returns(a=0. | units.K * units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_qtflux(g_i=0, g_j=0):
        returns(a=0. | units.mfu * units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_dudz(g_i=0, g_j=0):
        returns(a=0. | 1. / units.s)

    @remote_function(must_handle_array=True)
    def get_field_dvdz(g_i=0, g_j=0):
        returns(a=0. | 1. / units.s)

    @remote_function(must_handle_array=True)
    def get_field_dqtdz(g_i=0, g_j=0):
        returns(a=0. | units.mfu / units.m)

    @remote_function(must_handle_array=True)
    def get_field_dthldz(g_i=0, g_j=0):
        returns(a=0. | units.K / units.m)

    # setter functions for 3D fields usning index arrays
    @remote_function(must_handle_array=True)
    def set_field_U(g_i=0, g_j=0, g_k=0, a=0. | units.m / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_V(g_i=0, g_j=0, g_k=0, a=0. | units.m / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_W(g_i=0, g_j=0, g_k=0, a=0. | units.m / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_THL(g_i=0, g_j=0, g_k=0, a=0. | units.K):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_QT(g_i=0, g_j=0, g_k=0, a=0. | units.mfu):
        returns()

    #    @remote_function(must_handle_array=True)
    #    def set_field_E12(g_i=0,g_j=0,g_k=0,a=0.):
    #        returns()

    # getter functions for wtflux and qtflux
    @remote_function
    def get_wt_surf():
        returns(wtflux=0. | units.m * units.s ** -1 * units.K)

    @remote_function
    def get_wq_surf():
        returns(wqflux=0. | units.m / units.s)

    # setter functions for wtflux and qtflux
    @remote_function
    def set_wt_surf(wtflux=0. | units.m * units.s ** -1 * units.K):
        returns()

    @remote_function
    def set_wq_surf(wqflux=0. | units.m / units.s):
        returns()

    # setter functions for average momentum and heat roughness
    @remote_function
    def get_z0m_surf():
        returns(z0m=0. | units.m)

    @remote_function
    def get_z0h_surf():
        returns(z0h=0. | units.m)

    # setter functions for momentum and heat roughness
    @remote_function
    def set_z0m_surf(z0=0. | units.m):
        returns()

    @remote_function
    def set_z0h_surf(z0=0. | units.m):
        returns()

    @remote_function()
    def get_params_grid():
        returns(i=1, j=1, k=1, xsize=1.0 | units.m, ysize=1.0 | units.m)

    # exactEnd: if true, step exactly to tend,
    # otherwise, step past tend with normal time steps (default)
    # the end time form namoptions is effective only if exactEnd is false
    @remote_function
    def evolve_model(tend=0. | units.s, exactEnd=0):
        returns(walltime=0. | units.s)

    @remote_function
    def write_restart():
        returns()

    @remote_function
    def write_restart():
        returns()

class Dales(CommonCode, CodeWithNamelistParameters):
    QT_FORCING_GLOBAL = 0
    QT_FORCING_LOCAL = 1
    QT_FORCING_VARIANCE = 2
    QT_FORCING_STRONG = 3

    input_file_options = {"namelist": "namoptions", "profiles": "prof.inp", "forcings": "lscale.inp",
                          "scalars": "scalar.inp"}

    def __init__(self, **options):
        # Namelist file
        CodeWithNamelistParameters.__init__(self, namelist_parameters)
        CommonCode.__init__(self, DalesInterface(**options), **options)
        self.stopping_conditions = StoppingConditions(self)

        inputdir = None

        # Set working directory
        self.workdir = "."
        if "workdir" in options:
            self.workdir = options["workdir"]
            if os.path.exists(self.workdir): # This is a restart
                inputdir = self.workdir

        # Set experiment name
        if "exp" in options:
            self.parameters_RUN.iexpnr = int(options.get("exp"))

        # Set input directory
        if "inputdir" in options:
            if inputdir is None:
                inputdir = options["inputdir"]
            else:
                print "Input directory specification %s ignored, because it was already set to %s" % \
                      (options["inputdir"], inputdir)

        # Look up input directory
        if "case" in options:
            if inputdir is None:
                inputdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dales-repo", "cases",
                                        options["case"])
            else:
                print "Dales case directory specification %s ignored, because it was already set to %s" % \
                      (options["case"], inputdir)

        # Register input data files
        input_files = {k: None for k in Dales.input_file_options.keys()}
        if inputdir is not None:
            for opt in Dales.input_file_options:
                default_name = Dales.input_file_options[opt]
                input_files[opt] = os.path.join(inputdir,
                                                options.get(opt, default_name + ".%3.3i" % self.parameters_RUN.iexpnr))
        else:
            for opt in Dales.input_file_options:
                input_files[opt] = options.get(opt, None)

        self.extra_input_files = []
        if inputdir is not None:
            for f in os.listdir(inputdir):
                filepath = os.path.join(inputdir, f)
                if os.path.isfile(filepath) and f.endswith(".%3.3i" % self.parameters_RUN.iexpnr) and filepath not in \
                        input_files.values():
                    self.extra_input_files.append(filepath)

        namelist = input_files.get("namelist", None)
        self.parameters.input_file = namelist
        if namelist is not None and os.path.isfile(namelist):
            self.read_namelist_parameters(namelist)
            self._nml_file = namelist
        else:
            print "Could not find file %s, falling back to default options" % namelist
            self._nml_file = None

        self.initial_profile_grid = self.read_initial_profiles(input_files.get("profiles", None))
        self.initial_large_scale_forcings_grid = self.read_initial_forcings(input_files.get("forcings", None),
                                                                            self.initial_profile_grid.z)
        self.initial_scalar_grid = self.read_initial_scalars(input_files.get("scalars", None),
                                                             self.parameters_RUN.nsv,
                                                             self.initial_profile_grid.z)

    @staticmethod
    def read_initial_profiles(filepath):
        if filepath is None:
            zf = numpy.arange(0., 5000., 50.) | units.m
            grid = new_rectilinear_grid((len(zf),), cell_centers=[zf], axes_names=["z"])
            grid.qt = numpy.full(zf.shape, 0.005) | units.shu
            grid.thl = numpy.full(zf.shape, 292.5) | units.K
            grid.u = numpy.full(zf.shape, 1. / numpy.sqrt(2.)) | units.m / units.s
            grid.v = numpy.full(zf.shape, 1. / numpy.sqrt(2.)) | units.m / units.s
            grid.e12 = numpy.full(zf.shape, 1.) | units.m / units.s
        else:
            reader = make_file_reader(filepath)
            var = reader.variables[0] if len(reader.variables) > 0 else "zf"
            zf = reader.get_heights(var) | units.m
            grid = new_rectilinear_grid((len(zf),), cell_centers=[zf], axes_names=["z"])

            def get_vars(varname, value=0.):
                return reader[varname][:] if varname in reader.variables else numpy.full(zf.shape, value)

            grid.qt = get_vars("qt") | units.shu
            grid.thl = get_vars("thl", 288.) | units.K
            grid.u = get_vars("u") | units.m / units.s
            grid.v = get_vars("v") | units.m / units.s
            grid.e12 = get_vars("tke", 1.) | units.m / units.s
        return grid

    @staticmethod
    def read_initial_forcings(filepath, zf):
        if filepath is None or not os.path.isfile(filepath):
            grid = new_rectilinear_grid((len(zf),), cell_centers=[zf], axes_names=["z"])
            grid.ug = numpy.full(zf.shape, 1. / numpy.sqrt(2.)) | units.m / units.s
            grid.vg = numpy.full(zf.shape, 1. / numpy.sqrt(2.)) | units.m / units.s
            grid.wfls = numpy.full(zf.shape, 0.) | units.m / units.s
            grid.dqtdxls = numpy.full(zf.shape, 0.) | units.shu / units.s
            grid.dqtdyls = numpy.full(zf.shape, 0.) | units.shu / units.s
            grid.dqtdtls = numpy.full(zf.shape, 0.) | units.shu * units.m / units.s ** 2
            grid.dthlrad = numpy.full(zf.shape, 0.) | units.shu / units.s
        else:
            reader = make_file_reader(filepath)
            grid = new_rectilinear_grid((len(zf),), cell_centers=[zf], axes_names=["z"])

            def get_vars(varname, z_in, z_out):
                if varname not in reader.variables:
                    return numpy.zeros(z_out.shape)
                vals = reader[varname][:]
                if numpy.array_equal(z_in, z_out):
                    return vals
                else:
                    return numpy.interp(z_out, z_in, vals)

            var = reader.variables[0] if len(reader.variables) > 0 else "zf"
            z = reader.get_heights(var)
            grid.ug = get_vars("ug", z, zf.value_in(units.m)) | units.m / units.s
            grid.vg = get_vars("vg", z, zf.value_in(units.m)) | units.m / units.s
            grid.wfls = get_vars("wfls", z, zf.value_in(units.m)) | units.m / units.s
            grid.dqtdxls = get_vars("dqtdxls", z, zf.value_in(units.m)) | units.shu / units.s
            grid.dqtdyls = get_vars("dqtdyls", z, zf.value_in(units.m)) | units.shu / units.s
            grid.dqtdtls = get_vars("dqtdtls", z, zf.value_in(units.m)) | units.shu * units.m / units.s ** 2
            grid.dthlrad = get_vars("dthlrad", z, zf.value_in(units.m)) | units.shu / units.s
        return grid

    @staticmethod
    def read_initial_scalars(filepath, num_scalars, zf):
        if filepath is None or not os.path.isfile(filepath):
            grid = new_rectilinear_grid((len(zf),), cell_centers=[zf], axes_names=["z"])
            for i in range(1, num_scalars + 1):
                setattr(grid, "sv" + str(i), numpy.zeros(zf.shape) | units.shu)
        else:
            reader = make_file_reader(filepath)
            grid = new_rectilinear_grid((len(zf),), cell_centers=[zf], axes_names=["z"])
            for v in reader.variables:
                if v != "height":
                    setattr(grid, v, reader[v][:] | units.shu)
        return grid

    @staticmethod
    def write_initial_profile_file(grid, filepath):
        header = """ profile
    zf         thl        qt         u          v          tke             """
        numpy.savetxt(filepath, numpy.column_stack((
            grid.z.value_in(units.m),
            grid.thl.value_in(units.K),
            grid.qt.value_in(units.shu),
            grid.u.value_in(units.m / units.s),
            grid.v.value_in(units.m / units.s),
            grid.e12.value_in(units.m / units.s))), header=header)

    @staticmethod
    def write_large_scale_forcing_file(grid, filepath):
        header = """ large scale forcing
    height     ug         vg         wfls       dqtdxls    dqtqyls    dqtdtls    dthlrad"""
        numpy.savetxt(filepath, numpy.column_stack((
            grid.z.value_in(units.m),
            grid.ug.value_in(units.m / units.s),
            grid.vg.value_in(units.m / units.s),
            grid.wfls.value_in(units.m / units.s),
            grid.dqtdxls.value_in(units.shu / units.s),
            grid.dqtdyls.value_in(units.shu / units.s),
            grid.dqtdtls.value_in(units.shu * units.m / units.s ** 2),
            grid.dthlrad.value_in(units.shu / units.s),
        )), header=header)

    @staticmethod
    def write_initial_scalars(grid, filepath):
        header = """ input scalar profiles
    height"""
        cols = [grid.z.value_in(units.m)]
        for att in grid.get_attribute_names_defined_in_store():
            if att.startswith("sv"):
                header += "             " + att
                cols.append(grid.get_all_values_of_attribute_in_store(att).value_in(units.shu))
        numpy.savetxt(filepath, numpy.column_stack(tuple(cols)), header=header)

    def commit_parameters(self):
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)
        iexpnr = int(self.parameters_RUN.iexpnr)
        if self._nml_file:
            filename = os.path.basename(self._nml_file).split('.')[0]
        else:
            filename = "namoptions"
        filename += ".%3.3i" % iexpnr
        self.change_dir(os.path.abspath(self.workdir))
        self.set_input_file(filename)
        filepath = os.path.join(self.workdir, filename)
        replace = False
        if filepath == self._nml_file:
            filepath += '~'
            replace = True
        CodeWithNamelistParameters.write_namelist_parameters(self, filepath, do_patch=self._nml_file)
        if replace:
            os.rename(filepath, filepath[:-1])
        Dales.write_initial_profile_file(self.initial_profile_grid,
                                         filepath=os.path.join(self.workdir, "prof.inp.%3.3i" % iexpnr))
        Dales.write_initial_scalars(self.initial_scalar_grid,
                                    filepath=os.path.join(self.workdir, "scalar.inp.%3.3i" % iexpnr))
        Dales.write_large_scale_forcing_file(self.initial_large_scale_forcings_grid,
                                             filepath=os.path.join(self.workdir, "lscale.inp.%3.3i" % iexpnr))
        self.workdir = os.path.abspath(self.workdir)
        for f in self.extra_input_files:
            if os.path.abspath(f) != os.path.join(self.workdir, os.path.basename(f)):
                shutil.copy2(f, self.workdir)

        # print "code options written to %s"%dalesinputfile

        # check this...
        dt = self.parameters.starttime
        if dt != 0:
            self.set_start_date(10000 * dt.year + 100 * dt.month + dt.day)
            self.set_start_time(10000 * dt.hour + 100 * dt.minute + dt.second)
        self.set_qt_forcing(self.parameters.qt_forcing)

        self.overridden().commit_parameters()

    def define_parameters(self, obj):
        CodeWithNamelistParameters.define_parameters(self, obj)

        obj.add_interface_parameter(
            "input_file",
            "the input file name",
            None,
            "before_set_interface_parameter"
        )
        obj.add_alias_parameter(
            "restart_flag",
            "lwarmstart",
            "warm start from restart file if True",
            alias_set="parameters_RUN"
        )
        obj.add_alias_parameter(
            "restart_file",
            "startfile",
            "restart file name",
            alias_set="parameters_RUN"
        )
        obj.add_alias_parameter(
            "trestart",
            "trestart",
            "(simulation) time between writing restart files",
            alias_set="parameters_RUN"
        )
        obj.add_method_parameter(
            "get_itot",
            None,
            "itot",
            "number of cells in the x direction",
            None
        )
        obj.add_method_parameter(
            "get_jtot",
            None,
            "jtot",
            "number of cells in the y direction",
            None
        )
        obj.add_method_parameter(
            "get_xsize",
            None,
            "xsize",
            "size of the grid the x direction",
            None
        )
        obj.add_method_parameter(
            "get_ysize",
            None,
            "ysize",
            "size of the grid the y direction",
            None
        )
        obj.add_method_parameter(
            "get_dx",
            None,
            "dx",
            "grid spacing in x direction",
            None
        )
        obj.add_method_parameter(
            "get_dy",
            None,
            "dy",
            "grid spacing in y direction",
            None
        )
        obj.add_boolean_parameter(
            "get_exact_end",
            "set_exact_end",
            "exactEndFlag",
            "parameter specifying whether evolve should end exactly at target time ",
            False
        )
        obj.add_interface_parameter(
            "starttime",
            "absolute model start datetime",
            0,
            "before_set_interface_parameter"
        )
        obj.add_interface_parameter(
            "qt_forcing",
            "qt forcing type",
            0,
            "before_set_interface_parameter"
        )

    def define_properties(self, obj):
        obj.add_property('get_model_time', public_name="model_time")
        obj.add_property('get_timestep', public_name="timestep")

    def define_state(self, obj):
        obj.set_initial_state("UNINITIALIZED")
        obj.add_transition("UNINITIALIZED", "INITIALIZED", "initialize_code")
        obj.add_method("!UNINITIALIZED", "before_get_parameter")
        obj.add_transition("!UNINITIALIZED!STOPPED!INITIALIZED", "END", "cleanup_code")
        obj.add_transition("END", "STOPPED", "stop", False)
        obj.add_method("STOPPED", 'stop')
        obj.add_method("UNINITIALIZED", 'stop')
        obj.add_method("INITIALIZED", 'stop')
        obj.add_transition("INITIALIZED", "EDIT", "commit_parameters")
        obj.add_transition("EDIT", "RUN", "commit_grid")
        obj.add_method("INITIALIZED", "set_workdir")
        obj.add_method("INITIALIZED", "set_input_file")
        obj.add_method("INITIALIZED", "set_start_date_time")

        for state in ["EDIT", "RUN", "EVOLVED"]:
            obj.add_method(state, "get_model_time")
            obj.add_method(state, "get_timestep")

        obj.add_transition("RUN", "EVOLVED", "evolve_model", False)
        obj.add_method("EVOLVED", "evolve_model")

    # wrapping functions for hiding the dummy array passed to getter functions
    def get_profile_U(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_U_(indices)

    def get_profile_V(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_V_(indices)

    def get_profile_W(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_W_(indices)

    def get_profile_THL(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_THL_(indices)

    def get_profile_QT(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_QT_(indices)

    def get_profile_QL(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_QL_(indices)

    def get_profile_QL_ice(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_QL_ice_(indices)

    def get_profile_QR(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_QR_(indices)

    def get_profile_E12(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_E12_(indices)

    def get_profile_T(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_profile_T_(indices)

    def get_zf(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_zf_(indices)

    def get_zh(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_zh_(indices)

    def get_presh(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_presh_(indices)

    def get_presf(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_presf_(indices)
    
    def get_rhof(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_rhof_(indices)

    def get_rhobf(self, k=None):
        if k is None:
            kmin, kmax = self.get_z_grid_range()
            indices = numpy.arange(kmin, kmax + 1)
        else:
            indices = k
        return self.get_rhobf_(indices)

    # retrieve a 3D field
    # indices are one-based, for zero-based index access use the grids
    # field is 'U', 'V', 'W', 'THL', 'QT', 'QL', 'E12', 'T'
    def get_field(self, field, imin=1, imax=None, jmin=1, jmax=None, kmin=1, kmax=None):

        grid_range = ()
        if imax is None or jmax is None or kmax is None:
            grid_range = self.get_grid_range()
        if imax is None:
            imax = grid_range[1] + 1
        if jmax is None:
            jmax = grid_range[3] + 1
        if kmax is None:
            kmax = grid_range[5] + 1

        # build index arrays
        if field in ('LWP', 'RWP', 'TWP'):  # 2D field
            points = numpy.mgrid[imin:imax, jmin:jmax]
            points = points.reshape(2, -1)
            i, j = points
            k = 0
        else:  # 3D field
            points = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
            points = points.reshape(3, -1)
            i, j, k = points

        if field == 'U':
            field = self.get_field_U(i, j, k)
        elif field == 'V':
            field = self.get_field_V(i, j, k)
        elif field == 'W':
            field = self.get_field_W(i, j, k)
        elif field == 'THL':
            field = self.get_field_THL(i, j, k)
        elif field == 'QT':
            field = self.get_field_QT(i, j, k)
        elif field == 'QL':
            field = self.get_field_QL(i, j, k)
        elif field == 'Qsat':
            field = self.get_field_Qsat(i, j, k)
        elif field == 'E12':
            field = self.get_field_E12(i, j, k)
        elif field == 'T':
            field = self.get_field_T(i, j, k)
        elif field == 'LWP':  # LWP - 2D field, k ignored
            field = self.get_field_LWP(i, j)
            return field.reshape((imax - imin, jmax - jmin))  # separate return here, to reshape to 2D
        elif field == 'RWP':  # RWP - 2D field, k ignored
            field = self.get_field_RWP(i, j)
            return field.reshape((imax - imin, jmax - jmin))  # separate return here, to reshape to 2D
        elif field == 'TWP':  # TWP - 2D field, k ignored
            field = self.get_field_TWP(i, j)
            return field.reshape((imax - imin, jmax - jmin))  # separate return here, to reshape to 2D
        else:
            raise Exception('get_field called with undefined variable name %s' % field)

        return field.reshape((imax - imin, jmax - jmin, kmax - kmin))

    # set a 3D field
    # indices are one-based, for zero-based index access use the grids
    # field is 'U', 'V', 'W', 'THL', 'QT'
    def set_field(self, field, a, imin=1, jmin=1, kmin=1):

        # set max indices from the size of a 
        try:
            imax = imin + a.shape[0]
            jmax = jmin + a.shape[1]
            kmax = kmin + a.shape[2]
            i, j, k = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
            i = i.flatten()
            j = j.flatten()
            k = k.flatten()
            a = a.flatten()
        except:
            i = imin
            j = jmin
            k = kmin

        if field == 'U':
            self.set_field_U(i, j, k, a)
        elif field == 'V':
            self.set_field_V(i, j, k, a)
        elif field == 'W':
            self.set_field_W(i, j, k, a)
        elif field == 'THL':
            self.set_field_THL(i, j, k, a)
        elif field == 'QT':
            self.set_field_QT(i, j, k, a)
        else:
            raise Exception('set_field called with undefined variable name %s' % field)

    # retrieve a 1D vertical profile - wrapper function consistent with get_field
    # field is 'U', 'V', 'W', 'THL', 'QT', 'QL', 'E12', 'T'
    def get_profile(self, field):
        profile = None
        kmin, kmax = self.get_z_grid_range()
        indices = numpy.arange(kmin, kmax + 1)
        print kmin, kmax
        if field == 'U':
            profile = self.get_profile_U_(indices)
        elif field == 'V':
            profile = self.get_profile_V_(indices)
        elif field == 'W':
            profile = self.get_profile_W_(indices)
        elif field == 'THL':
            profile = self.get_profile_THL_(indices)
        elif field == 'QT':
            profile = self.get_profile_QT_(indices)
        elif field == 'QL':
            profile = self.get_profile_QL_(indices)
        elif field == 'E12':
            profile = self.get_profile_E12_(indices)
        elif field == 'T':
            profile = self.get_profile_T_(indices)
        else:
            print('get_profile called with undefined field', field)
        return profile

    def get_itot(self):
        return self.get_params_grid()[0]

    def get_jtot(self):
        return self.get_params_grid()[1]

    def get_ktot(self):
        return self.get_params_grid()[2]

    def get_xsize(self):
        return self.get_params_grid()[3]

    def get_ysize(self):
        return self.get_params_grid()[4]

    def get_dx(self):
        itot, jtot, ktot, x, y = self.get_params_grid()
        return x / itot

    def get_dy(self):
        itot, jtot, ktot, x, y = self.get_params_grid()
        return y / jtot

    def get_grid_range(self):
        itot, jtot, ktot, x, y = self.get_params_grid()
        return 1, itot, 1, jtot, 1, ktot

    def get_grid_position(self, i, j, k):
        itot, jtot, ktot, x, y = self.get_params_grid()
        return i * x / itot, j * y / jtot, self.get_zf(k)

    def get_xy_grid_range(self):
        imin, imax, jmin, jmax, kmin, kmax = self.get_grid_range()
        return imin, imax, jmin, jmax

    def get_xy_grid_position(self, i, j):
        x, y, z = self.get_grid_position(i, j, 0)
        return x, y

    def get_z_grid_range(self):
        imin, imax, jmin, jmax, kmin, kmax = self.get_grid_range()
        return kmin, kmax

    def get_z_grid_position(self, k):
        x, y, z = self.get_grid_position(0, 0, k)
        return z

    def get_scalar_grid_range(self):
        return ()

    def define_grids(self, obj):
        obj.define_grid('grid', axes_names="xyz", grid_class=datamodel.RectilinearGrid,
                        state_guard="before_new_set_instance")
        obj.set_grid_range('grid', 'get_grid_range')
        obj.add_getter('grid', 'get_grid_position', names="xyz")
        for x in ['U', 'V', 'W', 'THL', 'QT', 'QL', 'E12', 'T', 'rswd', 'rswdir', 'rswdif', 'rswu', 'rlwd', 'rlwu',
                  'rswdcs', 'rswucs', 'rlwdcs', 'rlwucs']:
            obj.add_getter('grid', 'get_field_' + x, names=[x])
        for x in ['U', 'V', 'W', 'THL', 'QT']:
            obj.add_setter('grid', 'set_field_' + x, names=[x])

        obj.define_grid('profiles', axes_names="z", grid_class=datamodel.RectilinearGrid,
                        state_guard="before_new_set_instance")
        obj.set_grid_range('profiles', 'get_z_grid_range')
        obj.add_getter('profiles', 'get_z_grid_position', names="z")
        for x in ['U', 'V', 'W', 'THL', 'QT', 'QL', 'QL_ice', 'E12', 'T']:
            obj.add_getter('profiles', 'get_profile_' + x + '_', names=[x])

        # nudge grid  -experimental-
        obj.define_grid('nudging', axes_names="z", grid_class=datamodel.RectilinearGrid,
                        state_guard="before_new_set_instance")
        obj.set_grid_range('nudging', 'get_z_grid_range')
        obj.add_getter('nudging', 'get_z_grid_position', names="z")
        for x in ['U', 'V', 'THL', 'QT']:
            obj.add_getter('nudging', 'get_nudge_' + x, names=[x])
            obj.add_setter('nudging', 'set_nudge_' + x, names=[x])

        obj.define_grid('forcings', axes_names="z", grid_class=datamodel.RectilinearGrid,
                        state_guard="before_new_set_instance")
        obj.set_grid_range('forcings', 'get_z_grid_range')
        obj.add_getter('forcings', 'get_z_grid_position', names="z")
        for x in ['U', 'V', 'THL', 'QT']:
            obj.add_getter('forcings', 'get_tendency_' + x + '_', names=['tendency_' + x])
            obj.add_setter('forcings', 'set_tendency_' + x + '_', names=['tendency_' + x])

        obj.define_grid('scalars', grid_class=datamodel.RectilinearGrid, state_guard="before_new_set_instance")
        obj.set_grid_range('scalars', 'get_scalar_grid_range')
        obj.add_getter('scalars', 'get_surface_pressure', names=['pressure'])
        # obj.add_getter('scalars', 'get_rain', names=['rain'])
        for name in ["wt", "wq", "z0m", "z0h"]:
            obj.add_getter('scalars', 'get_' + name + "_surf", names=[name])
            obj.add_setter('scalars', 'set_' + name + "_surf", names=[name])

        obj.define_grid('surface_fields', grid_class=datamodel.RectilinearGrid, state_guard="before_new_set_instance")
        obj.set_grid_range('surface_fields', 'get_xy_grid_range')
        obj.add_getter('surface_fields', 'get_xy_grid_position', names="xy")
        for name in ["LWP", "RWP", "TWP", "ustar", "z0m", "z0h", "tskin", "qskin", "LE", "H", "obl", "thlflux",
                     "qtflux", "dudz", "dvdz", "dqtdz", "dthldz"]:
            obj.add_getter('surface_fields', 'get_field_' + name, names=[name])
