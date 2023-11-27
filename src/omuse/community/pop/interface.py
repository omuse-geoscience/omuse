import os

import numpy
from amuse.community.interface.common import CommonCode, CommonCodeInterface
from amuse.community.interface.stopping_conditions import (
    StoppingConditionInterface,
    StoppingConditions,
)
from amuse.datamodel import StructuredGrid
from amuse.datamodel.staggeredgrid import StaggeredGrid
from amuse.rfi.core import CodeInterface, legacy_function, remote_function
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.support.parameter_tools import CodeWithNamelistParameters
from omuse.units import units


def compute_cell_corners(nodes=None, u_lon=None, u_lat=None):
    if u_lon is None:
        u_lon = nodes.lon.value_in(units.rad)
    if u_lat is None:
        u_lat = nodes.lat.value_in(units.rad)

    size = u_lon.shape

    corners = numpy.zeros((2, size[0] + 1, size[1] + 1), dtype=numpy.double)

    corners[1, 1:, 1:] = u_lat  # u_lat is north east corner of t-cell
    corners[1, 0, 1:] = u_lat[-1, :]  # copy last column to first column of corners
    tiny = 1.0e-14
    corners[1, :, 0] = (-numpy.pi * 0.5) + tiny  # mock-up latitudes of row 0

    corners[0, 1:, 1:] = u_lon  # u_lon is north east corner of t-cell
    corners[0, 0, 1:] = u_lon[-1, :]  # copy last column to first column of corners
    corners[0, :, 0] = corners[0, :, 1]  # copy row 1 of corners down to mocked-up row 0

    return corners


class POPInterface(CodeInterface, LiteratureReferencesMixIn):

    """
    POP - Parallel Ocean Program

    .. [#] http://www.cesm.ucar.edu/models/cesm1.0/pop2/doc/sci/POPRefManual.pdf
    .. [#] https://github.com/NLeSC/eSalsa-POP

    """

    include_headers = ["worker_code.h"]
    use_modules = ["pop_interface"]

    MODE_NORMAL = "320x384x40"
    MODE_HIGH = "3600x2400x42"
    MODE_320x384x40 = "320x384x40"
    MODE_3600x2400x42 = "3600x2400x42"
    OTHER_MODES = ["96x40x12", "96x120x12", "120x56x12", "240x110x12", "test"]

    def __init__(self, mode=MODE_NORMAL, **keyword_arguments):
        self.mode = mode

        keyword_arguments.setdefault("number_of_workers", 8)
        CodeInterface.__init__(
            self, name_of_the_worker=self.name_of_the_worker(mode), **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)

        # ~ if 'AMUSE_DIR' in os.environ:
        # ~ amuse_dir = os.environ['AMUSE_DIR']
        # ~ path = amuse_dir + '/src/omuse/community/pop/'
        # ~ self.change_directory(path)

        if mode in self.OTHER_MODES:
            pass
        elif mode in [self.MODE_NORMAL, self.MODE_320x384x40]:
            pass
        elif mode in [self.MODE_HIGH, self.MODE_3600x2400x42]:
            pass
        else:
            raise Exception("Unknown mode")

    def name_of_the_worker(self, mode):
        return "pop_worker_" + mode

    ##forcings getters and setters
    @remote_function(must_handle_array=True)
    def get_node_wind_stress(i=0, j=0):
        returns(tau_x=0.0 | units.Pa, tau_y=0.0 | units.Pa)

    @remote_function(must_handle_array=True)
    def set_node_wind_stress(i=0, j=0, tau_x=0.0 | units.Pa, tau_y=0.0 | units.Pa):
        returns()

    @remote_function(must_handle_array=True)
    def get_node_coriolis_f(i=0, j=0):
        returns(coriolis_f=0.0 | units.s**-1)

    @remote_function(must_handle_array=True)
    def set_node_coriolis_f(i=0, j=0, coriolis_f=0.0 | units.s**-1):
        returns()

    ##getters for node and element state
    @remote_function(must_handle_array=True)
    def get_node_surface_state(i=0, j=0):
        returns(
            vx=0.0 | units.cm / units.s,
            vy=0.0 | units.cm / units.s,
            gradx=0.0 | units.s**-1,
            grady=0.0 | units.s**-1,
        )

    @remote_function(must_handle_array=True)
    def set_node_surface_state(
        i=0, j=0, gradx=0.0 | units.s**-1, grady=0.0 | units.s**-1
    ):
        returns()

    @remote_function(must_handle_array=True)
    def get_node_barotropic_vel_curtime(i=0, j=0):
        returns(
            vx_barotropic=0.0 | units.cm / units.s,
            vy_barotropic=0.0 | units.cm / units.s,
        )

    @remote_function(must_handle_array=True)
    def get_node_barotropic_vel_oldtime(i=0, j=0):
        returns(
            vx_barotropic=0.0 | units.cm / units.s,
            vy_barotropic=0.0 | units.cm / units.s,
        )


    @remote_function(must_handle_array=True)
    def set_node_barotropic_vel_alltime(
        i=0,
        j=0,
        vx_barotropic=0.0 | units.cm / units.s,
        vy_barotropic=0.0 | units.cm / units.s,
    ):
        returns()

    @remote_function(must_handle_array=True)
    def set_node_barotropic_vel_curtime(
        i=0,
        j=0,
        vx_barotropic=0.0 | units.cm / units.s,
        vy_barotropic=0.0 | units.cm / units.s,
    ):
        returns()

    @remote_function(must_handle_array=True)
    def set_node_barotropic_vel_oldtime(
        i=0,
        j=0,
        vx_barotropic=0.0 | units.cm / units.s,
        vy_barotropic=0.0 | units.cm / units.s,
    ):
        returns()

    @remote_function(must_handle_array=True)
    def get_element_surface_state(i=0, j=0):
        returns(
            temp=0.0 | units.Celsius, salt=0.0 | units.g / units.g
        )  # salinity in gram per kg

    @remote_function(must_handle_array=True)
    def get_element_surface_heat_flux(i=0, j=0):
        returns(
            surface_heat_flux=0.0 | units.W / units.m**2
        )  # surface heat flux in W/m**2

    @remote_function(must_handle_array=True)
    def get_element_surface_forcing_temp(i=0, j=0):
        returns(restoring_temp=0.0 | units.Celsius)

    @remote_function(must_handle_array=True)
    def set_element_surface_forcing_temp(i=0, j=0, restoring_temp=0.0 | units.Celsius):
        returns()

    @remote_function(must_handle_array=True)
    def get_element_surface_forcing_salt(i=0, j=0):
        returns(restoring_salt=0.0 | units.g / units.g)  # to be checked

    @remote_function(must_handle_array=True)
    def set_element_surface_forcing_salt(
        i=0, j=0, restoring_salt=0.0 | units.g / units.g
    ):
        returns()

    @remote_function(must_handle_array=True)
    def get_element_ssh_curtime(i=0, j=0):
        returns(ssh=0.0 | units.cm)

    @remote_function(must_handle_array=True)
    def get_element_ssh_oldtime(i=0, j=0):
        returns(ssh=0.0 | units.cm)

    @remote_function(must_handle_array=True)
    def get_element_ssh_guess(i=0, j=0):
        returns(ssh=0.0 | units.cm)

    @remote_function(must_handle_array=True)
    def set_element_ssh_alltime(i=0, j=0, ssh=0.0 | units.cm):
        returns()

    @remote_function(must_handle_array=True)
    def set_element_ssh_curtime(i=0, j=0, ssh=0.0 | units.cm):
        returns()

    @remote_function(must_handle_array=True)
    def set_element_ssh_oldtime(i=0, j=0, ssh=0.0 | units.cm):
        returns()

    @remote_function(must_handle_array=True)
    def set_element_ssh_guess(i=0, j=0, ssh=0.0 | units.cm):
        returns()

    ##these two return in units.m in adcirc but here they return latitude longitude in degrees
    @remote_function(must_handle_array=True)
    def get_node_position(i=0, j=0):
        returns(lat=0.0 | units.rad, lon=0.0 | units.rad)

    @remote_function(must_handle_array=True)
    def get_element_position(i=0, j=0):
        returns(lat=0.0 | units.rad, lon=0.0 | units.rad)

    # vertical position, may vary horizontally when using partial bottom cells
    @remote_function(must_handle_array=True)
    def get_node_vposition(i=0, j=0, k=0):
        returns(depth=0.0 | units.cm)

    @remote_function(must_handle_array=True)
    def get_element_vposition(i=0, j=0, k=0):
        returns(depth=0.0 | units.cm)

    # distance from the surface to mid point of layer in vertical grid
    @remote_function(must_handle_array=True)
    def get_zt(k=0):
        returns(z=0.0 | units.cm)

    @remote_function
    def get_number_of_nodes():
        returns(n_nodes=0)

    # get number of gridpoints in x-direction (west-east) and y-direction (south-north)
    @remote_function
    def get_domain_size():
        returns(nx=0, ny=0)

    # get the total number of vertical levels used
    @remote_function
    def get_number_of_vertical_levels():
        returns(km=0)

    # get&set the vertical grid spacings
    @remote_function(can_handle_array=True)
    def get_dz(k=0):
        returns(dz=0.0 | units.cm)

    @remote_function(can_handle_array=True)
    def set_dz(k=0, dz=0.0 | units.cm):
        returns()

    # returns the maximal depth at this position
    @remote_function(must_handle_array=True)
    def get_node_depth(i=0, j=0):
        returns(depth=0.0 | units.cm)

    @remote_function(must_handle_array=True)
    def get_element_depth(i=0, j=0):
        returns(depth=0.0 | units.cm)

    @remote_function(must_handle_array=True)
    def set_KMT(i=0, j=0, depth_level=0):
        returns()

    @remote_function(must_handle_array=True)
    def get_KMT(i=0, j=0):
        returns(depth_level=0)

    # returns x velocity for 3D grid
    @remote_function(must_handle_array=True)
    def get_node3d_velocity_xvel_curtime(i=0, j=0, k=0):
        returns(xvel=0.0 | units.cm / units.s)

    @remote_function(must_handle_array=True)
    def get_node3d_velocity_xvel_oldtime(i=0, j=0, k=0):
        returns(xvel=0.0 | units.cm / units.s)

    @remote_function(must_handle_array=True)
    def set_node3d_velocity_xvel_alltime(i=0, j=0, k=0, xvel=0.0 | units.cm / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_node3d_velocity_xvel_curtime(i=0, j=0, k=0, xvel=0.0 | units.cm / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_node3d_velocity_xvel_oldtime(i=0, j=0, k=0, xvel=0.0 | units.cm / units.s):
        returns()

    # returns y velocity for 3D grid
    @remote_function(must_handle_array=True)
    def get_node3d_velocity_yvel_curtime(i=0, j=0, k=0):
        returns(yvel=0.0 | units.cm / units.s)

    @remote_function(must_handle_array=True)
    def get_node3d_velocity_yvel_oldtime(i=0, j=0, k=0):
        returns(yvel=0.0 | units.cm / units.s)

    @remote_function(must_handle_array=True)
    def set_node3d_velocity_yvel_alltime(i=0, j=0, k=0, yvel=0.0 | units.cm / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_node3d_velocity_yvel_curtime(i=0, j=0, k=0, yvel=0.0 | units.cm / units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_node3d_velocity_yvel_oldtime(i=0, j=0, k=0, yvel=0.0 | units.cm / units.s):
        returns()

    # returns z velocity for 3D grid
    @remote_function(must_handle_array=True)
    def get_node3d_velocity_zvel(i=0, j=0, k=0):
        returns(zvel=0.0 | units.cm / units.s)

    # returns temperature for 3D grid
    @remote_function(must_handle_array=True)
    def get_element3d_temperature(i=0, j=0, k=0):
        returns(temp=0.0 | units.Celsius)

    @remote_function(must_handle_array=True)
    def set_element3d_temperature(i=0, j=0, k=0, temp=0.0 | units.Celsius):
        returns()

    # returns salinity for 3D grid
    @remote_function(must_handle_array=True)
    def get_element3d_salinity(i=0, j=0, k=0):
        returns(salt=0.0 | units.g / units.g)

    @remote_function(must_handle_array=True)
    def set_element3d_salinity(i=0, j=0, k=0, salt=0.0 | units.g / units.g):
        returns()

    # returns densiity for 3D grid
    @remote_function(must_handle_array=True)
    def get_element3d_density(i=0, j=0, k=0):
        returns(rho=0.0 | units.g / units.cm**3)

    @remote_function(must_handle_array=True)
    def set_element3d_density(i=0, j=0, k=0, rho=0.0 | units.g / units.cm**3):
        returns()

    @remote_function
    def initialize_code():
        returns()

    @remote_function
    def evolve_model(tend=0.0 | units.s):
        pass

    @remote_function
    def cleanup_code():
        returns()

    @remote_function
    def get_model_time():
        returns(time=0.0 | units.s)

    @remote_function
    def get_timestep():
        returns(dt=0.0 | units.s)

    @remote_function
    def get_timestep_next():
        returns(dt=0.0 | units.s)

    # facilitate state transitions
    @remote_function
    def prepare_parameters():
        returns()

    @remote_function
    def commit_parameters():
        returns()

    @remote_function
    def recommit_forcings_():
        returns()

    def prepare_forcings(self):
        pass

    @remote_function
    def recommit_grids():
        returns()

    def recommit_forcings(self):
        self.recommit_forcings_()
        self.prepare_parameters()  # to update the windstress to the one just set

    @remote_function
    def get_horiz_grid_option():
        returns(option="s")

    @remote_function
    def set_horiz_grid_option(option="s"):
        returns()

    @remote_function
    def get_horiz_grid_file():
        returns(filename="s")

    @remote_function
    def set_horiz_grid_file(filename="s"):
        returns()

    @remote_function
    def get_vert_grid_option():
        returns(option="s")

    @remote_function
    def set_vert_grid_option(option="s"):
        returns()

    @remote_function
    def get_vert_grid_file():
        returns(filename="s")

    @remote_function
    def set_vert_grid_file(filename="s"):
        returns()

    @remote_function
    def get_topography_option():
        returns(option="s")

    @remote_function
    def set_topography_option(option="s"):
        returns()

    @remote_function
    def get_topography_file():
        returns(filename="s")

    @remote_function
    def set_topography_file(filename="s"):
        returns()

    @remote_function
    def get_bottom_cell_file():
        returns(filename="s")

    @remote_function
    def set_bottom_cell_file(filename="s"):
        returns()

    @remote_function
    def get_ts_option():
        returns(option="s")

    @remote_function
    def set_ts_option(option="s"):
        returns()

    @remote_function
    def get_ts_file():
        returns(filename="s")

    @remote_function
    def set_ts_file(filename="s"):
        returns()

    @remote_function
    def get_ts_file_format():
        returns(filename="s")

    @remote_function
    def set_ts_file_format(filename="s"):
        returns()

    @remote_function
    def set_nprocs(nprocs=0):
        returns()

    @remote_function
    def get_distribution():
        returns(option="s")

    @remote_function
    def set_distribution(option="s"):
        returns()

    @remote_function
    def get_distribution_file():
        returns(filename="s")

    @remote_function
    def set_distribution_file(filename="s"):
        returns()

    @remote_function
    def get_ew_boundary_type():
        returns(option="s")

    @remote_function
    def set_ew_boundary_type(option="s"):
        returns()

    @remote_function
    def get_ns_boundary_type():
        returns(option="s")

    @remote_function
    def set_ns_boundary_type(option="s"):
        returns()

    @remote_function
    def get_lonmin():
        returns(lonmin=0.0 | units.deg)

    @remote_function
    def set_lonmin(lonmin=0.0 | units.deg):
        returns()

    @remote_function
    def get_lonmax():
        returns(lonmax=0.0 | units.deg)

    @remote_function
    def set_lonmax(lonmax=0.0 | units.deg):
        returns()

    @remote_function
    def get_latmin():
        returns(latmin=0.0 | units.deg)

    @remote_function
    def set_latmin(latmin=0.0 | units.deg):
        returns()

    @remote_function
    def get_latmax():
        returns(latmax=0.0 | units.deg)

    @remote_function
    def set_latmax(latmax=0.0 | units.deg):
        returns()

    @remote_function
    def get_reinit_gradp():
        returns(reinit_gradp=False)

    @remote_function
    def set_reinit_gradp(reinit_gradp=False):
        returns()

    @remote_function
    def get_reinit_rho():
        returns(reinit_rho=False)

    @remote_function
    def set_reinit_rho(reinit_rho=False):
        returns()

    @remote_function
    def get_restart_option():
        returns(option="s")

    @remote_function
    def set_restart_option(option="s"):
        returns()

    @remote_function
    def get_restart_freq_option():
        returns(option=0)

    @remote_function
    def set_restart_freq_option(option=0):
        returns()

    @remote_function
    def get_restart_file():
        returns(filename="s")

    @remote_function
    def set_restart_file(filename="s"):
        returns()

    @remote_function
    def get_tavg_option():
        returns(option="s")

    @remote_function
    def set_tavg_option(option="s"):
        returns()

    @remote_function
    def get_tavg_freq_option():
        returns(option=0)

    @remote_function
    def set_tavg_freq_option(option=0):
        returns()

    @remote_function
    def get_tavg_file():
        returns(filename="s")

    @remote_function
    def set_tavg_file(filename="s"):
        returns()

    @remote_function
    def get_movie_option():
        returns(option="s")

    @remote_function
    def set_movie_option(option="s"):
        returns()

    @remote_function
    def get_movie_freq_option():
        returns(option=0)

    @remote_function
    def set_movie_freq_option(option=0):
        returns()

    @remote_function
    def get_movie_file():
        returns(filename="s")

    @remote_function
    def set_movie_file(filename="s"):
        returns()

    @remote_function
    def get_runid():
        returns(option="s")

    @remote_function
    def set_runid(option="s"):
        returns()

    @remote_function
    def get_dt_option():
        returns(filename="s")

    @remote_function
    def set_dt_option(filename="s"):
        returns()

    @remote_function
    def get_dt_count():
        returns(option=0)

    @remote_function
    def set_dt_count(option=0):
        returns()

    @remote_function
    def change_directory(pathname="s"):
        returns()

    # forcing options
    @remote_function
    def set_shf_data_type(option="s"):
        returns()

    @remote_function
    def get_shf_filename():
        returns(filename="s")

    @remote_function
    def get_shf_data_type():
        returns(type="s")

    @remote_function
    def set_shf_monthly_file(filename="s"):
        returns()

    @remote_function
    def set_sfwf_data_type(option="s"):
        returns()

    @remote_function
    def get_sfwf_filename():
        returns(filename="s")

    @remote_function
    def get_sfwf_data_type():
        returns(type="s")

    @remote_function
    def set_sfwf_monthly_file(filename="s"):
        returns()

    @remote_function
    def get_ws_filename():
        returns(filename="s")

    @remote_function
    def get_ws_data_type():
        returns(type_="s")

    @remote_function
    def set_ws_data_type(type_="s"):
        returns()

    @remote_function
    def set_ws_monthly_file(filename="s"):
        returns()

    @remote_function
    def get_fwf_imposed():
        returns(fwf=0.0 | units.Sv)

    @remote_function
    def set_fwf_imposed(fwf=0.0 | units.Sv):
        returns()

    @remote_function
    def get_namelist_filename():
        returns(filename="s")

    @remote_function
    def set_namelist_filename(filename="s"):
        returns()


class POP(CommonCode, CodeWithNamelistParameters):
    nprocs = 0

    def __init__(
        self, mode=POPInterface.MODE_NORMAL, namelist_file="pop_in", **options
    ):
        self._nml_file = namelist_file
        self.nprocs = options.setdefault("number_of_workers", 8)
        CodeWithNamelistParameters.__init__(self, options.get("parameters", dict()))
        CommonCode.__init__(self, POPInterface(mode=mode, **options), **options)
        self.parameters.namelist_file = self._nml_file
        # ~ self.parameters.namelist_filename = "amuse.nml"

    def configuration_file_set(self):
        self.read_namelist_parameters(
            self.parameters.namelist_file, add_missing_parameters=True
        )
        handler = self.get_handler("PARAMETER")
        CodeWithNamelistParameters.define_parameters(self, handler)

    def define_properties(self, object):
        object.add_property("get_model_time", public_name="model_time")
        object.add_property("get_timestep", public_name="timestep")
        object.add_property("get_timestep_next", public_name="timestep_next")

    def define_state(self, object):
        object.set_initial_state("UNINITIALIZED")
        object.add_method("UNINITIALIZED", "change_directory")
        object.add_method("UNINITIALIZED", "set_namelist_filename")

        object.add_transition(
            "!UNINITIALIZED!INITIALIZED!STOPPED", "END", "cleanup_code"
        )
        object.add_transition("UNINITIALIZED", "INITIALIZED", "initialize_code")
        object.add_transition("END", "STOPPED", "stop", False)
        object.add_method("STOPPED", "stop")

        # some parameters have to be derived before they can be read, therefore
        # reading anything should be preceded by a call to prepare_parameters
        # in the two-phase init, prepare_parameters is called by commit_parameters()
        # commit_parameters must be called to finish the initialization
        object.add_transition("INITIALIZED", "EDIT", "commit_parameters")
        object.add_transition("INITIALIZED", "RUN", "commit_parameters")

        object.add_method("INITIALIZED", "set_horiz_grid_option")
        object.add_method("INITIALIZED", "set_horiz_grid_file")
        object.add_method("INITIALIZED", "set_vert_grid_option")
        object.add_method("INITIALIZED", "set_vert_grid_file")
        object.add_method("INITIALIZED", "set_topography_option")
        object.add_method("INITIALIZED", "set_topography_file")

        object.add_method("INITIALIZED", "set_shf_data_type")
        object.add_method("INITIALIZED", "set_shf_monthly_file")
        object.add_method("INITIALIZED", "set_sfwf_data_type")
        object.add_method("INITIALIZED", "set_sfwf_monthly_file")
        object.add_method("INITIALIZED", "set_ws_monthly_file")

        # ~ object.add_method('INITIALIZED', 'get_ts_option')
        object.add_method("INITIALIZED", "set_ts_option")
        object.add_method("INITIALIZED", "set_dz")
        # ~ object.add_method('INITIALIZED', 'get_ts_file')
        object.add_method("INITIALIZED", "set_ts_file")
        # ~ object.add_method('INITIALIZED', 'get_ts_file_format')
        object.add_method("INITIALIZED", "set_ts_file_format")
        object.add_method("INITIALIZED", "set_nprocs")
        # ~ object.add_method('INITIALIZED', 'get_distribution')
        object.add_method("INITIALIZED", "set_distribution")
        # ~ object.add_method('INITIALIZED', 'get_distribution_file')
        object.add_method("INITIALIZED", "set_distribution_file")
        # ~ object.add_method('INITIALIZED', 'get_ew_boundary_type')
        object.add_method("INITIALIZED", "set_ew_boundary_type")
        # ~ object.add_method('INITIALIZED', 'get_ns_boundary_type')
        object.add_method("INITIALIZED", "set_ns_boundary_type")

        # ~ object.add_method('INITIALIZED', 'get_restart_option')
        object.add_method("INITIALIZED", "set_restart_option")
        # ~ object.add_method('INITIALIZED', 'get_restart_freq_option')
        object.add_method("INITIALIZED", "set_restart_freq_option")
        # ~ object.add_method('INITIALIZED', 'get_restart_file')
        object.add_method("INITIALIZED", "set_restart_file")
        # ~ object.add_method('INITIALIZED', 'get_tavg_option')
        object.add_method("INITIALIZED", "set_tavg_option")
        # ~ object.add_method('INITIALIZED', 'get_tavg_freq_option')
        object.add_method("INITIALIZED", "set_tavg_freq_option")
        # ~ object.add_method('INITIALIZED', 'get_tavg_file')
        object.add_method("INITIALIZED", "set_tavg_file")
        # ~ object.add_method('INITIALIZED', 'get_movie_option')
        object.add_method("INITIALIZED", "set_movie_option")
        # ~ object.add_method('INITIALIZED', 'get_movie_freq_option')
        object.add_method("INITIALIZED", "set_movie_freq_option")
        # ~ object.add_method('INITIALIZED', 'get_movie_file')
        object.add_method("INITIALIZED", "set_movie_file")

        # ~ object.add_method('INITIALIZED', 'get_runid')
        object.add_method("INITIALIZED", "set_runid")
        # ~ object.add_method('INITIALIZED', 'get_dt_option')
        object.add_method("INITIALIZED", "set_dt_option")
        # ~ object.add_method('INITIALIZED', 'get_dt_count')
        object.add_method("INITIALIZED", "set_dt_count")

        object.add_method("INITIALIZED", "set_KMT")

        # you can only edit stuff in state EDIT
        # ~ object.add_method('INITIALIZED','before_set_parameter')
        # you can only read stuff in states RUN and EDIT
        # ~ object.add_method('EDIT','before_get_parameter')
        # ~ object.add_method('RUN','before_get_parameter')
        # ~ object.add_method('INITIALIZED','before_get_parameter')
        # setting a parameter is the only way to endup in state EDIT from RUN
        # ~ object.add_transition('RUN','EDIT','before_set_parameter')
        # ~ object.add_transition('RUN','EDIT_FORCINGS','before_set_parameter')

        for state in ["RUN", "EDIT"]:
            object.add_method(state, "get_dz")
            object.add_method(state, "get_node_coriolis_f")
            object.add_method(state, "get_node_wind_stress")
            object.add_method(state, "get_node_position")
            object.add_method(state, "get_element_position")
            object.add_method(state, "get_node_position")
            object.add_method(state, "get_element_surface_state")
            object.add_method(state, "get_element_surface_heat_flux")
            object.add_method(state, "get_element_surface_forcing_temp")
            object.add_method(state, "get_element_surface_forcing_salt")
            object.add_method(state, "get_node_surface_state")
            object.add_method(state, "get_node_barotropic_vel_curtime")
            object.add_method(state, "get_node_barotropic_vel_oldtime")
            object.add_method(state, "get_node3d_velocity_xvel_curtime")
            object.add_method(state, "get_node3d_velocity_xvel_oldtime")
            object.add_method(state, "get_node3d_velocity_yvel_curtime")
            object.add_method(state, "get_node3d_velocity_yvel_oldtime")
            object.add_method(state, "get_element3d_temperature")
            object.add_method(state, "get_element3d_salinity")
            object.add_method(state, "get_element3d_density")
            object.add_method(state, "get_element_ssh")
            object.add_method(state, "get_element_ssh_curtime")
            object.add_method(state, "get_element_ssh_oldtime")
            object.add_method(state, "get_element_ssh_guess")

        object.add_method(
            "RUN", "get_node3d_velocity_zvel"
        )  # because depends on xvel, yvel; needs prepare first

        object.add_method("EDIT", "set_node_barotropic_vel_alltime")
        object.add_method("EDIT", "set_node_barotropic_vel_curtime")
        object.add_method("EDIT", "set_node_barotropic_vel_oldtime")
        object.add_method("EDIT", "set_node_surface_state")
        object.add_method("EDIT", "set_node3d_velocity_xvel_alltime")
        object.add_method("EDIT", "set_node3d_velocity_xvel_curtime")
        object.add_method("EDIT", "set_node3d_velocity_xvel_oldtime")
        object.add_method("EDIT", "set_node3d_velocity_yvel_alltime")
        object.add_method("EDIT", "set_node3d_velocity_yvel_curtime")
        object.add_method("EDIT", "set_node3d_velocity_yvel_oldtime")
        object.add_method("EDIT", "set_element3d_temperature")
        object.add_method("EDIT", "set_element3d_salinity")
        object.add_method("EDIT", "set_element3d_density")
        object.add_method("EDIT", "set_element_ssh")
        object.add_method("EDIT", "set_element_ssh_curtime")
        object.add_method("EDIT", "set_element_ssh_oldtime")
        object.add_method("EDIT", "set_element_ssh_guess")

        object.add_method("EDIT_FORCINGS", "set_node_wind_stress")
        object.add_method("EDIT_FORCINGS", "set_node_coriolis_f")
        object.add_method("EDIT_FORCINGS", "set_element_surface_forcing_temp")
        object.add_method("EDIT_FORCINGS", "set_element_surface_forcing_salt")

        # before we can run the model we need to recommit_parameters
        object.add_transition("EDIT", "RUN", "recommit_grids")
        object.add_transition("EDIT_FORCINGS", "RUN", "recommit_forcings")
        object.add_transition("EDIT_FORCINGS", "EDIT", "recommit_forcings")
        object.add_transition(
            "EDIT", "EDIT_FORCINGS", "prepare_forcings"
        )  # prepare_forcings is a dummy
        object.add_transition("RUN", "EDIT_FORCINGS", "prepare_forcings")
        object.add_transition("RUN", "EDIT", "prepare_forcings")

        # can only evolve from run
        object.add_transition("RUN", "EVOLVED", "evolve_model", False)

        # from state EVOLVED you may evolve again if you like
        object.add_method("EVOLVED", "evolve_model")

        # prepare model for reading/writing after evolve
        object.add_transition("EVOLVED", "RUN", "prepare_parameters")
        object.add_transition("EVOLVED", "EDIT", "prepare_parameters")

        # ensure that get_model_time only returns valid values
        object.add_method("RUN", "get_model_time")
        object.add_method("EDIT", "get_model_time")
        object.add_method("EDIT_FORCINGS", "get_model_time")
        object.add_method("EVOLVED", "get_model_time")

    def initialize_code(self):
        self.write_namelist_parameters("amuse.nml", do_patch=True)
        self.overridden().initialize_code()

    def commit_parameters(self):
        self.set_nprocs(self.nprocs)

        try:
            if self.parameters_HMIX_ANISO_NML.var_viscosity_outfile:
                dirname = os.path.dirname(
                    self.parameters_HMIX_ANISO_NML.var_viscosity_outfile
                )
                os.makedirs(dirname)
        except:
            pass

        if self.parameters.vert_grid_option == "amuse":
            kmax = self.get_number_of_vertical_levels()
            if len(self.parameters.vertical_layer_thicknesses) == kmax:
                self.set_dz(
                    range(1, kmax + 1), self.parameters.vertical_layer_thicknesses
                )
            else:
                raise Exception(
                    "length of parameter vertical_layer_thicknesses needs to be {0}".format(
                        kmax
                    )
                )
        if self.parameters.topography_option == "amuse":
            d = numpy.zeros(self.get_domain_size())
            try:
                d = self.parameters.depth_index
                indices = numpy.indices(d.shape)
                self.set_KMT(indices[0].flatten() + 1, indices[1].flatten() + 1, d.flat)
            except:
                raise Exception("invalid depth index")
        self.overridden().commit_parameters()

    def get_firstlast_node(self):
        size = self.get_domain_size()
        return 1, size[0], 1, size[1]

    def get_firstlast_grid3d(self):
        size = self.get_domain_size()
        km = self.get_number_of_vertical_levels()
        return 1, size[0], 1, size[1], 1, km

    def get_ugrid_latlon_range(self):
        start = self.get_node_position(1, 1)
        size = self.get_domain_size()
        end = self.get_node_position(size[0], size[1])
        return start[0], end[0], start[1], end[1]

    def get_tgrid_latlon_range(self):
        start = self.get_element_position(1, 1)
        size = self.get_domain_size()
        end = self.get_element_position(size[0], size[1])
        return start[0], end[0], start[1], end[1]

    def get_element3d_position(self, i, j, k):
        lat, lon = self.get_element_position(i, j)
        z = self.get_element_vposition(i, j, k)
        return lat, lon, z

    def get_node3d_position(self, i, j, k):
        lat, lon = self.get_node_position(i, j)
        z = self.get_node_vposition(i, j, k)
        return lat, lon, z

    def _compute_cell_corners(self, u_lon=None, u_lat=None):
        return compute_cell_corners(self.nodes, u_lon, u_lat)

    def get_grid(self):
        return StaggeredGrid(self.elements, self.nodes, self._compute_cell_corners)

    def define_particle_sets(self, object):
        # for now we refer to the U grid as the nodes and T grid as the elements
        # the forcings can be on either grid, depends on what is preferred for coupling with adcirc I guess
        axes_names = ["lon", "lat"]

        object.define_grid("nodes", axes_names=axes_names, grid_class=StructuredGrid)
        object.set_grid_range("nodes", "get_firstlast_node")
        object.add_getter("nodes", "get_node_position", names=("lat", "lon"))
        object.add_getter("nodes", "get_node_depth", names=("depth",))
        object.add_getter(
            "nodes", "get_node_surface_state", names=("vx", "vy", "gradx", "grady")
        )
        object.add_getter(
            "nodes", "get_node_barotropic_vel_curtime", names=("vx_barotropic_all", "vy_barotropic_all")
        )
        object.add_getter(
            "nodes", "get_node_barotropic_vel_curtime", names=("vx_barotropic_cur", "vy_barotropic_cur")
        )
        object.add_getter(
            "nodes", "get_node_barotropic_vel_oldtime", names=("vx_barotropic_old", "vy_barotropic_old")
        )
        object.add_setter("nodes", "set_node_surface_state", names=("gradx", "grady"))
        object.add_setter(
            "nodes", "set_node_barotropic_vel_alltime", names=("vx_barotropic_all", "vy_barotropic_all")
        )
        object.add_setter(
            "nodes", "set_node_barotropic_vel_curtime", names=("vx_barotropic_cur", "vy_barotropic_cur")
        )
        object.add_setter(
            "nodes", "set_node_barotropic_vel_oldtime", names=("vx_barotropic_old", "vy_barotropic_old")
        )

        object.define_grid("nodes3d")
        object.set_grid_range("nodes3d", "get_firstlast_grid3d")
        object.add_getter("nodes3d", "get_node3d_position", names=("lat", "lon", "z"))
        object.add_getter("nodes3d", "get_node3d_velocity_xvel_curtime", names=("xvel_all",))
        object.add_getter("nodes3d", "get_node3d_velocity_xvel_curtime", names=("xvel_cur",))
        object.add_getter("nodes3d", "get_node3d_velocity_xvel_oldtime", names=("xvel_old",))
        object.add_getter("nodes3d", "get_node3d_velocity_yvel_curtime", names=("yvel_all",))
        object.add_getter("nodes3d", "get_node3d_velocity_yvel_curtime", names=("yvel_cur",))
        object.add_getter("nodes3d", "get_node3d_velocity_yvel_oldtime", names=("yvel_old",))
        object.add_getter("nodes3d", "get_node3d_velocity_zvel", names=("zvel",))
        object.add_setter("nodes3d", "set_node3d_velocity_xvel_alltime", names=("xvel_all",))
        object.add_setter("nodes3d", "set_node3d_velocity_xvel_curtime", names=("xvel_cur",))
        object.add_setter("nodes3d", "set_node3d_velocity_xvel_oldtime", names=("xvel_old",))
        object.add_setter("nodes3d", "set_node3d_velocity_yvel_alltime", names=("yvel_all",))
        object.add_setter("nodes3d", "set_node3d_velocity_yvel_curtime", names=("yvel_cur",))
        object.add_setter("nodes3d", "set_node3d_velocity_yvel_oldtime", names=("yvel_old",))

        # these are all on the U-grid
        object.define_grid("forcings", axes_names=axes_names, grid_class=StructuredGrid)
        object.set_grid_range("forcings", "get_firstlast_node")
        object.add_getter("forcings", "get_node_coriolis_f", names=("coriolis_f",))
        object.add_setter("forcings", "set_node_coriolis_f", names=("coriolis_f",))
        object.add_getter("forcings", "get_node_wind_stress", names=("tau_x", "tau_y"))
        object.add_setter("forcings", "set_node_wind_stress", names=("tau_x", "tau_y"))
        object.add_getter("forcings", "get_node_position", names=("lat", "lon"))

        # elements are on the T-grid
        object.define_grid("elements", axes_names=axes_names, grid_class=StructuredGrid)
        object.set_grid_range("elements", "get_firstlast_node")
        object.add_getter("elements", "get_element_position", names=("lat", "lon"))
        object.add_getter("elements", "get_element_depth", names=("depth",))
        object.add_getter(
            "elements", "get_element_surface_state", names=("temperature", "salinity")
        )
        object.add_getter(
            "elements", "get_element_surface_heat_flux", names=("surface_heat_flux",)
        )
        object.add_getter("elements", "get_element_ssh", names=("ssh",))
        object.add_getter("elements", "get_element_ssh_curtime", names=("cur_ssh",))
        object.add_getter("elements", "get_element_ssh_oldtime", names=("old_ssh",))
        object.add_getter("elements", "get_element_ssh_guess", names=("guess_ssh",))
        object.add_setter("elements", "set_element_ssh", names=("ssh",))
        object.add_setter("elements", "set_element_ssh_curtime", names=("cur_ssh",))
        object.add_setter("elements", "set_element_ssh_oldtime", names=("old_ssh",))
        object.add_setter("elements", "set_element_ssh_guess", names=("guess_ssh",))

        # elements are on the T-grid
        object.define_grid("elements3d")
        object.set_grid_range("elements3d", "get_firstlast_grid3d")
        object.add_getter(
            "elements3d", "get_element3d_position", names=("lat", "lon", "z")
        )
        object.add_getter(
            "elements3d", "get_element3d_temperature", names=("temperature",)
        )
        object.add_getter("elements3d", "get_element3d_salinity", names=("salinity",))
        object.add_getter("elements3d", "get_element3d_density", names=("rho",))
        object.add_setter(
            "elements3d", "set_element3d_temperature", names=("temperature",)
        )
        object.add_setter("elements3d", "set_element3d_salinity", names=("salinity",))
        object.add_setter("elements3d", "set_element3d_density", names=("rho",))

        # these are all on the T-grid
        object.define_grid(
            "element_forcings", axes_names=axes_names, grid_class=StructuredGrid
        )
        object.set_grid_range("element_forcings", "get_firstlast_node")
        object.add_getter(
            "element_forcings",
            "get_element_surface_forcing_temp",
            names=("restoring_temp",),
        )
        object.add_setter(
            "element_forcings",
            "set_element_surface_forcing_temp",
            names=("restoring_temp",),
        )
        object.add_getter(
            "element_forcings",
            "get_element_surface_forcing_salt",
            names=("restoring_salt",),
        )
        object.add_setter(
            "element_forcings",
            "set_element_surface_forcing_salt",
            names=("restoring_salt",),
        )
        object.add_getter(
            "element_forcings", "get_element_position", names=("lat", "lon")
        )

    @property
    def Tgrid(self):
        return self.elements3d

    @property
    def Ugrid(self):
        return self.nodes3d

    def define_errorcodes(self, handler):
        handler.add_errorcode(-1, "interface function returned -1, general error")
        handler.add_errorcode(-2, "POP internal error")
        handler.add_errorcode(
            -11, "evolve model aborted due to runaway kinetic energy, Ek>100"
        )

    def define_parameters(self, object):
        CodeWithNamelistParameters.define_parameters(self, object)

        object.add_interface_parameter(
            "namelist_file",
            "configuration (namelist) file with simulation setup",
            self._nml_file,
            state_guard="configuration_file_set",
        )

        object.add_method_parameter(
            "get_horiz_grid_option",
            "set_horiz_grid_option",
            "horiz_grid_option",
            "Option for horizontal grid should be either 'internal', 'amuse' or 'file'",
            default_value="internal",
        )
        object.add_method_parameter(
            "get_horiz_grid_file",
            "set_horiz_grid_file",
            "horiz_grid_file",
            "Filename for specifying the horizontal grid latitudes and longitudes",
            default_value="",
        )

        object.add_method_parameter(
            "get_vert_grid_option",
            "set_vert_grid_option",
            "vert_grid_option",
            "Option for vertical grid should be either 'internal', 'amuse' or 'file'",
            default_value="internal",
        )
        object.add_method_parameter(
            "get_vert_grid_file",
            "set_vert_grid_file",
            "vert_grid_file",
            "Filename for specifying vertical grid layer depths",
            default_value="",
        )

        object.add_method_parameter(
            "get_topography_option",
            "set_topography_option",
            "topography_option",
            "Option for topography should be either 'amuse','internal' or 'file'",
            default_value="internal",
        )

        object.add_interface_parameter(
            "depth_index",
            "single value or array of depth indices, determines bathymetry in case topography_option=amuse together with vertical_layer_thicknesses",
            default_value=0,
        )

        object.add_method_parameter(
            "get_topography_file",
            "set_topography_file",
            "topography_file",
            "Filename for specifying the topography in number of depths levels per gridpoint",
            default_value="",
        )

        object.add_method_parameter(
            "get_bottom_cell_file",
            "set_bottom_cell_file",
            "bottom_cell_file",
            "Filename for specifying the partial bottom cells, setting this enables partial bottom cells",
            default_value="",
        )

        object.add_method_parameter(
            "get_ts_option",
            "set_ts_option",
            "ts_option",
            "Option for the initial source for temperature and salinity",
            default_value="internal",
        )

        object.add_method_parameter(
            "get_ts_file",
            "set_ts_file",
            "ts_file",
            "Path to filename for initial source of temperature and salinity",
            default_value="",
        )
        object.add_method_parameter(
            "get_ts_file_format",
            "set_ts_file_format",
            "ts_file_format",
            "File format for the file specified in ts_file, can be either 'bin' or 'nc'",
            default_value="nc",
        )

        object.add_method_parameter(
            "get_distribution",
            "set_distribution",
            "distribution",
            "The distribution to use to distribution the decomposed domain among MPI workers, supported now: 'cartesian' and 'predefined'",
            default_value="cartesian",
        )
        object.add_method_parameter(
            "get_distribution_file",
            "set_distribution_file",
            "distribution_file",
            "The path to the filename for the predefined distribution",
            default_value="",
        )
        object.add_method_parameter(
            "get_ew_boundary_type",
            "set_ew_boundary_type",
            "ew_boundary_type",
            "The type of the east-west boundary, can be either 'closed' or 'cyclic'",
            default_value="cyclic",
        )
        object.add_method_parameter(
            "get_ns_boundary_type",
            "set_ns_boundary_type",
            "ns_boundary_type",
            "The type of the north-south boundary, can be any of 'closed', 'cyclic', or 'tripole'",
            default_value="closed",
        )

        object.add_method_parameter(
            "get_restart_option",
            "set_restart_option",
            "restart_option",
            "Specifies the quantity used in freq_option, can be any of 'never', 'nhour', 'nday', 'nmonth', 'nyear'",
            default_value="never",
        )
        object.add_method_parameter(
            "get_restart_freq_option",
            "set_restart_freq_option",
            "restart_freq_option",
            "The frequency with which restart files are written in days, or months, etc. depending on restart_option",
            default_value=1,
        )
        object.add_method_parameter(
            "get_restart_file",
            "set_restart_file",
            "restart_file",
            "The filename of where the restart files should be written, the runid and data will be appended to the filename",
            default_value="",
        )
        object.add_method_parameter(
            "get_tavg_option",
            "set_tavg_option",
            "tavg_option",
            "Specifies the quantity used in freq_option, can be any of 'never', 'nhour', 'nday', 'nmonth', 'nyear'",
            default_value="never",
        )
        object.add_method_parameter(
            "get_tavg_freq_option",
            "set_tavg_freq_option",
            "tavg_freq_option",
            "The frequency with which tavg output files are written in days, or months, etc. depending on restart_option",
            default_value=1,
        )
        object.add_method_parameter(
            "get_tavg_file",
            "set_tavg_file",
            "tavg_file",
            "The filename of where the tavg output files should be written, the runid and data will be appended to the filename",
            default_value="",
        )
        object.add_method_parameter(
            "get_movie_option",
            "set_movie_option",
            "movie_option",
            "Specifies the quantity used in freq_option, can be any of 'never', 'nhour', 'nday', 'nmonth', 'nyear'",
            default_value="never",
        )
        object.add_method_parameter(
            "get_movie_freq_option",
            "set_movie_freq_option",
            "movie_freq_option",
            "The frequency with which movie output files are written in days, or months, etc. depending on restart_option",
            default_value=1,
        )
        object.add_method_parameter(
            "get_movie_file",
            "set_movie_file",
            "movie_file",
            "The filename of where the movie output files should be written, the runid and data will be appended to the filename",
            default_value="",
        )

        object.add_method_parameter(
            "get_runid",
            "set_runid",
            "runid",
            "The name of the current run used for naming the output files, default is 'AMUSE'",
            default_value="AMUSE",
        )
        object.add_method_parameter(
            "get_dt_option",
            "set_dt_option",
            "dt_option",
            "The option used for setting the dt_count, can be 'auto_dt', 'steps_per_year', 'steps_per_day', 'seconds', 'hours'",
            default_value="steps_per_day",
        )
        object.add_method_parameter(
            "get_dt_count",
            "set_dt_count",
            "dt_count",
            "The amount of seconds, hours, or steps per day or year, depending on the value of dt_option",
            default_value=45,
        )

        object.add_method_parameter(
            "get_shf_data_type",
            "set_shf_data_type",
            "surface_heat_flux_forcing",
            "Option for surface heat flux should be either 'analytic', 'amuse' or 'file'",
            ##"Setting for surface heat flux",
            default_value="none",
        )
        object.add_method_parameter(
            "get_sfwf_data_type",
            "set_sfwf_data_type",
            "surface_freshwater_flux_forcing",
            "Options for surface salinity flux should be either 'analytic', 'amuse' or 'file'",
            # "Setting for surface freshwater flux forcing",
            default_value="none",
        )
        object.add_method_parameter(
            "get_ws_data_type",
            "set_ws_data_type",
            "windstress_forcing",
            "Setting for surface wind stress forcing",
            default_value="none",
        )
        object.add_method_parameter(
            "get_shf_filename",
            "set_shf_monthly_file",
            "surface_heat_flux_monthly_file",
            "Input filename for surface heat flux forcing",
            default_value="",
        )
        object.add_method_parameter(
            "get_sfwf_filename",
            "set_sfwf_monthly_file",
            "surface_freshwater_flux_monthly_file",
            "Input filename for surface freshwater flux forcing",
            default_value="",
        )
        object.add_method_parameter(
            "get_ws_filename",
            "set_ws_monthly_file",
            "windstress_monthly_file",
            "Input filename for wind stress forcing",
            default_value="",
        )
        object.add_method_parameter(
            "get_fwf_imposed",
            "set_fwf_imposed",
            "fwf_imposed",
            "Specifies the annual amount of imposed fresh water flux (in Sverdrups)",
            default_value=0.0 | units.Sv,
        )

        # ~ object.add_method_parameter(
        # ~ "get_namelist_filename",
        # ~ "set_namelist_filename",
        # ~ "namelist_filename",
        # ~ "Input filename used by code for reading the settings",
        # ~ default_value = 'amuse.nml'
        # ~ )

        object.add_interface_parameter(
            "vertical_layer_thicknesses",
            "input layer thicknesses (in case topography_opt==amuse) ",
            [] | units.cm,
        )

        object.add_method_parameter(
            "get_lonmin",
            "set_lonmin",
            "lonmin",
            "west boundary of domain (deg) in case horiz_grid_option=amuse",
            default_value=10.0 | units.deg,
        )

        object.add_method_parameter(
            "get_lonmax",
            "set_lonmax",
            "lonmax",
            "east boundary of domain (deg) in case horiz_grid_option=amuse",
            default_value=70.0 | units.deg,
        )

        object.add_method_parameter(
            "get_latmin",
            "set_latmin",
            "latmin",
            "south boundary of domain (deg) in case horiz_grid_option=amuse",
            default_value=10.0 | units.deg,
        )

        object.add_method_parameter(
            "get_latmax",
            "set_latmax",
            "latmax",
            "north boundary of domain (deg) in case horiz_grid_option=amuse",
            default_value=70.0 | units.deg,
        )

        object.add_method_parameter(
            "get_reinit_gradp",
            "set_reinit_gradp",
            "reinit_gradp",
            "flag indicates whether gradp terms are calculated from input Psurf on recommit_grids",
            default_value=False,
        )

        object.add_method_parameter(
            "get_reinit_rho",
            "set_reinit_rho",
            "reinit_rho",
            "flag indicates whether rho are calculated from input salt and temperature on recommit_grids",
            default_value=False,
        )
