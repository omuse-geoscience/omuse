import os
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.rfi.core import CodeInterface
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.parameter_tools import CodeWithIniFileParameters
from amuse.rfi.core import legacy_function,remote_function,LegacyFunctionSpecification
from amuse.rfi.tools.fortran_tools import FortranCodeGenerator
from amuse.support.literature import LiteratureReferencesMixIn

from amuse import datamodel

from omuse.units import units

class DFlowFMInterface(CodeInterface, LiteratureReferencesMixIn):
    """
    DFlowFM - D-Flow Flexible Mesh
    
    .. [#] https://www.deltares.nl/en/software/module/d-flow-flexible-mesh/
    
    """
    
    use_modules=["dflowfm_omuse"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="dflowfm_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
    
    @remote_function
    def initialize():
        returns ()

    @remote_function
    def commit_parameters():
        returns ()
        
    @remote_function
    def evolve_model(tend=0. | units.s):
        returns ()
    
    @remote_function
    def get_model_time():
        returns (model_time = 0. | units.s)

    @remote_function
    def get_use_patm():
        returns (use_patm=False)

    @remote_function
    def get_use_wind():
        returns (use_wind=False)

    @remote_function
    def get_use_waterlevel():
        returns (use_waterevel=False)

    @remote_function
    def set_use_patm(use_patm=False):
        returns ()

    @remote_function
    def set_use_wind(use_wind=False):
        returns ()

    @remote_function
    def set_use_waterlevel(use_waterlevel=False):
        returns ()
    
    @remote_function
    def get_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_1d_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_boundary_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_1d_boundary_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_global_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_global_internal_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_global_boundary_flow_nodes_range():
        returns (imin=0, imax=0)

    @remote_function(must_handle_array=True)
    def get_x_position(index=0):
        returns (x=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_y_position(index=0):
        returns (y=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_water_level(index=0):
        returns (h=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_ucx(index=0):
        returns (vx=0. | units.m/units.s)

    @remote_function(must_handle_array=True)
    def get_ucy(index=0):
        returns (vy=0. | units.m/units.s)

    @remote_function(must_handle_array=True)
    def get_patm(index=0):
        returns (p=0. | units.N/units.m**2)

    @remote_function(must_handle_array=True)
    def set_patm(index=0,p=0. | units.N/units.m**2):
        returns ()

    @remote_function
    def get_net_nodes_range():
        returns (imin=0, imax=0)

    @remote_function(must_handle_array=True)
    def get_x_position_net_nodes(index=0):
        returns (x=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_y_position_net_nodes(index=0):
        returns (y=0.) # units deg or m, set late in define_methods

    #~ @remote_function
    #~ def get_internal_flow_links_range():
        #~ returns (imin=0, imax=0)

    #~ @remote_function
    #~ def get_boundary_flow_links_range():
        #~ returns (imin=0, imax=0)

    @remote_function
    def get_global_flow_links_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_global_internal_flow_links_range():
        returns (imin=0, imax=0)

    @remote_function
    def get_global_boundary_flow_links_range():
        returns (imin=0, imax=0)

    @remote_function(must_handle_array=True)
    def get_x_position_flow_links(index=0):
        returns (x=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_y_position_flow_links(index=0):
        returns (y=0.) # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_wx(index=0):
        returns (wind_vx=0. | units.m/units.s) 

    @remote_function(must_handle_array=True)
    def get_wy(index=0):
        returns (wind_vy=0. | units.m/units.s) 

    @remote_function(must_handle_array=True)
    def get_zbndz(index=0):
        returns (water_level=0. | units.m) 

    @remote_function(must_handle_array=True)
    def set_wx(index=0,wind_vx=0. | units.m/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def set_wy(index=0,wind_vy=0. | units.m/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def set_zbndz(index=0,water_level=0. | units.m):
        returns ()

    @remote_function(must_handle_array=True)
    def get_is_waterlevel_bnd(index=0):
        returns (is_waterlevel_bnd="b") 

    @remote_function(must_handle_array=True)
    def get_xbndz(index=0):
        returns (xbndz=0. )  # units deg or m, set late in define_methods

    @remote_function(must_handle_array=True)
    def get_ybndz(index=0):
        returns (ybndz=0. )  # units deg or m, set late in define_methods


class DFlowFM(InCodeComponentImplementation, CodeWithIniFileParameters):

    def __init__(self, **options):
        self._workdir=options.get("workdir", ".")
        self._ini_file=os.path.join(self._workdir, options.get("ini_file",""))
        CodeWithIniFileParameters.__init__(self, options.get("parameters", dict()) ) 
        self._coordinates=options.get("coordinates", "cartesian") # should be automatic
        _cwd=os.getcwd()
        os.chdir(self._workdir)
        InCodeComponentImplementation.__init__(self,  DFlowFMInterface(**options), **options)
        os.chdir(_cwd)
        if self._ini_file:
            self.parameters.ini_file=self._ini_file
        
    def define_properties(self, handler):
        handler.add_property('get_model_time', public_name = "model_time")

    def configuration_file_set(self):
        self.read_inifile_parameters(self.parameters.ini_file, add_missing_parameters=True)
        handler=self.get_handler('PARAMETER')
        CodeWithIniFileParameters.define_parameters(self,handler)


    def define_parameters(self,handler):
        CodeWithIniFileParameters.define_parameters(self, handler)
      
        handler.add_interface_parameter(
            "ini_file",
            "configuration file with simulation setup",
            self._ini_file,
            state_guard="configuration_file_set"
        )

        handler.add_boolean_parameter(
            "get_use_wind", 
            "set_use_wind",
            "use_interface_wind", 
            "set wind field thorugh interface (True) or not (False)", 
            default_value = False
        )

        handler.add_boolean_parameter(
            "get_use_patm", 
            "set_use_patm",
            "use_interface_patm", 
            "set atmospheric pressure field thorugh interface (True) or not (False)", 
            default_value = False
        )

        handler.add_boolean_parameter(
            "get_use_waterlevel", 
            "set_use_waterlevel",
            "use_interface_waterlevel_boundary", 
            "set waterlevel boundaries through interface (True) or not (False) (still needs input ext describing boundaries)", 
            default_value = False
        )

    def define_grids(self, handler):
        if self._coordinates=="cartesian":
            axes_names=['x','y']
            coordinates="position"
        elif self._coordinates=="spherical":
            axes_names=['lon','lat']
            coordinates="lonlat"

        handler.define_grid('flow_nodes',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('flow_nodes', 'get_global_internal_flow_nodes_range')
        handler.add_getter('flow_nodes', 'get_x_position', names=axes_names[0:1])
        handler.add_getter('flow_nodes', 'get_y_position', names=axes_names[1:2])
        handler.add_getter('flow_nodes', 'get_water_level', names=["water_level"])
        handler.add_getter('flow_nodes', 'get_ucx', names=["vx"])
        handler.add_getter('flow_nodes', 'get_ucy', names=["vy"])

        handler.define_grid('boundary_nodes',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('boundary_nodes', 'get_global_boundary_flow_nodes_range')
        handler.add_getter('boundary_nodes', 'get_x_position', names=axes_names[0:1])
        handler.add_getter('boundary_nodes', 'get_y_position', names=axes_names[1:2])
        handler.add_getter('boundary_nodes', 'get_water_level', names=["water_level"])
        handler.add_getter('boundary_nodes', 'get_ucx', names=["vx"])
        handler.add_getter('boundary_nodes', 'get_ucy', names=["vy"])

        handler.define_grid('net_nodes',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('net_nodes', 'get_net_nodes_range')
        handler.add_getter('net_nodes', 'get_x_position_net_nodes', names=axes_names[0:1])
        handler.add_getter('net_nodes', 'get_y_position_net_nodes', names=axes_names[1:2])

        handler.define_grid('flow_links',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('flow_links', 'get_global_internal_flow_links_range')
        handler.add_getter('flow_links', 'get_x_position_flow_links', names=axes_names[0:1])
        handler.add_getter('flow_links', 'get_y_position_flow_links', names=axes_names[1:2])

        handler.define_grid('boundary_flow_links',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('boundary_flow_links', 'get_global_boundary_flow_links_range')
        handler.add_getter('boundary_flow_links', 'get_x_position_flow_links', names=axes_names[0:1])
        handler.add_getter('boundary_flow_links', 'get_y_position_flow_links', names=axes_names[1:2])

        # should be set later, when its known if these are enabled..
        handler.define_grid('flow_nodes_forcing',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('flow_nodes_forcing', 'get_global_flow_nodes_range')
        handler.add_getter('flow_nodes_forcing', 'get_x_position', names=axes_names[0:1])
        handler.add_getter('flow_nodes_forcing', 'get_y_position', names=axes_names[1:2])
        handler.add_getter('flow_nodes_forcing', 'get_patm', names=["atmospheric_pressure"])
        handler.add_setter('flow_nodes_forcing', 'set_patm', names=["atmospheric_pressure"])

        handler.define_grid('flow_links_forcing',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('flow_links_forcing', 'get_global_flow_links_range')
        handler.add_getter('flow_links_forcing', 'get_x_position_flow_links', names=axes_names[0:1])
        handler.add_getter('flow_links_forcing', 'get_y_position_flow_links', names=axes_names[1:2])
        handler.add_getter('flow_links_forcing', 'get_wx', names=["wind_vx"])
        handler.add_setter('flow_links_forcing', 'set_wx', names=["wind_vx"])
        handler.add_getter('flow_links_forcing', 'get_wy', names=["wind_vy"])
        handler.add_setter('flow_links_forcing', 'set_wy', names=["wind_vy"])

        handler.define_grid('boundary_links_forcing',axes_names = axes_names, 
                state_guard="before_new_set_instance", grid_class=datamodel.UnstructuredGrid)
        handler.set_grid_range('boundary_links_forcing', 'get_global_boundary_flow_links_range')
        handler.add_getter('boundary_links_forcing', 'get_x_position_flow_links', names=axes_names[0:1])
        handler.add_getter('boundary_links_forcing', 'get_y_position_flow_links', names=axes_names[1:2])
        handler.add_getter('boundary_links_forcing', 'get_zbndz', names=["water_level"])
        handler.add_setter('boundary_links_forcing', 'set_zbndz', names=["water_level"])
# adhoc attribute to know if a boundary point is waterlevel boundary
        handler.add_getter('boundary_links_forcing', 'get_is_waterlevel_bnd', names=["is_waterlevel_bnd"])
        handler.add_getter('boundary_links_forcing', 'get_xbndz', names=[axes_names[0]+"_bndz"])
        handler.add_getter('boundary_links_forcing', 'get_ybndz', names=[axes_names[1]+"_bndz"])

    def define_methods(self, handler):
        if self._coordinates=="spherical":
          handler.add_method( "get_x_position", (handler.INDEX,), ( units.deg, handler.ERROR_CODE,) )
          handler.add_method( "get_y_position", (handler.INDEX,), ( units.deg, handler.ERROR_CODE,) )
          handler.add_method( "get_x_position_flow_links", (handler.INDEX,), ( units.deg, handler.ERROR_CODE,) )
          handler.add_method( "get_y_position_flow_links", (handler.INDEX,), ( units.deg, handler.ERROR_CODE,) )
          handler.add_method( "get_xbndz", (handler.INDEX,), ( units.deg, handler.ERROR_CODE,) )
          handler.add_method( "get_ybndz", (handler.INDEX,), ( units.deg, handler.ERROR_CODE,) )
        elif self._coordinates=="cartesian":
          handler.add_method( "get_x_position", (handler.INDEX,), ( units.m, handler.ERROR_CODE,) )
          handler.add_method( "get_y_position", (handler.INDEX,), ( units.m, handler.ERROR_CODE,) )
          handler.add_method( "get_x_position_flow_links", (handler.INDEX,), ( units.m, handler.ERROR_CODE,) )
          handler.add_method( "get_y_position_flow_links", (handler.INDEX,), ( units.m, handler.ERROR_CODE,) )
          handler.add_method( "get_xbndz", (handler.INDEX,), ( units.m, handler.ERROR_CODE,) )
          handler.add_method( "get_ybndz", (handler.INDEX,), ( units.m, handler.ERROR_CODE,) )
        else:
          raise Exception("unknown coordinates")

    def commit_parameters(self):
        if self.channel.number_of_workers==1:
            self.write_inifile_parameters(os.path.join(self._workdir, "omuse.mdu"))
        else:
            self.write_multiple_inifile_parameters(os.path.join(self._workdir, "omuse.mdu"))
        self.overridden().commit_parameters()

    # convenience function to write multiple mdu files for parallel runs
    # this is a bit ad-hoc, better make something in parameter_tools?
    def write_multiple_inifile_parameters(self, outputfile):
        if not self.ini_geometry.PartitionFile:
                raise Exception("please set parameter ini_geometry.PartitionFile")
        if self.ini_numerics.Icgsolver not in [6,7]:
                raise Exception("please set ini_numerics.Icgsolver to 6 (PETSC) or 7 (GS)")
 
        orig_netfile=self.ini_geometry.NetFile
        orig_snapshotdir=self.ini_output.SnapshotDir
        orig_restartfile=self.ini_restart.RestartFile

        basename=outputfile.split('.')[0]
        netbase=orig_netfile.rsplit("_net.nc")[0]
        for i in range(self.channel.number_of_workers):
            n = "{:04d}".format(i)
            filename=basename + '_' + n + ".mdu"

            self.ini_geometry.NetFile=netbase + '_' + n + '_net.nc'
            self.ini_output.SnapshotDir='snapshots_' + n
            #~ self.ini_restart.RestartFile=basename + '_' + n + '_rst.nc'

            self.write_inifile_parameters(filename)
        
        self.ini_geometry.NetFile=orig_netfile
        self.ini_output.SnapshotDir=orig_snapshotdir
        self.ini_restart.RestartFile=orig_restartfile

    def define_state(self, handler):
        handler.set_initial_state('UNINITIALIZED')
        handler.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize')
        handler.add_transition('INITIALIZED', 'PARAM', 'commit_parameters')

        for method in ["get_global_flow_nodes_range", "get_global_internal_flow_nodes_range",
                       "get_global_boundary_flow_nodes_range",
                        "evolve_model",
                      ]:
            handler.add_method('PARAM', method)



partitioner_parameters={
    "filename" : dict(short="fnam", dtype="string", default="omuse", description="filename", ptype="simple"),
    "number_of_subdomains" : dict(short="md_Ndomains", dtype="int32", default=3, description="number of subdomains, Metis (>0) or polygon (0)", ptype="simple"),
    "contiguous_domains" : dict(short="md_jacontiguous", dtype="int32", default=0, description="contiguous domains, Metis (1) or not (0)", ptype="simple"),
    "intended_solver" : dict(short="md_icgsolver", dtype="int32", default=6, description="intended solver 6: PETSc (recommended), 7: parallel GS+CG", ptype="simple"),
    "partition_method" : dict(short="md_pmethod", dtype="int32", default=0, description="Partition method, 0=recursive bisection, 1=K-way", ptype="simple"),
    "dry_points_file" : dict(short="md_dryptsfile", dtype="string", default="none", description="optional dry points file", ptype="simple"),
    "enclosure_file" : dict(short="md_encfile", dtype="string", default="none", description="optional enclosure file to clip outer parts from the grid .pol", ptype="simple"),
    "generate_polygon_file" : dict(short="md_genpolygon", dtype="int32", default=0, description="make partition file (1) or not (0)", ptype="simple"),
    }

code_generator_partitioner=FortranCodeGenerator(partitioner_parameters)

class PartitionerInterface(CommonCodeInterface, CodeInterface):

    use_modules=["partitioner_omuse"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="partitioner_worker", **keyword_arguments)

    exec(code_generator_partitioner.generate_interface_functions())

    @remote_function
    def do_partition():
      returns ()

class Partitioner(CommonCode, InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  PartitionerInterface(**options), **options)

    def define_parameters(self, handler):
        code_generator_partitioner.generate_parameter_definitions(handler)

    def define_state(self, handler):
        CommonCode.define_state(self,handler)
        handler.add_method("RUN", "do_partition")
        handler.add_transition('INITIALIZED', 'RUN', 'commit_parameters')
