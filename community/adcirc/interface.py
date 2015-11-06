import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

from amuse.units import trigo

from write_grid import adcirc_grid_writer,adcirc_parameter_writer

class AdcircInterface(CodeInterface, 
                      CommonCodeInterface,
                      StoppingConditionInterface,
                      LiteratureReferencesMixIn):
    """
    
    ADCIRC - ADvanced CIRCulation model

    .. [#] http://adcirc.org/
    
    """
    use_modules=['StoppingConditions','adcirc_interface']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'adcirc_worker'

    @remote_function
    def get_rootdir():
        returns (rootdir="s")
    @remote_function
    def set_rootdir(rootdir="."):
        returns ()

    @remote_function
    def commit_grid():
        pass

    @remote_function
    def get_model_time():
        returns (time=0.| units.s)
        
    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)
        
    @remote_function
    def evolve_model(tend=0. | units.s):
        pass

    @remote_function(can_handle_array=True)
    def get_node_eta(index=0):
        returns (eta=0. | units.m)
    @remote_function(can_handle_array=True)
    def set_node_eta(index=0,eta=0. | units.m):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_vx(index=0):
        returns (vx=0.| units.m/units.s)
    @remote_function(can_handle_array=True)
    def set_node_vx(index=0,vx=0.| units.m/units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_vy(index=0):
        returns (vy=0.| units.m/units.s)
    @remote_function(can_handle_array=True)
    def set_node_vy(index=0,vy=0.| units.m/units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_eta_prev(index=0):
        returns (eta_prev=0. | units.m)
    @remote_function(can_handle_array=True)
    def set_node_eta_prev(index=0,eta_prev=0. | units.m):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_status(index=0):
        returns (status='s')
    @remote_function(can_handle_array=True)
    def set_node_status(index=0,status='wet'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_coriolis_f(index=0):
        returns (coriolis_f=0. | units.s**-1)
    @remote_function(can_handle_array=True)
    def set_node_coriolis_f(index=0,coriolis_f=0. | units.s**-1):
        returns ()
        
    @remote_function(can_handle_array=True)
    def get_node_wind_stress(index=0):
        returns (tau_x=0. | units.Pa,tau_y=0. | units.Pa)
    @remote_function(can_handle_array=True)
    def set_node_wind_stress(index=0,tau_x=0.| units.Pa,tau_y=0. | units.Pa):
        returns ()        

    @remote_function(can_handle_array=True)
    def get_node_position(index='i'):
        returns (x=0.| units.m,y=0. | units.m)
    @remote_function(can_handle_array=True)
    def get_node_coordinates(index='i'):
        returns (lon=0.| units.rad,lat=0. | units.rad)

    @remote_function(can_handle_array=True)
    def get_node_depth(index='i'):
        returns (depth=0.| units.m)

    @remote_function(can_handle_array=True)
    def get_node_sigma(index='i',zindex='i'):
        returns (sigma='d',z='d' | units.m)

    @remote_function(can_handle_array=True)
    def get_element_position(index='i'):
        returns (x=0.| units.m,y=0. | units.m)
    @remote_function(can_handle_array=True)
    def get_element_coordinates(index='i'):
        returns (lon=0.| units.rad,lat=0. | units.rad)

    @remote_function(can_handle_array=True)
    def get_node_velocities_3d(index='i',zindex='i'):
        returns (vx=0.| units.m/units.s, vy=0. | units.m/units.s, vz=0. | units.m/units.s)

    @remote_function(can_handle_array=True)
    def get_element_nodes(index=0):
        returns (n1=0,n2=0,n3=0)

    @remote_function(can_handle_array=True)
    def get_element_status(index=0):
        returns (status='s')
    @remote_function(can_handle_array=True)
    def set_element_status(index=0,status='wet'):
        returns ()

    @remote_function
    def get_number_of_nodes():
        returns (n_nodes=0)

    @remote_function
    def get_number_of_elements():
        returns (n_elements=0)
    
    @remote_function
    def get_number_of_elevation_boundary_segments():
        returns (n_elev_boundaries=0)
   
    @remote_function
    def get_number_of_nodes_in_elevation_boundary_segment(index_of_segment=0):
        returns (n_nodes=0)   

    @remote_function
    def get_number_of_vertical_nodes():
        returns (n_znodes=0)

    @remote_function(can_handle_array=True)
    def get_elevation_boundary_node(index=0,index_of_segment=0):
        returns (node_index=0)

    @remote_function(can_handle_array=True)
    def set_elevation_boundary_eta(index=0,eta=0. | units.m,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_elevation_boundary_eta(index=0,index_of_segment=0):
        returns (eta=0. | units.m)

    @remote_function
    def get_number_of_flow_boundary_segments():
        returns (n_elev_boundaries=0)
   
    @remote_function
    def get_number_of_nodes_in_flow_boundary_segment(index_of_segment=0):
        returns (n_nodes=0)   

    @remote_function(can_handle_array=True)
    def get_flow_boundary_node(index=0,index_of_segment=0):
        returns (node_index=0)

    @remote_function(can_handle_array=True)
    def get_flow_boundary_type(index=0,index_of_segment=0):
        returns (node_type=0)

    @remote_function
    def get_use_interface_elevation_boundary():
        returns (use_interface_elevation_boundary='b')
    @remote_function
    def set_use_interface_elevation_boundary(use_interface_elevation_boundary='b'):
        returns ()

    @remote_function
    def get_use_interface_met_forcing():
        returns (use_interface_met_forcing='b')
    @remote_function
    def set_use_interface_met_forcing(use_interface_met_forcing='b'):
        returns ()

class Adcirc(CommonCode):
  
    MODE_2D = "2D"
    MODE_3D = "3D"
    
    def __init__(self, mode=MODE_2D, coordinates="cartesian", **options):
        self.mode=mode
        self.coordinates=coordinates
        CommonCode.__init__(self,  AdcircInterface(**options), **options)
        self._nodes=None
        self._elements=None
        self._elev_boundary=None
        self._flow_boundary=None
        self._parameters=None
        self.stopping_conditions = StoppingConditions(self)
       
    def assign_grid_and_boundary(self,nodes,elements,elev_boundary, flow_boundary):
        self._nodes=nodes
        self._elements=elements
        self._elev_boundary=elev_boundary
        self._flow_boundary=flow_boundary

    def commit_parameters(self):        
        if self.parameters.use_interface_parameters:
          param=adcirc_parameter_writer()
          if self._parameters is not None:
            param.parameters=self._parameters
          else:
            self._parameters=param.parameters
          if self.parameters.bottom_friction_law not in ["linear","quadratic","hybrid"]:
            raise Exception("invalid/ unimplemented bottom friction law: %s"%bottom_friction_law)
          param.update( ESLM=self.parameters.A_H.value_in(units.m**2/units.s),
                        SLAM0=trigo.in_deg(self.parameters.central_longitude),
                        SFEA0=trigo.in_deg(self.parameters.central_latitude),
                        ICS=1 if self.coordinates=="cartesian" else 2,
                        DTDP=(-1 if self.parameters.use_predictor_corrector else 1)*self.parameters.timestep.value_in(units.s),
                        NOLIBF=["linear","quadratic","hybrid"].index(self.parameters.bottom_friction_law),
                        TAU0=self.parameters.GWCE_weighting_factor,
                        TAU=self.parameters.linear_bottom_friction_coeff.value_in(units.s**-1),
                        CF=self.parameters.quadratic_bottom_friction_coeff,
                        STATIM=self.parameters.begin_time.value_in(units.day) )
          param.write()
        if self.parameters.use_interface_grid:
          adcirc_grid_writer(coordinates=self.coordinates).write_grid(self._nodes,self._elements, 
            self._elev_boundary,self._flow_boundary)
        self.overridden().commit_parameters()

    def get_firstlast_node(self):
        return 1,self.get_number_of_nodes()
    def get_firstlast_element(self):
        return 1,self.get_number_of_elements()
    def get_firstlast_node_of_elevation_boundary_segment(self,index_of_segment):
        return 1,self.get_number_of_nodes_in_elevation_boundary_segment(index_of_segment)
    def get_firstlast_node_of_flow_boundary_segment(self,index_of_segment):
        return 1,self.get_number_of_nodes_in_flow_boundary_segment(index_of_segment)
    def get_firstlast_vertical_index(self):
        return 1,self.get_number_of_vertical_nodes()
    def get_firstlast_grid3d(self):
        return 1,self.get_number_of_nodes(),1,self.get_number_of_vertical_nodes()

    #~ def get_node_position_3d(self,index,zindex):
        #~ sigma,z=self.get_node_sigma(index,zindex)
        #~ x,y=self.get_node_position(index)
        #~ return x,y,z

    def define_parameters(self, object):
        object.add_default_form_parameter(
            "rootdir", 
            "set the root directory", 
            "."
        )

        object.add_default_form_parameter(
            "use_interface_elevation_boundary", 
            "toggle the use of interface boundary conditions for the elevation specified boundaries", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_met_forcing", 
            "toggle the use of interface meteorological forcings", 
            False
        )
        object.add_interface_parameter(
            "A_H",
            "turbulent lateral friction coefficient",
            100.0 | units.m**2/units.s,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "timestep",
            "ADCIRC timestep",
            360.0 | units.s,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "use_predictor_corrector",
            "flag for use of predictor corrector integrator",
            True,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "use_interface_parameters",
            "flag for use of interface parameters (i.e. write fort.15)",
            False,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "use_interface_grid",
            "flag for use of interface grid (i.e. write fort.14)",
            False,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "central_latitude",
            "central latitude used for CPP projection between x,y and lon,lat",
            0. | units.deg,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "central_longitude",
            "central longitude used for CPP projection between x,y and lon,lat",
            0. | units.deg,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "bottom_friction_law",
            "type of bottom friction law [linear,quadratic, hybrid]",
            "linear",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "GWCE_weighting_factor",
            "factor in GWCE Weighting primitive and wave contributions [-5,-4,-3,-2,-1, or >0 suggested: 0.005-0.1]",
            0.005,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "linear_bottom_friction_coeff",
            "linear bottom friction coefficient [1/s], only used for linear bottom friction law",
            0. | units.s**-1,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "quadratic_bottom_friction_coeff",
            "quadratic bottom friction coefficient [dimensionless], only used for quadratic bottom friction law",
            0.,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "begin_time",
            "begin time of the simulation",
            0. | units.day,
            "before_set_interface_parameter"
        )
           
    def define_properties(self, object):
        object.add_property('get_model_time', public_name = "model_time")


    def get_grid(self):
        from amuse.datamodel.grids import new_unstructured_grid

        num_elems = self.get_number_of_elements()
        num_nodes = self.get_number_of_nodes()

        index_i = range(1,num_elems+1)
        triangles = self.get_element_nodes(index_i)
        triangles = zip(triangles[0],triangles[1],triangles[2])
        triangles = numpy.array(triangles)-1

        grid_corner_lon = numpy.zeros(num_elems * 3, dtype=numpy.double)
        grid_corner_lat = numpy.zeros(num_elems * 3, dtype=numpy.double)

        pos = self.get_node_coordinates(range(1,num_nodes+1))
        lon = pos[0].value_in(units.rad)
        lat = pos[1].value_in(units.rad)

        i=0
        for triangle in triangles:
            n1,n2,n3 = triangle

            grid_corner_lat[i*3+0] = lat[n1]
            grid_corner_lat[i*3+1] = lat[n2]
            grid_corner_lat[i*3+2] = lat[n3]

            grid_corner_lon[i*3+0] = lon[n1]
            grid_corner_lon[i*3+1] = lon[n2]
            grid_corner_lon[i*3+2] = lon[n3]
            i+=1

        corners = [grid_corner_lon, grid_corner_lat]
        corners = numpy.array(corners)
        adcirc_grid = new_unstructured_grid(num_elems, 3, corners, axes_names=("lon","lat"))

        return adcirc_grid


    def define_particle_sets(self, object):
        object.define_grid('nodes',axes_names = ['x','y'], grid_class=datamodel.UnstructuredGrid)
        object.set_grid_range('nodes', 'get_firstlast_node')
        object.add_getter('nodes', 'get_node_position', names=('x','y'))
        object.add_getter('nodes', 'get_node_coordinates', names=('lon','lat'))          
        object.add_getter('nodes', 'get_node_depth', names=('depth',))
        object.add_getter('nodes', 'get_node_eta', names=('eta',))
        object.add_getter('nodes', 'get_node_vx', names=('vx',))
        object.add_getter('nodes', 'get_node_vy', names=('vy',))
        object.add_getter('nodes', 'get_node_eta_prev', names=('eta_prev',))
        object.add_getter('nodes', 'get_node_status', names=('status',))
        object.add_setter('nodes', 'set_node_eta', names=('eta',))
        object.add_setter('nodes', 'set_node_vx', names=('vx',))
        object.add_setter('nodes', 'set_node_vy', names=('vy',))
        object.add_setter('nodes', 'set_node_eta_prev', names=('eta_prev',))
        object.add_setter('nodes', 'set_node_status', names=('status',))

                
        if self.mode in [self.MODE_3D]:
            object.define_grid('grid3d',axes_names=['x','y','z'])
            object.set_grid_range('grid3d', 'get_firstlast_grid3d')
            object.add_getter('grid3d', 'get_node_sigma', names = ('sigma','z'))
            object.add_getter('grid3d', 'get_node_velocities_3d', names = ('wx','wy','wz'))
            object.add_getter('grid3d', 'get_node_position', names = ('x','y'))
            object.add_getter('grid3d', 'get_node_coordinates', names = ('lon','lat'))

            object.add_gridded_getter('nodes', 'get_node_sigma','get_firstlast_vertical_index', names = ('sigma','z'))
            object.add_gridded_getter('nodes', 'get_node_velocities_3d','get_firstlast_vertical_index', names = ('wx','wy','wz'))

        object.define_grid('forcings',axes_names = ['x','y'], grid_class=datamodel.UnstructuredGrid)
        object.set_grid_range('forcings', 'get_firstlast_node')
        object.add_getter('forcings', 'get_node_coriolis_f', names=('coriolis_f',))
        object.add_setter('forcings', 'set_node_coriolis_f', names=('coriolis_f',))
        object.add_getter('forcings', 'get_node_wind_stress', names=('tau_x','tau_y'))
        object.add_setter('forcings', 'set_node_wind_stress', names=('tau_x','tau_y'))
        object.add_getter('forcings', 'get_node_position', names=('x','y'))
        object.add_getter('forcings', 'get_node_coordinates', names=('lon','lat'))

        object.define_grid('elements',axes_names = ['x','y'], grid_class=datamodel.UnstructuredGrid)
        object.set_grid_range('elements', 'get_firstlast_element')    
        object.add_getter('elements', 'get_element_nodes', names=('n1','n2','n3'))
        object.add_getter('elements', 'get_element_position', names=('x','y'))
        object.add_getter('elements', 'get_element_coordinates', names=('lon','lat'))
        object.add_getter('elements', 'get_element_status', names=('status',))

          
    def elevation_boundaries(self):
        n=self.get_number_of_elevation_boundary_segments()
        for i in range(1,n+1):
            yield self._create_new_grid(self.specify_elevation_boundary_grid, index = i)

    def specify_elevation_boundary_grid(self, definition, index=1):
        definition.set_grid_range('get_firstlast_node_of_elevation_boundary_segment') 
        definition.add_getter('get_elevation_boundary_node', names=('node',))
        definition.add_getter('get_elevation_boundary_eta', names=('eta',))
        definition.add_setter('set_elevation_boundary_eta', names=('eta',))
        definition.define_extra_keywords({'index_of_segment':index})

    def flow_boundaries(self):
        n=self.get_number_of_flow_boundary_segments()
        for i in range(1,n+1):
            yield self._create_new_grid(self.specify_flow_boundary_grid, index = i)

    def specify_flow_boundary_grid(self, definition, index=1):
        definition.set_grid_range('get_firstlast_node_of_flow_boundary_segment') 
        definition.add_getter('get_flow_boundary_node', names=('node',))
        definition.add_getter('get_flow_boundary_type', names=('type',))
        definition.define_extra_keywords({'index_of_segment':index})

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_method('!UNINITIALIZED', 'before_get_parameter')
        object.add_method('!UNINITIALIZED', 'before_set_parameter')
        object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        object.add_transition('INITIALIZED','EDIT','commit_parameters')

        #~ object.set_initial_state('UNINITIALIZED')
        #~ object.add_transition('!STOPPED', 'END', 'cleanup_code')
        #~ object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        #~ object.add_transition('END', 'STOPPED', 'stop', False)
        #~ object.add_method('STOPPED', 'stop')        
 
        object.add_method('INITIALIZED', 'assign_grid_and_boundary')
        object.add_method('INITIALIZED', 'before_set_interface_parameter')
        object.add_method('INITIALIZED', 'set_rootdir')

        for state in ["RUN","EDIT","EVOLVED"]:
          object.add_method(state,"get_model_time")
          object.add_method(state,"get_timestep")          
          object.add_method(state,"get_number_of_nodes")
          object.add_method(state,"get_number_of_elements")
          object.add_method(state,"get_number_of_nodes_in_elevation_boundary_segment")
          object.add_method(state,"get_number_of_nodes_flow_boundary_segment")
          object.add_method(state,"get_number_of_vertical_nodes")
          object.add_method(state,"get_grid")

        object.add_transition('EDIT', 'RUN', 'commit_grid')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_method('EVOLVED', 'evolve_model')
