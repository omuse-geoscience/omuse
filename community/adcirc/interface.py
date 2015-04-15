import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

from write_grid import adcirc_grid_writer,adcirc_parameter_writer

class AdcircInterface(CodeInterface, 
                      CommonCodeInterface,
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
    def get_model_time():
        returns (time=0.| units.s)
        
    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)
        
    @remote_function
    def evolve_model(tend=0. | units.s):
        pass

    @remote_function(can_handle_array=True)
    def get_node_state(index=0):
        returns (eta=0. | units.m,vx=0.| units.m/units.s,vy=0.| units.m/units.s)

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
    def get_element_position(index='i'):
        returns (x=0.| units.m,y=0. | units.m)

    @remote_function(can_handle_array=True)
    def get_element_nodes(index=0):
        returns (n1=0,n2=0,n3=0)

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
    def __init__(self, **options):
        # option for carthesian/ spherical coordinates
        CommonCode.__init__(self,  AdcircInterface(**options), **options)
        self._nodes=None
        self._elements=None
        self._elev_boundary=None
        self._flow_boundary=None
        self._parameters=None

    def assign_grid_and_boundary(self,nodes,elements,elev_boundary, flow_boundary):
        self._nodes=nodes
        self._elements=elements
        self._elev_boundary=elev_boundary
        self._flow_boundary=flow_boundary

    def commit_parameters(self):        
        if self.parameters.use_interface_parameters:
          A_H=self.parameters.A_H
          dt=self.parameters.timestep
          use_precor=self.parameters.use_predictor_corrector
          param=adcirc_parameter_writer()
          if self._parameters is not None:
            param.parameters=self._parameters
          param.set_non_default(A_H, dt, use_precor)
          param.write()
        if self.parameters.use_interface_grid:
          adcirc_grid_writer().write_grid(self._nodes,self._elements, 
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

    def define_parameters(self, object):
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
            default_value = 100.0 | units.m**2/units.s
        ) 
        object.add_interface_parameter(
            "timestep",
            "ADCIRC timestep",
            default_value = 360.0 | units.s
        ) 
        object.add_interface_parameter(
            "use_predictor_corrector",
            "flag for use of predictor corrector integrator",
            default_value = True
        ) 
        object.add_interface_parameter(
            "use_interface_parameters",
            "flag for use of interface parameters (i.e. write fort.15)",
            default_value = False
        ) 
        object.add_interface_parameter(
            "use_interface_grid",
            "flag for use of interface grid (i.e. write fort.14)",
            default_value = False
        )

    def define_particle_sets(self, object):
        object.define_grid('nodes',axes_names = ['x','y'])
        object.set_grid_range('nodes', 'get_firstlast_node')
        object.add_getter('nodes', 'get_node_state', names=('eta','vx','vy'))
        object.add_getter('nodes', 'get_node_position', names=('x','y'))

        object.define_grid('forcings',axes_names = ['x','y'])
        object.set_grid_range('forcings', 'get_firstlast_node')
        object.add_getter('forcings', 'get_node_coriolis_f', names=('coriolis_f',))
        object.add_setter('forcings', 'set_node_coriolis_f', names=('coriolis_f',))
        object.add_getter('forcings', 'get_node_wind_stress', names=('tau_x','tau_y'))
        object.add_setter('forcings', 'set_node_wind_stress', names=('tau_x','tau_y'))
        object.add_getter('forcings', 'get_node_position', names=('x','y'))

        object.define_grid('elements',axes_names = ['x','y'])
        object.set_grid_range('elements', 'get_firstlast_element')    
        object.add_getter('elements', 'get_element_nodes', names=('n1','n2','n3'))
        object.add_getter('elements', 'get_element_position', names=('x','y'))
          
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

