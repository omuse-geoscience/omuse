import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option


class AdcircInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
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


class Adcirc(CommonCode):
    def __init__(self, **options):
        # option for carthesian/ spherical coordinates
        CommonCode.__init__(self,  AdcircInterface(**options), **options)

    def get_firstlast_node(self):
        return 1,self.get_number_of_nodes()
    def get_firstlast_element(self):
        return 1,self.get_number_of_elements()
    def get_firstlast_node_of_elevation_boundary_segment(self,index_of_segment):
        return 1,self.get_number_of_nodes_in_elevation_boundary_segment(index_of_segment)
    def get_firstlast_node_of_flow_boundary_segment(self,index_of_segment):
        return 1,self.get_number_of_nodes_in_flow_boundary_segment(index_of_segment)
  
    def define_particle_sets(self, object):
        object.define_grid('nodes',axes_names = ['x','y'])
        object.set_grid_range('nodes', 'get_firstlast_node')    
        object.add_getter('nodes', 'get_node_state', names=('eta','vx','vy'))
        object.add_getter('nodes', 'get_node_position', names=('x','y'))
  
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

