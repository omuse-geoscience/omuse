from amuse.community import *

from amuse.community.interface.common import CommonCodeInterface, CommonCode

class POPInterface(CodeInterface):
    
    """

    POP - Parallel Ocean Program

    .. [#] https://github.com/NLeSC/eSalsa-POP

    """
    include_headers = ['worker_code.h']
    use_modules = ['pop_interface']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker=self.name_of_the_worker(), number_of_workers = 8, **keyword_arguments)

    def name_of_the_worker(self):
        return 'pop_worker'


    ##forcings getters and setters
    @remote_function(must_handle_array=True)
    def get_node_wind_stress(i=0,j=0):
        returns (tau_x=0. | units.Pa,tau_y=0. | units.Pa)
    @remote_function(must_handle_array=True)
    def set_node_wind_stress(i=0,j=0,tau_x=0.| units.Pa,tau_y=0. | units.Pa):
        returns ()

    @remote_function(must_handle_array=True)
    def get_node_coriolis_f(i=0,j=0):
        returns (coriolis_f=0. | units.s**-1)
    @remote_function(must_handle_array=True)
    def set_node_coriolis_f(i=0,j=0,coriolis_f=0. | units.s**-1):
        returns ()


    ##getters for node and element state
    @remote_function(must_handle_array=True)
    def get_node_surface_state(i=0,j=0):
        returns (ssh=0. | units.cm, vx=0. | units.m/units.s, vy=0. | units.m/units.s)

    @remote_function(must_handle_array=True)
    def get_element_surface_state(i=0,j=0):
        returns (temp=0. | units.C, salt=0. | units.g/units.kg) #salinity in gram per kg 
#        returns (salt=0. | units.g/units.kg, temp=0. | units.C) #salinity in gram per kg 
    


    ##these two return in units.m in adcirc but here they return latitude longitude in degrees
    @remote_function(must_handle_array=True)
    def get_node_position(i=0,j=0):
        returns (lat=0. | units.deg, lon=0. | units.deg)

    @remote_function(must_handle_array=True)
    def get_element_position(i=0,j=0):
        returns (lat=0. | units.deg, lon=0. | units.deg)

    @remote_function
    def get_number_of_nodes():
        returns (n_nodes=0)

    #get number of gridpoints in x-direction (west-east) and y-direction (south-north)
    @remote_function
    def get_domain_size():
        returns (nx=0, ny=0)

    #get the total number of vertical levels used
    @remote_function
    def get_number_of_vertical_levels():
        returns (km=0)    

    #returns the maximal depth at this position
    @remote_function(must_handle_array=True)
    def get_node_depth(i=0,j=0):
        returns (depth=0. | units.m)
    @remote_function(must_handle_array=True)
    def get_element_depth(i=0,j=0):
        returns (depth=0. | units.m)


    
        
    @remote_function
    def initialize_code():
        returns ()

    @remote_function
    def evolve_model(tend=0. | units.s):
        pass

    @remote_function
    def cleanup_code():
        returns ()

    @remote_function
    def get_model_time():
        returns (time=0.| units.s)

    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)

    @remote_function
    def get_timestep_next():
        returns (dt=0. | units.s)

    #facilitate state transitions
    @remote_function
    def prepare_parameters():
        returns ()
    @remote_function
    def recommit_parameters():
        returns ()

    def synchronize_model(self):
        self.recommit_parameters()
        self.prepare_parameters()
        









    
class POP(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  POPInterface(**options), **options)
    
    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('!UNINITIALIZED!INITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        #some parameters have to be derived before they can be read, therefore
        #reading anything should be preceded by a call to prepare_parameters
        object.add_transition('INITIALIZED', 'EDIT', 'prepare_parameters')
        object.add_transition('INITIALIZED', 'RUN', 'prepare_parameters')

        #you can only edit stuff in state EDIT
        object.add_method('EDIT','before_set_parameter')
        #you can only read stuff in states RUN and EDIT
        object.add_method('EDIT','before_get_parameter')
        object.add_method('RUN','before_get_parameter')
        #setting a parameter is the only way to endup in state EDIT from RUN
        object.add_transition('RUN','EDIT','before_set_parameter')
        object.add_transition('RUN','EDIT_FORCINGS','before_set_parameter')

        for state in ["RUN","EDIT"]:
            object.add_method(state, 'get_node_coriolis_f')
            object.add_method(state, 'get_node_wind_stress')
            object.add_method(state, 'get_node_position')
            object.add_method(state, 'get_element_position')
            object.add_method(state, 'get_node_position')
            object.add_method(state, 'get_element_surface_state')
            object.add_method(state, 'get_node_surface_state')
        object.add_method('EDIT', 'set_node_coriolis_f')
        object.add_method('EDIT_FORCINGS', 'set_node_wind_stress')

        #before we can run the model we need to recommit_parameters 
        object.add_transition('EDIT', 'RUN', 'recommit_parameters')
        object.add_transition('EDIT_FORCINGS', 'RUN', 'synchronize_model')

        #can only evolve from run
        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)

        #from state EVOLVED you may evolve again if you like
        object.add_method('EVOLVED','evolve_model')

        #prepare model for reading/writing after evolve
        object.add_transition('EVOLVED', 'RUN', 'prepare_parameters')
        object.add_transition('EVOLVED', 'EDIT', 'prepare_parameters')



    def get_firstlast_node(self):
        size = self.get_domain_size()
        return 1,size[0],1,size[1]
    def get_firstlast_element(self):
        return self.get_firstlast_node(self)

    def get_ugrid_latlon_range(self):
        start = self.get_node_position(1,1)
        size = self.get_domain_size()
        end = self.get_node_position(size[0],size[1])
        return start[0],end[0],start[1],end[1]

    def get_tgrid_latlon_range(self):
        start = self.get_element_position(1,1)
        size = self.get_domain_size()
        end = self.get_element_position(size[0],size[1])
        return start[0],end[0],start[1],end[1]


    def define_particle_sets(self, object):
        #for now we refer to the U grid as the nodes and T grid as the elements
        #the forcings can be on either grid, depends on what is preferred for coupling with adcirc I guess

        object.define_grid('nodes', axes_names = ['x','y'])
        object.set_grid_range('nodes', 'get_firstlast_node')
        object.add_getter('nodes', 'get_node_surface_state', names=('ssh','vx','vy'))
        object.add_getter('nodes', 'get_node_position', names=('lat','lon'))

        object.define_grid('forcings',axes_names = ['x','y'])
        object.set_grid_range('forcings', 'get_firstlast_node')
        object.add_getter('forcings', 'get_node_coriolis_f', names=('coriolis_f',))    #on the U-grid
        object.add_setter('forcings', 'set_node_coriolis_f', names=('coriolis_f',))
        object.add_getter('forcings', 'get_node_wind_stress', names=('tau_x','tau_y')) #currently on the T-grid, fix later
        object.add_setter('forcings', 'set_node_wind_stress', names=('tau_x','tau_y'))
        object.add_getter('forcings', 'get_node_position', names=('lat','lon'))

        object.define_grid('elements',axes_names = ['x','y'])
        object.set_grid_range('elements', 'get_firstlast_element')
     #   object.add_getter('elements', 'get_element_nodes', names=('n1','n2','n3'))   #adcirc specific
        object.add_getter('elements', 'get_element_position', names=('lat','lon'))
        object.add_getter('elements', 'get_element_surface_state', names=('temp','salt'))
















