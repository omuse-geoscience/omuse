from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.pop.interface import POPInterface
from omuse.community.pop.interface import POP

from amuse.units import units

default_options=dict(redirection="none")

class POPInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = POPInterface(**default_options)
        print instance
        instance.stop()

    def test2(self):
        instance = POPInterface(**default_options)
        instance.initialize_code()
        instance.stop()

    #test whether we can evolve the model 2 timesteps, calling evolve_model twice
    def test3(self):
        instance = POPInterface(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        time_start = instance.get_model_time()['time']
        dt = instance.get_timestep_next()['dt']
        print 'time=', time_start, 'dt=', dt

        tend = time_start + dt
        instance.evolve_model(tend)

        dt = instance.get_timestep_next()['dt'] #dt may be different each time step because of halfsteps

        time = instance.get_model_time()['time']
        print 'time_start=', time_start, ' time=', time, ' dt=', dt, ' tend=', tend
        self.assertTrue(time > time_start, msg='Model time has not advanced despite evolve_model() step 1')
        self.assertTrue(time > tend-dt*0.5, msg='Model time has not advanced sufficiently despite evolve_model() step 1')
        self.assertTrue(time <= tend+dt*0.5, msg='Model time has advanced too much due to evolve_model() step 1')

        time_start = instance.get_model_time()['time']

        tend = tend + dt
        instance.evolve_model(tend)

        dt = instance.get_timestep_next()['dt'] #dt may be different each time step because of halfsteps
        time = instance.get_model_time()['time']
        self.assertTrue(time > time_start, msg='Model time has not advanced despite evolve_model() step 1')
        self.assertTrue(time > tend-dt*0.5, msg='Model time has not advanced sufficiently despite evolve_model() step 1')
        self.assertTrue(time <= tend+dt*0.5, msg='Model time has advanced too much due to evolve_model() step 1')

        instance.cleanup_code()
        instance.stop()

    def test4(self):
        instance = POPInterface(**default_options)
        instance.initialize_code()
        instance.commit_parameters()

        #get domain size
        size = instance.get_domain_size()
        print size
        print 'nx_global=', size['nx'], 'ny_global=', size['ny']

        #get domain size
        num = instance.get_number_of_nodes()
        print num
        print 'n_nodes=', num['n_nodes']

        self.assertEquals(num['n_nodes'], size['nx']*size['ny'], msg='The number of grid points does not match domain size')

        instance.stop()
    

    def test5(self):
        instance = POPInterface(redirection="none")
        instance.initialize_code()
        instance.commit_parameters()

        """ In POP the T-grid and U grid are aligned as follows:
              U ------U-------U-------U-------U
              |       |       |       |       |
          T-------T-------T-------T-------T   |
          |   |   |   |   |   |   |   |   |   |
          |   U --|---U---|---U---|---U---|---U
          |       |       |       |       |   |
          T-------T-------T-------T-------T   |
          |       |       |       |       |   |
          |   U   |   U   |   U   |   U   |---U
          |       |       |       |       |   |
          T-------T-------T-------T-------T   |
          |       |       |       |       |   |
          |   U   |   U   |   U   |   U   |---U
          |       |       |       |       |   |
          T-------T-------T-------T-------T   |
          |       |       |       |       |   |
          |   U   |   U   |   U   |   U   |---U
          |       |       |       |       |
          T-------T-------T-------T-------T

          This test checks whether the get_node_position and 
          get_element_position work correctly and whether the
          T and U grid are correctly positioned.

          It is assumed that the grid starts in the bottom
          south-west corner at global indexes 1,1 and ends at
          the top north-east corder with global indexes
          nx_global,ny_global

        """

        #get the position of the T grid in the lower left box
        sw = instance.get_element_position(1, 1)
        se = instance.get_element_position(2, 1)
        nw = instance.get_element_position(1, 2)
        ne = instance.get_element_position(2, 2)
        
        #check if these are legal latitude longitude values in radians
        range_check = sw['lon'] > -2*numpy.pi and sw['lon'] < 2*numpy.pi
        self.assertTrue(range_check, msg='longitude of point 1,1 not within expected range for radians (-2*PI and 2*PI)')
        range_check = sw['lat'] > -2*numpy.pi and sw['lat'] < 2*numpy.pi
        self.assertTrue(range_check, msg='latitude of point 1,1 not within expected range for radians (-2*PI and 2*PI)')

        #check if grid is increasing from west to east and
        #south to north
        increasing_we = (se['lon'] - sw['lon']) > 0
        increasing_sn = (nw['lat'] - sw['lat']) > 0
        self.assertTrue(increasing_we, msg='longitude appears to be non-increasing in west-east direction')
        self.assertTrue(increasing_sn, msg='latitude appears to be non-increasing in south-north direction')

        #check if point 1,1 on U grid is within (1,1), (1,2), (2,1), and (2,2) on the T grid
        u = instance.get_node_position(1, 1)
        u_correct = (sw['lon'] < u['lon'] and se['lon'] > u['lon'])
        u_correct = u_correct and (nw['lon'] < u['lon'] and ne['lon'] > u['lon'])
        u_correct = u_correct and (sw['lat'] < u['lat'] and nw['lat'] > u['lat'])
        u_correct = u_correct and (se['lat'] < u['lat'] and ne['lat'] > u['lat'])

        self.assertTrue(increasing_sn, msg='u point with index 1,1 does not seem to be embedded by t cell at 1,1')

        instance.stop()



