from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.pop.interface import POPInterface
from amuse.community.pop.interface import POP

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
        time = instance.get_model_time()['time']
        dt = instance.get_timestep()['dt']
        print 'time=', time, 'dt=', dt

        tend = time + dt
        instance.evolve_model(tend)

        tend = tend + dt
        instance.evolve_model(tend)

        instance.cleanup_code()
        instance.stop()

    def test4(self):
        instance = POPInterface(**default_options)
        instance.initialize_code()

        #get domain size
        size = instance.get_domain_size()
        print size
        print 'nx_global=', size['nx'], 'ny_global=', size['ny']

        #get domain size
        num = instance.get_number_of_nodes()
        print num
        print 'n_nodes=', num['n_nodes']

        self.assertEquals(num['n_nodes'], size['nx']*size['ny'])

        instance.stop()
    
    
    def test5(self):
        instance = POPInterface(redirection="none")
        instance.initialize_code()

        #get the position of node with global index i=3, j=5
        print 'begin'
        for i in range(1, 10):
            for j in range(1, 10):
                pos = instance.get_node_position(i, j)
                print '(', pos['lat'], pos['lon'], ')',
            print ''
        print 'end'

        pos = instance.get_node_position(203, 53)
        print '(', pos['lat'], pos['lon'], ')'
        print 'lat=', pos['lat'], 'lon=', pos['lon']

        instance.stop()
