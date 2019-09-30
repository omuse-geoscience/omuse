from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from omuse.community.pop.interface import POPInterface
from omuse.community.pop.interface import POP

from amuse.units import units

#default_options=dict(channel_type="sockets",redirection="none")
default_options=dict(redirection="none")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from nose.tools import nottest



class POPTests(TestWithMPI):


    #test the behavior of the state machine
    @nottest
    def test1(self):
        instance = POP(**default_options)

        self.assertEquals(instance.state_machine._current_state.name, 'UNINITIALIZED')

        #a read of a parameter requires the state to be either EDIT or RUN, which means we pass through INITIALIZED
        fcor = instance.forcings.coriolis_f 

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'RUN')

        #check if we read a sensible value
        print 'fcor[1,1] = ', fcor[1,1]
        self.assertTrue(fcor[1,1] != 0 | units.s**-1, msg='Expected coriolis force to be not equal to zero for node 1,1')

        #proceed to evolve
        time = instance.get_model_time()
        tend = time + instance.get_timestep_next()
        instance.evolve_model(tend)

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'EVOLVED')

        #try to read something again, should cause a state transition to RUN
        fcor = instance.forcings.coriolis_f 

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'RUN')

        #check if we can write to coriolis_f
        instance.forcings.coriolis_f = fcor

        #check if we are in the expected state
        self.assertEquals(instance.state_machine._current_state.name, 'EDIT')
        
        instance.stop()



    @nottest
    def test2(self):
        instance = POP(**default_options)

        #proceed to evolve
        #time = instance.get_model_time()
        #tend = time + (1 | units.day)
        #instance.evolve_model(tend)

        #extract the tau_x and tau_y from POP
        tau_x = instance.forcings.tau_x
        tau_y = instance.forcings.tau_y

        print 'tau_x='
        print tau_x[range(1,5),1]
#        print tau_x
        print 'tau_y='
        print tau_y[range(1,5),1]
#        print tau_y

        #check if tau_x and tau_y are not only zeroes, this also fails for things like NaN and Inf
        self.assertTrue(numpy.sum(numpy.abs(tau_x)) > (0.0 | units.Pa), msg='Expected tau_x to contain some actual values')
        self.assertTrue(numpy.sum(numpy.abs(tau_y)) > (0.0 | units.Pa), msg='Expected tau_x to contain some actual values')

        #check to see if we can set and retrieve the same wind stress
        size = instance.get_domain_size()
        tau_x = numpy.random.random(size) | units.Pa
        tau_y = numpy.random.random(size) | units.Pa

        #cannot write to forcings directly unfortunately
        #instance.forcings.tau_x = tau_x
        #instance.forcings.tau_y = tau_y
        forcings=instance.forcings.empty_copy()
        forcings.tau_x=tau_x
        forcings.tau_y=tau_y
        forcings.new_channel_to(instance.forcings).copy_attributes(["tau_x","tau_y"])

        self.assertEquals(instance.state_machine._current_state.name, 'EDIT_FORCINGS')

        #now retrieve the wind stress from the model
        tau_x_returned = instance.forcings.tau_x
        tau_y_returned = instance.forcings.tau_y

        self.assertEquals(instance.state_machine._current_state.name, 'RUN')

        #almost equals
        self.assertAlmostRelativeEqual(tau_x, tau_x_returned, places=10)
        self.assertAlmostRelativeEqual(tau_y, tau_y_returned, places=10)

        #check to see if the forcings we have set are not overwritten by internal routines,
        #evolve the model 
        time = instance.get_model_time()
        tend = time + instance.get_timestep_next()
        instance.evolve_model(tend)

        #check if wind stress is still the same
        tau_x_returned = instance.forcings.tau_x
        tau_y_returned = instance.forcings.tau_y
        self.assertAlmostRelativeEqual(tau_x, tau_x_returned, places=10)
        self.assertAlmostRelativeEqual(tau_y, tau_y_returned, places=10)

        instance.stop()



    @nottest
    def test3(self):
        p = POP(**default_options)

        mystr = ''
        mynum = 0
        bogus_file = '/fake/path/to/file'

        mystr = p.get_ts_option()
        self.assertEquals(mystr, 'internal')
        p.set_ts_option('restart')
        mystr = p.get_ts_option()
        self.assertEquals(mystr, 'restart')
        mystr = p.get_ts_file()
        self.assertEquals(mystr, '')
        p.set_ts_file(bogus_file)
        mystr = p.get_ts_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_ts_file_format()
        self.assertEquals(mystr, 'bin')
        p.set_ts_file_format('nc')
        mystr = p.get_ts_file_format()
        self.assertEquals(mystr, 'nc')

        mystr = p.get_distribution()
        self.assertEquals(mystr, 'cartesian')

        p.set_distribution('predefined')
        mystr = p.get_distribution()
        self.assertEquals(mystr, 'predefined')

        mystr = p.get_distribution_file()
        self.assertEquals(mystr, '')
        p.set_distribution_file(bogus_file)
        mystr = p.get_distribution_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_ew_boundary_type()
        self.assertEquals(mystr, 'cyclic')
        p.set_ew_boundary_type('closed')
        mystr = p.get_ew_boundary_type()
        self.assertEquals(mystr, 'closed')
        mystr = p.get_ns_boundary_type()
        self.assertEquals(mystr, 'closed')
        p.set_ns_boundary_type('tripole')
        mystr = p.get_ns_boundary_type()
        self.assertEquals(mystr, 'tripole')

        mystr = p.get_restart_option()
        self.assertEquals(mystr, 'never')
        p.set_restart_option('nday')
        mystr = p.get_restart_option()
        self.assertEquals(mystr, 'nday')
        mynum = p.get_restart_freq_option()
        self.assertEquals(mynum, 1)
        p.set_restart_freq_option(5)
        mynum = p.get_restart_freq_option()
        self.assertEquals(mynum, 5)
        mystr = p.get_restart_file()
        self.assertEquals(mystr, '')
        p.set_restart_file(bogus_file)
        mystr = p.get_restart_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_tavg_option()
        self.assertEquals(mystr, 'never')
        p.set_tavg_option('nday')
        mystr = p.get_tavg_option()
        self.assertEquals(mystr, 'nday')
        mynum = p.get_tavg_freq_option()
        self.assertEquals(mynum, 1)
        p.set_tavg_freq_option(5)
        mynum = p.get_tavg_freq_option()
        self.assertEquals(mynum, 5)
        mystr = p.get_tavg_file()
        self.assertEquals(mystr, '')
        p.set_tavg_file(bogus_file)
        mystr = p.get_tavg_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_movie_option()
        self.assertEquals(mystr, 'never')
        p.set_movie_option('nday')
        mystr = p.get_movie_option()
        self.assertEquals(mystr, 'nday')
        mynum = p.get_movie_freq_option()
        self.assertEquals(mynum, 1)
        p.set_movie_freq_option(5)
        mynum = p.get_movie_freq_option()
        self.assertEquals(mynum, 5)
        mystr = p.get_movie_file()
        self.assertEquals(mystr, '')
        p.set_movie_file(bogus_file)
        mystr = p.get_movie_file()
        self.assertEquals(mystr, bogus_file)

        mystr = p.get_runid()
        self.assertEquals(mystr, 'AMUSE')
        p.set_runid('MYRUNID')
        mystr = p.get_runid()
        self.assertEquals(mystr, 'MYRUNID')

        mystr = p.get_dt_option()
        self.assertEquals(mystr, 'steps_per_day')
        p.set_dt_option('seconds')
        mystr = p.get_dt_option()
        self.assertEquals(mystr, 'seconds')
        p.set_dt_option('steps_per_year')
        mystr = p.get_dt_option()
        self.assertEquals(mystr, 'steps_per_year')

        mynum = p.get_dt_count()
        self.assertEquals(mynum, 45)
        p.set_dt_count(5)
        mynum = p.get_dt_count()
        self.assertEquals(mynum, 5)


    @nottest
    def test4(self):
        p = POP(**default_options)

        mystr = ''
        mynum = 0
        bogus_file = '/fake/path/to/file'

        mystr = p.parameters.ts_option
        self.assertEquals(mystr, 'internal')
        p.parameters.ts_option = 'restart'
        mystr = p.parameters.ts_option
        self.assertEquals(mystr, 'restart')
        mystr = p.parameters.ts_file
        self.assertEquals(mystr, '')
        p.parameters.ts_file = bogus_file
        mystr = p.parameters.ts_file
        self.assertEquals(mystr, bogus_file)

        mystr = p.parameters.ts_file_format
        self.assertEquals(mystr, 'bin')
        p.parameters.ts_file_format = 'nc'
        mystr = p.parameters.ts_file_format
        self.assertEquals(mystr, 'nc')

        mystr = p.parameters.distribution
        self.assertEquals(mystr, 'cartesian')

        p.parameters.distribution = 'predefined'
        mystr = p.parameters.distribution
        self.assertEquals(mystr, 'predefined')

        mystr = p.parameters.distribution_file
        self.assertEquals(mystr, '')
        p.parameters.distribution_file = bogus_file
        mystr = p.parameters.distribution_file
        self.assertEquals(mystr, bogus_file)

        mystr = p.parameters.ew_boundary_type
        self.assertEquals(mystr, 'cyclic')
        p.parameters.ew_boundary_type = 'closed'
        mystr = p.parameters.ew_boundary_type
        self.assertEquals(mystr, 'closed')

        mystr = p.parameters.ns_boundary_type
        self.assertEquals(mystr, 'closed')
        p.parameters.ns_boundary_type = 'tripole'
        mystr = p.parameters.ns_boundary_type
        self.assertEquals(mystr, 'tripole')

        mystr = p.parameters.restart_option
        self.assertEquals(mystr, 'never')
        p.parameters.restart_option = 'nday'
        mystr = p.parameters.restart_option
        self.assertEquals(mystr, 'nday')
        mynum = p.parameters.restart_freq_option
        self.assertEquals(mynum, 1)
        p.parameters.restart_freq_option = 5
        mynum = p.parameters.restart_freq_option
        self.assertEquals(mynum, 5)
        mystr = p.parameters.restart_file
        self.assertEquals(mystr, '')
        p.parameters.restart_file = bogus_file
        mystr = p.parameters.restart_file
        self.assertEquals(mystr, bogus_file)

        mystr = p.parameters.tavg_option
        self.assertEquals(mystr, 'never')
        p.parameters.tavg_option = 'nday'
        mystr = p.parameters.tavg_option
        self.assertEquals(mystr, 'nday')
        mynum = p.parameters.tavg_freq_option
        self.assertEquals(mynum, 1)
        p.parameters.tavg_freq_option = 5
        mynum = p.parameters.tavg_freq_option
        self.assertEquals(mynum, 5)
        mystr = p.parameters.tavg_file
        self.assertEquals(mystr, '')
        p.parameters.tavg_file = bogus_file
        mystr = p.parameters.tavg_file
        self.assertEquals(mystr, bogus_file)

        mystr = p.parameters.movie_option
        self.assertEquals(mystr, 'never')
        p.parameters.movie_option = 'nday'
        mystr = p.parameters.movie_option
        self.assertEquals(mystr, 'nday')
        mynum = p.parameters.movie_freq_option
        self.assertEquals(mynum, 1)
        p.parameters.movie_freq_option = 5
        mynum = p.parameters.movie_freq_option
        self.assertEquals(mynum, 5)
        mystr = p.parameters.movie_file
        self.assertEquals(mystr, '')
        p.parameters.movie_file = bogus_file
        mystr = p.parameters.movie_file
        self.assertEquals(mystr, bogus_file)

        mystr = p.parameters.runid
        self.assertEquals(mystr, 'AMUSE')
        p.parameters.runid = 'MYRUNID'
        mystr = p.parameters.runid
        self.assertEquals(mystr, 'MYRUNID')

        mystr = p.parameters.dt_option
        self.assertEquals(mystr, 'steps_per_day')
        p.parameters.dt_option = 'seconds'
        mystr = p.parameters.dt_option
        self.assertEquals(mystr, 'seconds')
        p.parameters.dt_option = 'steps_per_year'
        mystr = p.parameters.dt_option
        self.assertEquals(mystr, 'steps_per_year')

        mynum = p.parameters.dt_count
        self.assertEquals(mynum, 45)
        p.parameters.dt_count = 5
        mynum = p.parameters.dt_count
        self.assertEquals(mynum, 5)


    @nottest
    def test5(self):
        p = POP(**default_options)

        temp3d = p.elements3d.temp.value_in(units.C)
        self.assertTrue(temp3d.all() < 50., msg='Less than 50 degrees seems reasonable for temperature')
        self.assertTrue(temp3d.all() >-10., msg='More than -10 degrees seems reasonable for temperature')

        salt3d = p.elements3d.salt.value_in(units.g/units.kg)

        self.assertTrue(salt3d.all() < 1., msg='Less than one gram of salt per kg of water seems reasonable for salinity')
        self.assertTrue(salt3d.all() >= 0., msg='Salinity should be positive or 0')

        #check if we can write and again read what was written
        bogus3d = numpy.random.random(p.elements3d.shape) | units.C
        p.elements3d.temp = bogus3d
        temp3d = p.elements3d.temp
        self.assertEquals(temp3d, bogus3d)

        #check if we can write and again read what was written
        bogus3d = numpy.random.random(p.elements3d.shape) | units.g / units.kg
        p.elements3d.salt = bogus3d
        salt3d = p.elements3d.salt
        self.assertEquals(salt3d, bogus3d)

        #check if we can write and again read what was written
        bogus3d = numpy.random.random(p.elements3d.shape) | units.g / units.cm**3
        p.elements3d.rho = bogus3d
        rho3d = p.elements3d.rho
        self.assertEquals(rho3d, bogus3d)

        p.stop()


    @nottest
    def test6(self):
        p = POP(**default_options)

        #check if we can write and again read what was written
        bogus3d = numpy.random.random(p.nodes3d.shape) | units.cm / units.s
        p.nodes3d.xvel = bogus3d
        xvel3d = p.nodes3d.xvel
        self.assertEquals(xvel3d, bogus3d)

        #check if we can write and again read what was written
        bogus3d = numpy.random.random(p.nodes3d.shape) | units.cm / units.s
        p.nodes3d.yvel = bogus3d
        yvel3d = p.nodes3d.yvel
        self.assertEquals(yvel3d, bogus3d)

        p.stop()

    @nottest
    def test7(self):
        p = POP(**default_options)

        #test whether the getters for depth are working as expected
        dzt = p.elements3d.z
        dzu = p.nodes3d.z

        km = p.get_number_of_vertical_levels()

        for i in range(0,km-1):
            self.assertTrue(dzu[5, 5, i] <= dzt[5, 5, i], msg='expect dzu to be equal to or smaller than dzt')
            self.assertTrue(dzt[5, 5, i] > (0.0 | units.cm), msg='expect dzt to be larger than 0.0')
 
        p.stop()

    @nottest
    def test8(self):
        p = POP(**default_options)

        #run for some time to ensure UVEL and VVEL contain something
        time = p.get_model_time()
        tend = time + (1.0 | units.day)
        p.evolve_model(tend)

        #test whether the getter for vertical velocity is working correctly
        xvel = p.nodes3d.xvel
        yvel = p.nodes3d.yvel
        zvel = p.nodes3d.zvel

        km = p.get_number_of_vertical_levels()
        size = p.get_domain_size()
        depth = p.nodes.depth.value_in(units.cm)

        #look for a non-land ocean cell with some depth
        for i in range(0, size[0]-1):
            for j in range(0, size[1]-1):
                if (depth[i,j] > 1000.0):
                    break

        print "Printing info for gridpoint " + str(i) + "," + str(j) + " with depth " + str(depth[i,j])
#        i=1
#        j=127

        print "XVEL:"
        print xvel[i, j, 0:km-1]
        print "YVEL:"
        print yvel[i, j, 0:km-1]
        print "ZVEL:"
        print zvel[i, j, 0:km-1]

        self.assertTrue( any( zv != (0.0 | units.cm / units.s) for zv in zvel[i, j, 0:km-1] ), msg="Expect at least one value to be not equal to zero")

        p.stop()

    @nottest
    def test9(self):
        p = POP(**default_options)

        shf = p.elements.shf

        print "SHF:"
        print shf
        shf_values = shf.value_in(units.W / units.m**2).flatten()
        print shf_values
        self.assertTrue( any( s != 0.0 for s in shf_values), msg="Expect at least one value to be not equal to zero")

        p.stop()
 

    #test the forcing options setting functions
    def test10(self):

        p = POP(**default_options)

        self.assertTrue(p.parameters.windstress_monthly_file == '', msg="Default value for windstress_monthly_file should be empty string")
        self.assertTrue(p.parameters.surface_heat_flux_monthly_file == '', msg="Default value for surface_heat_flux_monthly_file should be empty string")
        self.assertTrue(p.parameters.surface_freshwater_flux_monthly_file == '', msg="Default value for surface_freshwater_flux_monthly_file should be empty string")

        self.assertTrue(p.parameters.windstress_forcing == 'none', msg="Default for windstress forcing should be none")
        self.assertTrue(p.parameters.surface_heat_flux_forcing == 'none', msg="Default for surface heat flux forcing should be none")
        self.assertTrue(p.parameters.surface_freshwater_flux_forcing == 'none', msg="Default for surface freshwater flux forcing should be none")

        bogus_file = '/path/to/non-existent/file'

        p.set_monthly_ws_file(bogus_file)
        p.set_monthly_shf_file(bogus_file)
        p.set_monthly_sfwf_file(bogus_file)

        self.assertTrue(p.parameters.windstress_monthly_file == bogus_file, msg="Error retrieving set value for windstress_monthly_file")
        self.assertTrue(p.parameters.surface_heat_flux_monthly_file == bogus_file, msg="Error retrieving set value for surface_heat_flux_monthly_file")
        self.assertTrue(p.parameters.surface_freshwater_flux_monthly_file == bogus_file, msg="Error retrieving set value for surface_freshwater_flux_monthly_file")

        self.assertTrue(p.parameters.windstress_forcing == 'monthly', msg="Setting for windstress forcing should be monthly")
        self.assertTrue(p.parameters.surface_heat_flux_forcing == 'monthly', msg="Setting for surface heat flux forcing should be monthly")
        self.assertTrue(p.parameters.surface_freshwater_flux_forcing == 'monthly', msg="Setting for surface freshwater flux forcing should be monthly")







