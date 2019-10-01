from __future__ import print_function

import numpy
from amuse.community import *
from amuse.community.interface.common import CommonCode, CommonCodeInterface

from amuse.units import units

import subprocess 

# singly gyre (cf Viebahn 2014)
def single_gyre_wind_model(x,y,Ly,tau):
  tau_x=-tau/numpy.pi*numpy.cos(numpy.pi*y/Ly)
  return tau_x,0.*tau_x

# original qgmodel wind
def jans_wind_model(x,y,Ly,tau):
  tau_x=tau*(numpy.cos(2*numpy.pi*(y/Ly-0.5)) + 
              2*numpy.sin(numpy.pi*(y/Ly-0.5)))
  return tau_x,0.*tau_x

# dijkstra 2005 paper wind
def dijkstra_wind_model(x,y,Ly,tau,sigma=0.):
  tau_x= - tau*( sigma * numpy.cos(numpy.pi*y/Ly)+(1.-sigma)*numpy.cos(2.*numpy.pi*y/Ly) )
  return tau_x,0.*tau_x

class QGmodelInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
        .. [#] Viebahn, J. and Dijkstra, H., International Journal of Bifurcation and Chaos, Vol. 24, No. 2 (2014) 1430007

    """
    
    include_headers = ['worker_code.h']
    use_modules = ['qgmodel',]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="qgmodel_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

#        stacksize=subprocess.check_output('mpirun sh -c "ulimit -s"', shell=True)
#        if int(stacksize) < 2**17:
#          raise Exception("remember to increase the stacksize for qgmodel!")


##    def recommit_parameters(self):
#        pass

    @legacy_function
    def initialize_grid():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function    
    def evolve_model():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='d', direction=function.IN, unit=units.s)
        function.result_type = 'i'
        return function

    def synchronize_model(self):
        pass

    @legacy_function    
    def get_counter():
        function = LegacyFunctionSpecification()  
        function.addParameter('counter', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('tnow', dtype='d', direction=function.OUT, unit=units.s)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('tnow', dtype='d', direction=function.IN, unit=units.s)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_begin_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('begin_time', dtype='d', direction=function.OUT, unit=units.s)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_begin_time():
        function = LegacyFunctionSpecification()  
        function.addParameter('begin_time', dtype='d', direction=function.IN, unit=units.s)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_wind_sigma():
        function = LegacyFunctionSpecification()  
        function.addParameter('wind_sigma', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_wind_sigma():
        function = LegacyFunctionSpecification()  
        function.addParameter('wind_sigma', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_ra_alpha():
        function = LegacyFunctionSpecification()  
        function.addParameter('ra_alpha', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_ra_alpha():
        function = LegacyFunctionSpecification()  
        function.addParameter('ra_alpha', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_raw_alpha():
        function = LegacyFunctionSpecification()  
        function.addParameter('raw_alpha', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_raw_alpha():
        function = LegacyFunctionSpecification()  
        function.addParameter('raw_alpha', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_timestep_method():
        function = LegacyFunctionSpecification()  
        function.addParameter("timestep_method", dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_timestep_method():
        function = LegacyFunctionSpecification()  
        function.addParameter("timestep_method", dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_dpsi_dt():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['dpsi_dt']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.m**2/units.s**2)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_psi1_state():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['psi']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.m**2/units.s)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def set_psi1_state():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['psi']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=units.m**2/units.s)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_wind_field():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['tau_x']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.Pa)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def set_wind_field():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['tau_x']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=units.Pa)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_psi2_state():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['psi_prev']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.m**2/units.s)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def set_psi2_state():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['psi_prev']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=units.m**2/units.s)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_position_of_index():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['x','y']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.m)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_wind_field_position_of_index():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['x','y']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.m)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_index_of_position():
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['x','y']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=units.m)
        for x in ['i','j']:
            function.addParameter(x, dtype='i', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_Lx():
        function = LegacyFunctionSpecification()  
        function.addParameter('Lx', dtype='d', direction=function.OUT, unit=units.m)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_Lx():
        function = LegacyFunctionSpecification()  
        function.addParameter('Lx', dtype='d', direction=function.IN, unit=units.m)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_Ly():
        function = LegacyFunctionSpecification()  
        function.addParameter('Ly', dtype='d', direction=function.OUT, unit=units.m)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_Ly():
        function = LegacyFunctionSpecification()  
        function.addParameter('Ly', dtype='d', direction=function.IN, unit=units.m)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_dy():
        function = LegacyFunctionSpecification()  
        function.addParameter('dy', dtype='d', direction=function.OUT, unit=units.m)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_dy():
        function = LegacyFunctionSpecification()  
        function.addParameter('dy', dtype='d', direction=function.IN, unit=units.m)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_dx():
        function = LegacyFunctionSpecification()  
        function.addParameter('dx', dtype='d', direction=function.OUT, unit=units.m)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_dx():
        function = LegacyFunctionSpecification()  
        function.addParameter('dx', dtype='d', direction=function.IN, unit=units.m)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_dt():
        function = LegacyFunctionSpecification()  
        function.addParameter('dt', dtype='d', direction=function.OUT, unit=units.s)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_dt():
        function = LegacyFunctionSpecification()  
        function.addParameter('dt', dtype='d', direction=function.IN, unit=units.s)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_T():
        function = LegacyFunctionSpecification()  
        function.addParameter('T', dtype='d', direction=function.OUT, unit=units.s)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_T():
        function = LegacyFunctionSpecification()  
        function.addParameter('T', dtype='d', direction=function.IN, unit=units.s)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_tau():
        function = LegacyFunctionSpecification()  
        function.addParameter('tau', dtype='d', direction=function.OUT, unit=units.Pa)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_tau():
        function = LegacyFunctionSpecification()  
        function.addParameter('tau', dtype='d', direction=function.IN, unit=units.Pa)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_err_tol():
        function = LegacyFunctionSpecification()  
        function.addParameter('err_tol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_err_tol():
        function = LegacyFunctionSpecification()  
        function.addParameter('err_tol', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_max_it():
        function = LegacyFunctionSpecification()  
        function.addParameter('max_it', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_max_it():
        function = LegacyFunctionSpecification()  
        function.addParameter('max_it', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_A_H():
        function = LegacyFunctionSpecification()  
        function.addParameter('A_H', dtype='d', direction=function.OUT,unit=units.m**2/units.s)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_A_H():
        function = LegacyFunctionSpecification()  
        function.addParameter('A_H', dtype='d', direction=function.IN,unit=units.m**2/units.s)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_R_H():
        function = LegacyFunctionSpecification()  
        function.addParameter('R_H', dtype='d', direction=function.OUT, unit=units.s**-1)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_R_H():
        function = LegacyFunctionSpecification()  
        function.addParameter('R_H', dtype='d', direction=function.IN, unit=units.s**-1)
        function.result_type = 'i'
        return function
    
    @legacy_function    
    def get_lambda0():
        function = LegacyFunctionSpecification()  
        function.addParameter('lambda0', dtype='d', direction=function.OUT, unit=units.m**-1)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_lambda0():
        function = LegacyFunctionSpecification()  
        function.addParameter('lambda0', dtype='d', direction=function.IN, unit=units.m**-1)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_lambda1():
        function = LegacyFunctionSpecification()  
        function.addParameter('lambda1', dtype='d', direction=function.OUT, unit=units.m**-1)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_lambda1():
        function = LegacyFunctionSpecification()  
        function.addParameter('lambda1', dtype='d', direction=function.IN, unit=units.m**-1)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_e111():
        function = LegacyFunctionSpecification()  
        function.addParameter('e111', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_e111():
        function = LegacyFunctionSpecification()  
        function.addParameter('e111', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_phi1z0():
        function = LegacyFunctionSpecification()  
        function.addParameter('phi1z0', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_phi1z0():
        function = LegacyFunctionSpecification()  
        function.addParameter('phi1z0', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_H():
        function = LegacyFunctionSpecification()  
        function.addParameter('H', dtype='d', direction=function.OUT, unit=units.m)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_H():
        function = LegacyFunctionSpecification()  
        function.addParameter('H', dtype='d', direction=function.IN, unit=units.m)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_rho():
        function = LegacyFunctionSpecification()  
        function.addParameter('rho', dtype='d', direction=function.OUT, unit=(units.kg/units.m**3))
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_rho():
        function = LegacyFunctionSpecification()  
        function.addParameter('rho', dtype='d', direction=function.IN, unit=(units.kg/units.m**3))
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_savecounter():
        function = LegacyFunctionSpecification()  
        function.addParameter('savecounter', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_savecounter():
        function = LegacyFunctionSpecification()  
        function.addParameter('savecounter', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_beta0():
        function = LegacyFunctionSpecification()  
        function.addParameter('beta0', dtype='d', direction=function.OUT, unit=(units.m*units.s)**-1)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_beta0():
        function = LegacyFunctionSpecification()  
        function.addParameter('beta0', dtype='d', direction=function.IN, unit=(units.m*units.s)**-1)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_relax_coef():
        function = LegacyFunctionSpecification()  
        function.addParameter('relax_coef', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_relax_coef():
        function = LegacyFunctionSpecification()  
        function.addParameter('relax_coef', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_Nx():
        function = LegacyFunctionSpecification()  
        function.addParameter('Nx', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_Ny():
        function = LegacyFunctionSpecification()  
        function.addParameter('Ny', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_Nt():
        function = LegacyFunctionSpecification()  
        function.addParameter('Nt', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_Nm():
        function = LegacyFunctionSpecification()  
        function.addParameter('Nm', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_Nm():
        function = LegacyFunctionSpecification()  
        function.addParameter('Nm', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    @legacy_function    
    def get_interface_wind():
        function = LegacyFunctionSpecification()  
        function.addParameter('interface_wind', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_interface_wind():
        function = LegacyFunctionSpecification()  
        function.addParameter('interface_wind', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function


    @legacy_function
    def get_boundary_index_range_inclusive():
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN)
        function.addParameter('minx', dtype='i', direction=function.OUT)
        function.addParameter('maxx', dtype='i', direction=function.OUT)
        function.addParameter('miny', dtype='i', direction=function.OUT)
        function.addParameter('maxy', dtype='i', direction=function.OUT)
        function.addParameter('minz', dtype='i', direction=function.OUT)
        function.addParameter('maxz', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_index_range_inclusive():
        function = LegacyFunctionSpecification()
        function.addParameter('minx', dtype='i', direction=function.OUT)
        function.addParameter('maxx', dtype='i', direction=function.OUT)
        function.addParameter('miny', dtype='i', direction=function.OUT)
        function.addParameter('maxy', dtype='i', direction=function.OUT)
        function.addParameter('minz', dtype='i', direction=function.OUT)
        function.addParameter('maxz', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_wind_field_index_range_inclusive():
        function = LegacyFunctionSpecification()
        function.addParameter('minx', dtype='i', direction=function.OUT)
        function.addParameter('maxx', dtype='i', direction=function.OUT)
        function.addParameter('miny', dtype='i', direction=function.OUT)
        function.addParameter('maxy', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

# new boundary stuff
    @legacy_function    
    def set_boundary_conditions():
        function = LegacyFunctionSpecification()  
        for x in ["boundary_west","boundary_east","boundary_south","boundary_north"]:
            function.addParameter(x, dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_boundary_conditions():
        function = LegacyFunctionSpecification()  
        for x in ["boundary_west","boundary_east","boundary_south","boundary_north"]:
            function.addParameter(x, dtype='string', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_boundary_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter("psi", dtype='d', direction=function.IN, unit=units.m**2/units.s)
        function.addParameter("dpsi_dt", dtype='d', direction=function.IN, unit=units.m**2/units.s**2)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_boundary_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN)
        function.addParameter("psi", dtype='d', direction=function.OUT, unit=units.m**2/units.s)
        function.addParameter("dpsi_dt", dtype='d', direction=function.OUT, unit=units.m**2/units.s**2)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'i'
        return function
        
    @legacy_function    
    def get_boundary_position_of_index():
        """
        Retrieves the x and y position of the center of
        the cell with coordinates i, j, k 
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_boundary', dtype='i', direction=function.IN)
        for x in ['x','y']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=units.m)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_psi_state_at_point():
        """
        Retrieves the psi and dpsi_dt for grid k at point x,y
        """
        function = LegacyFunctionSpecification()  
        function.must_handle_array = True
        function.addParameter("cell_size", dtype='d', direction=function.IN, unit=units.m)
        function.addParameter("x", dtype='d', direction=function.IN, unit=units.m)
        function.addParameter("y", dtype='d', direction=function.IN, unit=units.m)
        function.addParameter("k", dtype='i', direction=function.IN,default=1)
        function.addParameter("psi", dtype='d', direction=function.OUT, unit=units.m**2/units.s)
        function.addParameter("dpsi_dt", dtype='d', direction=function.OUT, unit=units.m**2/units.s**2)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'i'
        return function

class QGmodel(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  QGmodelInterface(**options), **options)
        self._position=[0,0] | units.m
    
    def get_position_x(self):
        return self._position[0]
    def get_position_y(self):
        return self._position[1]
    def set_position_x(self,x):
        self._position[0]=x
    def set_position_y(self,y):
        self._position[1]=y
    def get_index_of_position(self,x,y):
        return self.overridden().get_index_of_position(x-self._position[0],y-self._position[1])
    def get_position_of_index(self,i,j,k):
        x,y=self.overridden().get_position_of_index(i,j,k)
        return x+self._position[0],y+self._position[1]
    def get_wind_field_position_of_index(self,i,j):
        x,y=self.overridden().get_wind_field_position_of_index(i,j)
        return x+self._position[0],y+self._position[1]
    def get_psi_state_at_point(self,d,x,y,k=None):
        if k is None:
          return self.overridden().get_psi_state_at_point(d,x-self._position[0],y-self._position[1])
        else:
          return self.overridden().get_psi_state_at_point(d,x-self._position[0],y-self._position[1],k)
    def get_boundary_position_of_index(self,i,j,k,index_of_boundary):
        x,y=self.overridden().get_boundary_position_of_index(i,j,k,index_of_boundary)
        return x+self._position[0],y+self._position[1]

    @property
    def forcings(self):
        return self.wind_field
    
    def define_grids(self, object):
        object.define_grid('grid',axes_names = ['x','y'],grid_class=datamodel.CartesianGrid)
        object.set_grid_range('grid', 'get_index_range_inclusive')    
        object.add_getter('grid', 'get_dpsi_dt', names=('dpsi_dt',))
        object.add_getter('grid', 'get_psi1_state', names=('psi',))
        object.add_setter('grid', 'set_psi1_state', names=('psi',))
        object.add_getter('grid', 'get_psi2_state', names=('psi_prev',))
        object.add_setter('grid', 'set_psi2_state', names=('psi_prev',))
        object.add_getter('grid', 'get_position_of_index', names=('x','y'))

        object.define_grid('wind_field',axes_names = ['x','y'],grid_class=datamodel.CartesianGrid)
        object.set_grid_range('wind_field', 'get_wind_field_index_range_inclusive')    
        object.add_getter('wind_field', 'get_wind_field', names=('tau_x',))
        object.add_setter('wind_field', 'set_wind_field', names=('tau_x',))
        object.add_getter('wind_field', 'get_wind_field_position_of_index', names=('x','y'))

        self._boundaries={}
        for i,side in enumerate(["west","east","south","north"]):
          name="boundary_"+side
          object.define_grid(name,axes_names = ['x','y'],grid_class=datamodel.CartesianGrid)
          object.set_grid_range(name, 'get_boundary_index_range_inclusive')    
          object.add_getter(name, 'get_boundary_state', names=('psi','dpsi_dt'))
          object.add_setter(name, 'set_boundary_state', names=('psi','dpsi_dt'))
          object.add_getter(name, 'get_boundary_position_of_index', names=('x','y'))
          object.define_extra_keywords(name, {'index_of_boundary': i+1})
          self._boundaries[i+1]="self."+name
          self._boundaries[side]="self."+name

    def west_boundary_updater(self,boundary):
        raise Exception("define update function for interface boundaries")
    def east_boundary_updater(self,boundary):
        raise Exception("define update function for interface boundaries")
    def south_boundary_updater(self,boundary):
        raise Exception("define update function for interface boundaries")
    def north_boundary_updater(self,boundary):
        raise Exception("define update function for interface boundaries")
        
    def update_boundaries(self):
        if self.parameters.boundary_west=="interface": self.west_boundary_updater(self.boundary_west)
        if self.parameters.boundary_east=="interface": self.east_boundary_updater(self.boundary_east)
        if self.parameters.boundary_south=="interface": self.south_boundary_updater(self.boundary_south)
        if self.parameters.boundary_north=="interface": self.north_boundary_updater(self.boundary_north)

    def boundaries(self,x):
        return eval(self._boundaries[x])

    def set_wind(self,f): # convenience function to set wind and turn on interface_wind
        if not self.parameters.interface_wind:
          raise Exception("unexpected set_wind")
        tau_x,tau_y=f(self.wind_field.x,self.wind_field.y)
        self.wind_field.tau_x=tau_x

    def define_state(self, object): 
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('!UNINITIALIZED!INITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_transition('END', 'INITIALIZED', 'initialize_code',False)        
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')        
        
# may only be set before Initialization
        object.add_method('INITIALIZED', 'set_Lx')
        object.add_method('INITIALIZED', 'set_Ly')
        object.add_method('INITIALIZED', 'set_dx')
        object.add_method('INITIALIZED', 'set_dy')
        object.add_method('INITIALIZED', 'set_Nm')

# Nx and Ny need commit_parameter, all need initialize_code
        object.add_method('!UNINITIALIZED', 'before_get_parameter')
#        object.add_method('!UNINITIALIZED', "get_Nx")
#        object.add_method('!UNINITIALIZED', "get_Ny")

        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('EDIT', 'RUN', 'initialize_grid')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_transition('EVOLVED','RUN', 'synchronize_model')
        object.add_transition('EVOLVED','EDIT', 'synchronize_model')
        object.add_transition('RUN','EDIT', 'synchronize_model')
        object.add_transition('INITIALIZED','INITIALIZED','before_set_parameter', False)
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter')
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        
        
        object.add_method('INITIALIZED', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
 
        object.add_method('RUN', 'get_dpsi_dt')
        object.add_method("RUN", 'get_psi_state_at_point')
        object.add_method("RUN", 'get_wind_field')

        object.add_method('EDIT', 'set_psi1_state')
        object.add_method('EDIT', 'set_psi2_state')
        object.add_method('EDIT', 'set_wind_field')
        for state in ["RUN","EDIT"]:
          object.add_method(state, 'get_index_range_inclusive')
          object.add_method(state, 'get_wind_field_index_range_inclusive')
          object.add_method(state, 'get_boundary_index_range_inclusive')
          object.add_method(state, 'get_psi1_state')
          object.add_method(state, 'get_psi2_state')
          object.add_method(state, 'get_boundary_state')
          object.add_method(state, 'set_boundary_state')

        
    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")

    def commit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()

    def recommit_parameters(self):
        self.parameters.send_cached_parameters_to_code()
        self.overridden().recommit_parameters()

    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_position_x", 
            "set_position_x",
            "position_x", 
            "position of the grid in x direction", 
            default_value = 0. | units.m
        )

        object.add_method_parameter(
            "get_position_y", 
            "set_position_y",
            "position_y", 
            "position of the grid in y direction", 
            default_value = 0. | units.m
        )

        object.add_method_parameter(
            "get_Lx", 
            "set_Lx",
            "Lx", 
            "Lx is the x size of the grid", 
            default_value = 4.e6 | units.m
        )
        object.add_method_parameter(
            "get_begin_time", 
            "set_begin_time",
            "begin_time", 
            "time of starting time step", 
            default_value = 0.0 | units.s
        )

        object.add_method_parameter(
            "get_wind_sigma", 
            "set_wind_sigma",
            "wind_sigma", 
            "wind mix parameter (-1 - 1, <-1: jans wind)", 
            default_value = -99.0
        )

        object.add_method_parameter(
            "get_ra_alpha", 
            "set_ra_alpha",
            "robert_asselin_alpha", 
            "robert_asselin filter parameter (0.01-0.1)", 
            default_value = 0.1
        )

        object.add_method_parameter(
            "get_raw_alpha", 
            "set_raw_alpha",
            "raw_parameter", 
            "robert-asselin-williams filter parameter (0.0-1., 1.= normal ra leapfrog)", 
            default_value = 1.
        )

        object.add_method_parameter(
            "get_timestep_method", 
            "set_timestep_method",
            "timestep_method", 
            "timestep method (euler, leapfrog)", 
            default_value = "leapfrog"
        )

        object.add_boolean_parameter(
            "get_interface_wind", 
            "set_interface_wind",
            "interface_wind", 
            "set wind field thorugh interface (True) or let QGmodel generate wind (False)", 
            default_value = False
        )

        object.add_method_parameter(
            "get_Ly", 
            "set_Ly",
            "Ly", 
            "Ly is the y size of the grid", 
            default_value = 4.e6 | units.m
        )
        object.add_method_parameter(
            "get_dx", 
            "set_dx",
            "dx", 
            "dx is the mesh size (x)", 
            default_value = 1.e4 | units.m
        )
        object.add_method_parameter(
            "get_dy", 
            "set_dy",
            "dy", 
            "dy is the mesh size (y)", 
            default_value = 1.e4 | units.m
        )

        object.add_method_parameter(
            "get_Nx", 
            None,
            "Nx", 
            "Nx is the number of mesh points in the x (east-west) direction", 
            default_value = 401
        )
        object.add_method_parameter(
            "get_Ny", 
            None,
            "Ny", 
            "Nx is the number of mesh points in the y (north-south) direction", 
            default_value = 401
        )                
        
        object.add_method_parameter(
            "get_dt", 
            "set_dt",
            "dt", 
            "timestep of the evolve", 
            default_value = 3600.0 | units.s
        ) 
        object.add_method_parameter(
            "get_A_H", 
            "set_A_H",
            "A_H", 
            "turbulent lateral friction coefficient",
            default_value = 100.0 | units.m**2/units.s
        ) 
        object.add_method_parameter(
            "get_R_H", 
            "set_R_H",
            "R_H", 
            "bottom friction coefficient (f_0*delta_E/2/H)", 
            default_value = 0.0 | units.s**-1
        ) 
        object.add_method_parameter(
            "get_T", 
            "set_T",
            "T", 
            "T output parameter", # time output (not used?)
            default_value = 86400.0 | units.s
        ) 
        object.add_method_parameter(
            "get_savecounter", 
            "set_savecounter",
            "savecounter", 
            "savecounter", #(?) not used?
            default_value = 48
        ) 
        object.add_method_parameter(
            "get_lambda0", 
            "set_lambda0",
            "lambda0", 
            "rossby wavelength 0 (1/rossby radius)", 
            default_value = 0.0 | units.m**-1
        ) 
        object.add_method_parameter(
            "get_lambda1", 
            "set_lambda1",
            "lambda1", 
            "rossby wavelength 1 (1/rossby radius)", 
            default_value = 2.e-5 | units.m**-1
        )
        object.add_method_parameter(
            "get_e111", 
            "set_e111",
            "e111", 
            "e111 interaction parameter of diff modes", 
            default_value = 0.0
        )  
        object.add_method_parameter(
            "get_phi1z0", 
            "set_phi1z0",
            "phi1z0", 
            "phi1(z0) value of the first mode at the surface, sqrt(2) for constant N", 
            default_value = 1.4142135623731
        ) 
        object.add_method_parameter(
            "get_H", 
            "set_H",
            "ocean_depth", 
            "ocean depth", #(check) 
            default_value = 4000.0 | units.m
        )
        object.add_method_parameter(
            "get_rho", 
            "set_rho",
            "rho", 
            "density ", 
            default_value = 1000.0 | (units.kg/units.m**3)
        )
        object.add_method_parameter(
           "get_beta0", 
            "set_beta0",
            "beta0", 
            "beta0 coriolis coefficient", 
            default_value = 1.8616e-11 | (units.m*units.s)**-1
        )
        object.add_method_parameter(
            "get_err_tol", 
            "set_err_tol",
            "err_tol", 
            "error tolerance", 
            default_value = 1.e-6
        ) 
        object.add_method_parameter(
            "get_tau", 
            "set_tau",
            "tau", 
            "tau wind forcing parameter ", # (?)
            default_value = 0.05 | units.Pa 
        ) 
        object.add_method_parameter(
            "get_max_it", 
            "set_max_it",
            "max_it", 
            "maximum iterations", #(check) 
            default_value = 200
        ) 
        object.add_method_parameter(
            "get_relax_coef", 
            "set_relax_coef",
            "relax_coef", 
            "relaxation coefficient", # (meaning?)
            default_value = 1.7
        )    
        object.add_method_parameter(
            "get_Nm", 
            "set_Nm",
            "Nm", 
            "Number of depth levels ( base functions)", 
            default_value = 1 
        ) 

        object.add_caching_parameter(
            "set_boundary_conditions", 
            "boundary_west",
            "boundary_west", 
            "boundary conditions on first (west) X boundary", 
            "free_slip",
        )
        
        object.add_caching_parameter(
            "set_boundary_conditions", 
            "boundary_east",
            "boundary_east", 
            "boundary conditions on second (east) X boundary",
            "free_slip",
        )
        
        object.add_caching_parameter(
            "set_boundary_conditions", 
            "boundary_south",
            "boundary_south", 
            "boundary conditions on first (south) Y boundary", 
            "free_slip",
        )
        
        
        object.add_caching_parameter(
            "set_boundary_conditions", 
            "boundary_north",
            "boundary_north", 
            "boundary conditions on second (north) Y boundary",
            "free_slip",
        )
                
        object.add_vector_parameter(
            "x_boundary_conditions",
            "boundary conditions for the zonal (x) directorion",
            ("boundary_west", "boundary_east")
        )
                
        object.add_vector_parameter(
            "y_boundary_conditions",
            "boundary conditions for the meridional (y) directorion",
            ("boundary_south", "boundary_north")
        )


class QGmodelWithRefinements(QGmodel):
    def __init__(self,*args,**kwargs):
      self.refinements=[]
      self.verbose=False
      if 'verbose' in kwargs:
        self.verbose=kwargs["verbose"]
      self._args=args
      self._kwargs=kwargs
      QGmodel.__init__(self,*args,**kwargs)
          
    def interpolate_grid(self,grid):
      copy=grid.empty_copy()    
      x=grid.x.flatten()
      y=grid.y.flatten()
      dx=max(self.grid.cellsize()[0],grid.cellsize()[0])
      psi,dpsi_dt=self.get_psi_state_at_point(dx,x,y)
      copy.psi=psi.reshape(copy.shape)
      copy.dpsi_dt=dpsi_dt.reshape(copy.shape)    
      channel=copy.new_channel_to(grid)
      channel.copy_attributes(["psi","dpsi_dt"])
  
    def add_refinement(self,sys=None,position=[0.,0.] | units.m,offset=None,parameters=dict()):
      if sys is None:
        sys=self.__class__(*self._args,**self._kwargs)
      if callable(sys):
        sys=sys()
      sys.parameters.reset_from_memento(self.parameters)
      for param,val in parameters.items():
        setattr(sys.parameters, param, val)
      if offset is not None:
        position=[self.parameters.position_x+offset[0],
                  self.parameters.position_y+offset[1] ]
      sys.parameters.position_x=position[0]
      sys.parameters.position_y=position[1]

      self.refinements.append(sys)
      
      minpos=sys.grid.get_minimum_position()
      maxpos=sys.grid.get_maximum_position()
      
      parent_minpos=self.grid.get_minimum_position()
      parent_maxpos=self.grid.get_maximum_position()
      dx,dy=self.grid.cellsize()
      self._gridx=self.grid.x
      self._gridy=self.grid.y
      
      if minpos[0]<parent_minpos[0]+dx/2:
        print("linking west boundary")
        sys.parameters.boundary_west=self.parameters.boundary_west
        sys.west_boundary_updater=self.west_boundary_updater
      else:
        print("west boundary interpolated")
        sys.parameters.boundary_west="interface"
        sys.west_boundary_updater=self.interpolate_grid
      
      if maxpos[0]>parent_maxpos[0]-dx/2:
        print("linking east boundary")
        sys.parameters.boundary_east=self.parameters.boundary_east
        sys.east_boundary_updater=self.east_boundary_updater
      else:
        print("east boundary interpolated")
        sys.parameters.boundary_east="interface"
        sys.east_boundary_updater=self.interpolate_grid
  
      if minpos[1]<parent_minpos[1]+dy/2:
        print("linking south boundary")
        sys.parameters.boundary_south=self.parameters.boundary_south
        sys.south_boundary_updater=self.south_boundary_updater
      else:
        print("south boundary interpolated")
        sys.parameters.boundary_south="interface"
        sys.south_boundary_updater=self.interpolate_grid
        
      if maxpos[1]>parent_maxpos[1]-dy/2:
        print("linking north boundary")
        sys.parameters.boundary_north=self.parameters.boundary_north
        sys.north_boundary_updater=self.north_boundary_updater
      else:
        print("north boundary interpolated")
        sys.parameters.boundary_north="interface"
        sys.north_boundary_updater=self.interpolate_grid
      if not hasattr(sys,"get_psi_dpsidt"):
        sys.get_psi_dpsidt=sys.get_psi_state_at_point
  
      return sys

    def evolve_model(self,tend,dt=None,method="EBER",restriction_correction=None):
      if method=="EBER":
        self.evolve_model_EBER(tend,dt,restriction_correction=restriction_correction)
      elif method=="EREB":
        self.evolve_model_EREB(tend,dt)
      elif method=="BRE":
        self.evolve_model_BRE(tend,dt)
      elif method=="RBE":
        self.evolve_model_RBE(tend,dt)
      else:
        raise Exception("unknown method")
  
    def evolve_model_EBER(self,tend,dt=None,restriction_correction=None):
      if dt is None:
        dt=2*self.parameters.dt
      tnow=self.model_time
      while tnow<tend-dt/2:
        self.overridden().evolve_model(tnow+dt/2)
        for sys in self.refinements:
          if self.verbose:  print("update boundaries...", end=' ')
          sys.update_boundaries()
          if self.verbose:  print("done")
          sys.evolve_model(tnow+dt)
        self.overridden().evolve_model(tnow+dt)
        if len(self.refinements): 
          if self.verbose: print("update refined regions...", end=' ')
          self.update_refined_regions()
          if restriction_correction is not None:
            restriction_correction(self.grid)
          if self.verbose: print("done")
        tnow=self.model_time

    def evolve_model_EREB(self,tend,dt=None):
      if dt is None:
        dt=2*self.parameters.dt
        
      first=True
      tnow=self.model_time
      while tnow<tend-dt/2:
        if self.verbose:  print("update boundaries...", end=' ')
        if first:
          for sys in self.refinements:
            if self.verbose:  print("update boundaries...", end=' ')
            sys.update_boundaries()
          if self.verbose:  print("done")
          first=False
        for sys in self.refinements:
          sys.evolve_model(tnow+dt/2)
        self.overridden().evolve_model(tnow+dt/2)

        if self.verbose: print("update refined regions...", end=' ')
        if len(self.refinements): 
          self.update_refined_regions()
        if self.verbose: print("done")

        self.overridden().evolve_model(tnow+dt)
        if self.verbose:  print("update boundaries...", end=' ')
        for sys in self.refinements:
          if self.verbose:  print("update boundaries...", end=' ')
          sys.update_boundaries()
        for sys in self.refinements:
          sys.evolve_model(tnow+dt)

        tnow=self.model_time

    def evolve_model_BRE(self,tend,dt=None):
      if dt is None:
        dt=2*self.parameters.dt
      tnow=self.model_time
      while tnow<tend-dt/2:
        for sys in self.refinements:
          if self.verbose:  print("update boundaries...", end=' ')
          sys.update_boundaries()
          if self.verbose:  print("done")
        if len(self.refinements): 
          if self.verbose: print("update refined regions...", end=' ')
          self.update_refined_regions()
          if self.verbose: print("done")
        self.overridden().evolve_model(tnow+dt)
        for sys in self.refinements:
          sys.evolve_model(tnow+dt)
        tnow=self.model_time

    def evolve_model_RBE(self,tend,dt=None):
      if dt is None:
        dt=2*self.parameters.dt
      tnow=self.model_time
      while tnow<tend-dt/2:
        if len(self.refinements): 
          if self.verbose: print("update refined regions...", end=' ')
          self.update_refined_regions()
          if self.verbose: print("done")
        for sys in self.refinements:
          if self.verbose:  print("update boundaries...", end=' ')
          sys.update_boundaries()
          if self.verbose:  print("done")
        self.overridden().evolve_model(tnow+dt)
        for sys in self.refinements:
          sys.evolve_model(tnow+dt)
        tnow=self.model_time

  
    def get_psi_dpsidt(self,dx,x,y,k=None):
      minx=x-dx/2
      maxx=x+dx/2
      miny=y-dx/2
      maxy=y+dx/2
      done=numpy.zeros(minx.shape,dtype=bool)
      psi=numpy.zeros(minx.shape) | units.m**2/units.s
      dpsi=numpy.zeros(minx.shape) | units.m**2/units.s**2
      for sys in self.refinements+[self]:
        gridminx,gridminy=sys.grid.get_minimum_position()
        gridmaxx,gridmaxy=sys.grid.get_maximum_position()
        select=- ( (gridminx>minx)+(gridmaxx<maxx)+(gridminy>miny)+(gridmaxy<maxy)+done )
        if numpy.any(select) and sys is not self: 
          p,d=sys.get_psi_dpsidt(dx[select],x[select],y[select])
          psi[select]=p
          dpsi[select]=d
        done=done+select
      if numpy.any(select):
        p,d=self.get_psi_state_at_point(dx[select],x[select],y[select])
        psi[select]=p
        dpsi[select]=d
      return psi,dpsi
  
    def update_refined_regions(self):
      x=self._gridx
      y=self._gridy
      dx=0.*x+self.grid.cellsize()[0]

      minx=x-dx/2
      maxx=x+dx/2
      miny=y-dx/2
      maxy=y+dx/2
      psi=numpy.zeros(minx.shape) | units.m**2/units.s
      done=numpy.zeros(minx.shape,dtype=bool)
      for sys in self.refinements:
        gridminx,gridminy=sys.grid.get_minimum_position()
        gridmaxx,gridmaxy=sys.grid.get_maximum_position()
#        gridminx+=sys.grid.cellsize()[0]  # interior only...
#        gridmaxx-=sys.grid.cellsize()[0]
#        gridminy+=sys.grid.cellsize()[1]
#        gridmaxy-=sys.grid.cellsize()[1]                
#        select=-( (gridminx>minx)+(gridmaxx<maxx)+(gridminy>miny)+(gridmaxy<maxy)+done)
        select=(~done)*(minx>=gridminx)*(maxx<=gridmaxx)*(miny>=gridminy)*(maxy<=gridmaxy) 
        if numpy.any(select): 
          p,d=sys.get_psi_state_at_point(dx[select],x[select],y[select])
          psi[select]=p
        done=done+select
      if numpy.any(done):
        self.grid[done].psi=psi[done]

class QGmodel_(QGmodel):
    def evolve_model(self,tend,dt=None):
      if dt is None:
        dt=self.parameters.dt
      tnow=self.model_time
      while self.model_time<tend-dt/2:
        self.initialize_grid()
        self.overridden().evolve_model(tnow+dt)
        tnow=self.model_time
