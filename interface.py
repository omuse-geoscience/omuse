from amuse.community import *
from amuse.community.interface.common import CommonCode, CommonCodeInterface

from amuse.units import units

class QGmodelInterface(CodeInterface, CommonCodeInterface,LiteratureReferencesMixIn):
    """
        .. [#] Viebahn, J. and Dijkstra, H., International Journal of Bifurcation and Chaos, Vol. 24, No. 2 (2014) 1430007

    """
    
    include_headers = ['worker_code.h']
    use_modules = ['qgmodel',]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="qgmodel_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

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
    def get_free_slip():
        function = LegacyFunctionSpecification()  
        function.addParameter('free_slip', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function    
    def set_free_slip():
        function = LegacyFunctionSpecification()  
        function.addParameter('free_slip', dtype='i', direction=function.IN)
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
        
    def get_index_range_inclusive(self):
        return 1,self.get_Nx()['Nx'],1,self.get_Ny()['Ny'],1,self.get_Nm()['Nm']

class QGmodel(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  QGmodelInterface(**options), **options)
    
    def define_particle_sets(self, object):
        object.define_grid('grid')
        object.set_grid_range('grid', 'get_index_range_inclusive')    
        object.add_getter('grid', 'get_psi1_state', names=('psi',))
        object.add_setter('grid', 'set_psi1_state', names=('psi',))
        object.add_getter('grid', 'get_psi2_state', names=('psi_prev',))
        object.add_setter('grid', 'set_psi2_state', names=('psi_prev',))
        object.add_getter('grid', 'get_position_of_index', names=('x','y'))


    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")
    
    def define_parameters(self, object):
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
            "wind mix parameter (0-1)", 
            default_value = 1.0
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
            "get_dt", 
            "set_dt",
            "dt", 
            "timestep of the evolve", 
            default_value = 1800.0 | units.s
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
            "bottom friction coefficient", 
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
        object.add_boolean_parameter(
            "get_free_slip", 
            "set_free_slip",
            "free_slip", 
            "whether to use free_slip boundary conditions, True=free slip, False =no slip", 
            default_value = True
        ) 
        object.add_method_parameter(
            "get_Nm", 
            "set_Nm",
            "Nm", 
            "Number of depth levels", 
            default_value = 1 
        ) 
