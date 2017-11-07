import os.path
from omuse.units import units
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system
from amuse.community.interface.stopping_conditions import StoppingConditionInterface, StoppingConditions
from amuse import datamodel

from amuse.units import trigo

import numpy

class DalesInterface(CodeInterface,
                     CommonCodeInterface,
                     StoppingConditionInterface,
                     LiteratureReferencesMixIn):
    """

    DALES - Dutch Atmospheric Large Eddy Simulation

    """
    use_modules=["StoppingConditions","dales_interface"]

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)
        
    def name_of_the_worker(self):
        return "dales_worker"

    @remote_function
    def get_input_file():
        returns (input_file="s")

    @remote_function
    def set_input_file(input_file="namoptions"):
        pass

    @remote_function
    def set_workdir(directory="./"):
        pass
    
    @remote_function
    def get_model_time():
        returns (time=0.| units.s)

    @remote_function
    def get_timestep():
        returns (dt=0. | units.s)

    @remote_function
    def set_surface_pressure(p = 1.0e5 | units.Pa):
        pass

    @remote_function
    def set_tendency_surface_pressure(p = 1.0e5 | units.Pa):
        pass
    
    @remote_function
    def get_surface_pressure():
        returns (p=0.| units.Pa)
    
    @remote_function
    def commit_grid():        
        pass

# getter functions for vertical profiles - slab averages
# these take a dummy array as input, and return output of the same length
 

    @remote_function(must_handle_array=True)
    def get_profile_U_(k=0):
        returns (out=0.| units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_V_(k=0):
        returns (out=0.| units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_profile_W_(k=0):
        returns (out=0.| units.m / units.s)
        
    @remote_function(must_handle_array=True)
    def get_profile_THL_(k=0):
        returns (out=0.| units.K)
        
    @remote_function(must_handle_array=True)
    def get_profile_QT_(k=0):
        returns (out=0.)

    @remote_function(must_handle_array=True)
    def get_profile_QL_(k=0):
        returns (out=0.)

    @remote_function(must_handle_array=True)
    def get_profile_QL_ice_(k=0):
        returns (out=0.)

    @remote_function(must_handle_array=True)
    def get_profile_QR_(k=0):
        returns (out=0.)
        
    @remote_function(must_handle_array=True)
    def get_profile_E12_(k=0):
        returns (out=0.)

    @remote_function(must_handle_array=True)
    def get_profile_T_(k=0):
        returns (out=0. | units.K)

    # getter for cloud fraction. Uses the index array to define slabs.
    @remote_function(must_handle_array=True)
    def get_cloudfraction(k=0):
        returns (out=0.)


# getter functions for height levels
# these take a dummy array as input, and return output of the same length
    @remote_function(must_handle_array=True)
    def get_zf_(k=0):
        returns (out=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_zh_(k=0):
        returns (out=0. | units.m)

    @remote_function(must_handle_array=True)
    def get_presf_(k=0):
        returns (out=0. | units.Pa)

    @remote_function(must_handle_array=True)
    def get_presh_(k=0):
        returns (out=0. | units.Pa)


        
# setter functions for vertical tendencies / forcings
    @remote_function(must_handle_array=True)
    def set_tendency_U(a=0.):
        returns () 

    @remote_function(must_handle_array=True)
    def set_tendency_V(a=0.):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_THL(a=0.):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_QT(a=0.):
        returns ()

# getter functions for 3D fields usning index arrays
    @remote_function(must_handle_array=True)
    def get_field_U(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_V(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_W(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.m / units.s)

    @remote_function(must_handle_array=True)
    def get_field_THL(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.K)

    @remote_function(must_handle_array=True)
    def get_field_QT(g_i=0,g_j=0,g_k=0):
        returns (a=0.)

    @remote_function(must_handle_array=True)
    def get_field_QL(g_i=0,g_j=0,g_k=0):
        returns (a=0.)

    @remote_function(must_handle_array=True)
    def get_field_E12(g_i=0,g_j=0,g_k=0):
        returns (a=0.)

    @remote_function(must_handle_array=True)
    def get_field_T(g_i=0,g_j=0,g_k=0):
        returns (a=0.| units.K)

    @remote_function(must_handle_array=True)
    def get_field_LWP(g_i=0,g_j=0):
        returns (a=0.)

    @remote_function(must_handle_array=True)
    def get_field_RWP(g_i=0,g_j=0):
        returns (a=0.) 

    @remote_function(must_handle_array=True)
    def get_field_TWP(g_i=0,g_j=0):
        returns (a=0.)
        
    # setter functions for 3D fields usning index arrays
    @remote_function(must_handle_array=True)
    def set_field_U(g_i=0,g_j=0,g_k=0,a=0.):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_V(g_i=0,g_j=0,g_k=0,a=0.):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_W(g_i=0,g_j=0,g_k=0,a=0.):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_THL(g_i=0,g_j=0,g_k=0,a=0.):
        returns()
        
    @remote_function(must_handle_array=True)
    def set_field_QT(g_i=0,g_j=0,g_k=0,a=0.):
        returns()

#    @remote_function(must_handle_array=True)
#    def set_field_E12(g_i=0,g_j=0,g_k=0,a=0.):
#        returns()

    # setter functions for wtflux and qtflux
    @remote_function
    def set_wt_surf(wtflux=0.):
        returns()

    @remote_function
    def set_wq_surf(wqflux=0.):
        returns()
        
    @remote_function()
    def get_params_grid():
        returns (i=1, j=1, k=1, xsize=1.0 | units.m, ysize=1.0 | units.m)

    #exactEnd: if true, step exactly to tend,
    # otherwise, step past tend with normal time steps (default)
    # the end time form namoptions is effective only if exactEnd is false
    @remote_function
    def evolve_model(tend=0. | units.s, exactEnd=0):
        pass


class Dales(CommonCode):
    # grid size
    itot = None
    jtot = None
    k    = None
    xsize = None
    ysize = None
    dx = None
    dy = None
    
    def __init__(self,**options):
        CommonCode.__init__(self,  DalesInterface(**options), **options)
        self.stopping_conditions = StoppingConditions(self)

        print ("DalesInterface.__init__", options)
        if 'workdir' in options:
            # print('Dales.__init__() : setting workdir.')
            self.set_workdir(options['workdir'])

    def evolve_model(self, tend, exactEnd=None):
        if exactEnd is None:
            exactEnd=self.parameters.evolve_to_exact_time
        self.overridden().evolve_model(tend, exactEnd)
            
        
    def commit_parameters(self):
        self.overridden().commit_parameters()
        self.get_params()

    def define_parameters(self, object):
        object.add_default_form_parameter(
            "input_file",
            "set the input file path",
            "namoptions.001"
        )
        object.add_method_parameter(
            "get_itot",
            None,
            "itot",
            "number of cells in the x direction",
            None
        )
        object.add_method_parameter(
            "get_jtot",
            None,
            "jtot",
            "number of cells in the y direction",
            None
        )
        object.add_method_parameter(
            "get_xsize",
            None,
            "xsize",
            "size of the grid the x direction",
            None
        )
        object.add_method_parameter(
            "get_ysize",
            None,
            "ysize",
            "size of the grid the y direction",
            None
        )
        object.add_method_parameter(
            "get_dx",
            None,
            "dx",
            "grid spacing in x direction",
            None
        )
        object.add_method_parameter(
            "get_dy",
            None,
            "dy",
            "grid spacing in y direction",
            None
        )
        object.add_interface_parameter(
            "evolve_to_exact_time",
            "flag that determines whether to evolve to exact given end time",
            False,
            "before_set_interface_parameter"
        )

    def define_properties(self, object):
        object.add_property('get_model_time', public_name = "model_time")

    def commit_grid(self):        
        self.overridden().commit_grid()
        
    def define_state(self, object):
        object.set_initial_state("UNINITIALIZED")
        object.add_transition("UNINITIALIZED", "INITIALIZED", "initialize_code")
        #~ object.add_method("!UNINITIALIZED", "before_set_parameter")
        object.add_transition("!UNINITIALIZED!STOPPED!INITIALIZED", "END", "cleanup_code")
        object.add_transition("END", "STOPPED", "stop", False)
        object.add_method("STOPPED", 'stop')
        object.add_method("UNINITIALIZED", 'stop')
        object.add_method("INITIALIZED", 'stop')

        object.add_transition("INITIALIZED","EDIT","commit_parameters")
        object.add_transition("EDIT","RUN","commit_grid")

        object.add_method("INITIALIZED", "set_input_file")

        for state in ["RUN","EVOLVED"]:
          object.add_method(state,"get_model_time")
          object.add_method(state,"get_timestep")
        for state in ["RUN","EVOLVED"]:
          object.add_method(state,"get_itot")
          #~ object.add_method(state,"get_profile_field")
          #~ object.add_method(state,"get_layer_field")
          #~ object.add_method(state,"get_volume_field")

        object.add_transition("RUN", "EVOLVED", "evolve_model", False)
        object.add_method("EVOLVED", "evolve_model")

    # wrapping functions for hiding the dummy array passed to getter functions
    def get_profile_U(self):
        return self.get_profile_U_(numpy.zeros(self.k))

    def get_profile_V(self):
        return self.get_profile_V_(numpy.zeros(self.k))

    def get_profile_W(self):
        return self.get_profile_W_(numpy.zeros(self.k))
        
    def get_profile_THL(self):
        return self.get_profile_THL_(numpy.zeros(self.k))
        
    def get_profile_QT(self):
        return self.get_profile_QT_(numpy.zeros(self.k))

    def get_profile_QL(self):
        return self.get_profile_QL_(numpy.zeros(self.k))


    def get_profile_QL_ice(self):
        return self.get_profile_QL_ice_(numpy.zeros(self.k))

    def get_profile_QR(self):
        return self.get_profile_QR_(numpy.zeros(self.k))

    def get_profile_E12(self):
        return self.get_profile_E12_(numpy.zeros(self.k))

    def get_profile_T(self):
        return self.get_profile_T_(numpy.zeros(self.k))
    
    def get_zf(self):
        return self.get_zf_(numpy.zeros(self.k))

    def get_zh(self):
        return self.get_zh_(numpy.zeros(self.k))

    def get_presh(self):
        return self.get_presh_(numpy.zeros(self.k))

    def get_presf(self):
        return self.get_presf_(numpy.zeros(self.k))

    
    # retrieve a 3D field
    # field is 'U', 'V', 'W', 'THL', 'QT', 'QL', 'E12', 'T'
    def get_field(self, field, imin=0, imax=None, jmin=0, jmax=None, kmin=0, kmax=None):
        if imax is None:
            imax = self.itot
        if jmax is None:
            jmax = self.jtot
        if kmax is None:
            kmax = self.k    

        #build index arrays
        if field in ('LWP','RWP', 'TWP'):  # 2D field
            points = numpy.mgrid[imin:imax, jmin:jmax]
            points = points.reshape(2, -1)
            i,j = points
        else: # 3D field
            points = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
            points = points.reshape(3, -1)
            i,j,k = points
            
        if field == 'U':
            field = self.get_field_U(i,j,k)
        elif field == 'V':
            field = self.get_field_V(i,j,k)
        elif field == 'W':
            field = self.get_field_W(i,j,k)
        elif field == 'THL':
            field = self.get_field_THL(i,j,k)
        elif field == 'QT':
            field = self.get_field_QT(i,j,k)
        elif field == 'QL':
            field = self.get_field_QL(i,j,k)
        elif field == 'E12':
            field = self.get_field_E12(i,j,k)
        elif field == 'T':
            field = self.get_field_T(i,j,k)
        elif field == 'LWP':                              # LWP - 2D field, k ignored
            field = self.get_field_LWP(i,j)   
            return field.reshape ((imax-imin, jmax-jmin)) # separate return here, to reshape to 2D
        elif field == 'RWP':                              # RWP - 2D field, k ignored
            field = self.get_field_RWP(i,j)   
            return field.reshape ((imax-imin, jmax-jmin)) # separate return here, to reshape to 2D
        elif field == 'TWP':                              # TWP - 2D field, k ignored
            field = self.get_field_TWP(i,j)   
            return field.reshape ((imax-imin, jmax-jmin)) # separate return here, to reshape to 2D 
        else:
            print('get_field called with undefined field', field)
            
        return field.reshape ((imax-imin, jmax-jmin, kmax-kmin))

    # set a 3D field
    # field is 'U', 'V', 'W', 'THL', 'QT'
    def set_field(self, field, a, imin=0, jmin=0, kmin=0):
        if numpy.isscalar(a):
            a = numpy.array([[[a]]])

        # set max indices from the size of a 
        imax = imin + a.shape[0]
        jmax = jmin + a.shape[1]
        kmax = kmin + a.shape[2]

        #build index arrays
        points = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
        points = points.reshape(3, -1)
        i,j,k = points

        a = a.reshape(-1) # make a one-dimensional
                          # numpy uses C ordering by default - last index varies fastest
        #print('i', i)
        #print('j', j)
        #print('k', k)
        #print('a', a)
        
        if field == 'U':
            self.set_field_U(i,j,k,a)
        elif field == 'V':
            self.set_field_V(i,j,k,a)
        elif field == 'W':
            self.set_field_W(i,j,k,a)
        elif field == 'THL':
            self.set_field_THL(i,j,k,a)
        elif field == 'QT':
            self.set_field_QT(i,j,k,a)
        else:
            print('set_field called with undefined field', field)
                


    # retrieve a 1D vertical profile - wrapper function consistent with get_field
    # field is 'U', 'V', 'W', 'THL', 'QT', 'QL', 'E12', 'T'
    def get_profile(self, field):
        if field == 'U':
            profile = self.get_profile_U_(numpy.zeros(self.k))
        elif field == 'V':
            profile = self.get_profile_V_(numpy.zeros(self.k))
        elif field == 'W':
            profile = self.get_profile_W_(numpy.zeros(self.k))
        elif field == 'THL':
            profile = self.get_profile_THL_(numpy.zeros(self.k))
        elif field == 'QT':
            profile = self.get_profile_QT_(numpy.zeros(self.k))
        elif field == 'QL':
            profile = self.get_profile_QL_(numpy.zeros(self.k))
        elif field == 'E12':
            profile = self.get_profile_E12_(numpy.zeros(self.k))
        elif field == 'T':
            profile = self.get_profile_T_(numpy.zeros(self.k))
        else:
            print('get_profile called with undefined field', field)
            
        return profile


    # get parameters from the fortran code, store them in the Dales interface object
    # called from commit_parameters()
    def get_params(self):
        self.itot,self.jtot,self.k,self.xsize,self.ysize = self.get_params_grid()
        self.dx = self.xsize/self.itot
        self.dy = self.ysize/self.jtot
        self.zf = self.get_zf()
        self.zh = self.get_zh()
        
    def get_itot(self):
        return self.itot
    def get_jtot(self):
        return self.jtot
    def get_xsize(self):
        return self.xsize
    def get_ysize(self):
        return self.ysize
    def get_dx(self):
        return self.dx
    def get_dy(self):
        return self.dy
