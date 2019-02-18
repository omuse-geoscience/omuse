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
import f90nml

from amuse.units import trigo

import numpy

class DalesInterface(CodeInterface,
                     CommonCodeInterface,
                     StoppingConditionInterface,
                     LiteratureReferencesMixIn):
    """

    DALES - Dutch Atmospheric Large Eddy Simulation

    .. [#] Heus et al., Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications,  Geoscientific Model Development 3, 415, (2010)

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
    def set_qt_forcing(forcing_type=0):
        pass

    @remote_function
    def get_exact_end():
        returns (exactEndFlag=False)

    @remote_function
    def set_exact_end(exactEndFlag=False):
        pass

    @remote_function
    def set_workdir(directory="./"):
        pass

    @remote_function
    def get_workdir():
        returns (directory="s")

    @remote_function
    def set_start_date(date=0):
        pass

    @remote_function
    def set_start_time(time=0):
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
    def set_tendency_surface_pressure(p = 1.0e5 | units.Pa/units.s):
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
        returns (out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_QL_(k=0):
        returns (out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_QL_ice_(k=0):
        returns (out=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_profile_QR_(k=0):
        returns (out=0. | units.mfu)
        
    @remote_function(must_handle_array=True)
    def get_profile_E12_(k=0):
        returns (out=0. | units.m/units.s)

    @remote_function(must_handle_array=True)
    def get_profile_T_(k=0):
        returns (out=0. | units.K)

    # getter for cloud fraction. Uses the index array to define slabs.
    @remote_function(must_handle_array=True)
    def get_cloudfraction(k=0):
        returns (out=0. | units.m**2/units.m**2)

    # getter for accumulated surface rain flux
    @remote_function()
    def get_rain():
        returns (out=0. | units.kg / units.m**2)
        
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
    def set_tendency_U(a=0. | units.m/units.s**2):
        returns () 

    @remote_function(must_handle_array=True)
    def set_tendency_V(a=0. | units.m/units.s**2):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_THL(a=0. | units.K/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_QT(a=0. | units.mfu/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_QL(a=0. | units.mfu/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def set_ref_profile_QL(a=0.| units.mfu):
        returns ()

    @remote_function(must_handle_array=True)
    def set_qt_variability_factor(a=0. | 1 / units.s):
        returns ()

# indexed getter/setter functions for vertical tendencies / forcings
    @remote_function(must_handle_array=True)
    def set_tendency_U_(g_i=0, a=0. | units.m/units.s**2):
        returns () 

    @remote_function(must_handle_array=True)
    def set_tendency_V_(g_i=0, a=0. | units.m/units.s**2):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_THL_(g_i=0, a=0.| units.K/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def set_tendency_QT_(g_i=0,a=0.| units.mfu/units.s):
        returns ()

    @remote_function(must_handle_array=True)
    def get_tendency_U_(g_i=0):
        returns (a=0.| units.m/units.s**2) 

    @remote_function(must_handle_array=True)
    def get_tendency_V_(g_i=0):
        returns (a=0.| units.m/units.s**2)

    @remote_function(must_handle_array=True)
    def get_tendency_THL_(g_i=0):
        returns (a=0.| units.K/units.s)

    @remote_function(must_handle_array=True)
    def get_tendency_QT_(g_i=0):
        returns (a=0.| units.mfu/units.s)

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
        returns (a=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_field_QL(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_field_Qsat(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.mfu)

    @remote_function(must_handle_array=True)
    def get_field_E12(g_i=0,g_j=0,g_k=0):
        returns (a=0. | units.m/units.s)

    @remote_function(must_handle_array=True)
    def get_field_T(g_i=0,g_j=0,g_k=0):
        returns (a=0.| units.K)

    @remote_function(must_handle_array=True)
    def get_field_LWP(g_i=0,g_j=0):
        returns (a=0. | units.kg/units.m**2)

    @remote_function(must_handle_array=True)
    def get_field_RWP(g_i=0,g_j=0):
        returns (a=0. | units.kg/units.m**2) 

    @remote_function(must_handle_array=True)
    def get_field_TWP(g_i=0,g_j=0):
        returns (a=0. | units.kg/units.m**2)
        
    # setter functions for 3D fields usning index arrays
    @remote_function(must_handle_array=True)
    def set_field_U(g_i=0,g_j=0,g_k=0,a=0. | units.m/units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_V(g_i=0,g_j=0,g_k=0,a=0. | units.m/units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_W(g_i=0,g_j=0,g_k=0,a=0. | units.m/units.s):
        returns()

    @remote_function(must_handle_array=True)
    def set_field_THL(g_i=0,g_j=0,g_k=0,a=0. | units.K):
        returns()
        
    @remote_function(must_handle_array=True)
    def set_field_QT(g_i=0,g_j=0,g_k=0,a=0. | units.mfu):
        returns()

#    @remote_function(must_handle_array=True)
#    def set_field_E12(g_i=0,g_j=0,g_k=0,a=0.):
#        returns()

    # setter functions for wtflux and qtflux
    @remote_function
    def set_wt_surf(wtflux=0. | units.m * units.s**-1 * units.K):
        returns()

    @remote_function
    def set_wq_surf(wqflux=0. | units.m / units.s):
        returns()

    # setter functions for momentum and heat roughness
    @remote_function
    def set_z0m_surf(z0=0. | units.m):
        returns()

    @remote_function
    def set_z0h_surf(z0=0. | units.m):
        returns()
        
    @remote_function()
    def get_params_grid():
        returns (i=1, j=1, k=1, xsize=1.0 | units.m, ysize=1.0 | units.m)

    #exactEnd: if true, step exactly to tend,
    # otherwise, step past tend with normal time steps (default)
    # the end time form namoptions is effective only if exactEnd is false
    @remote_function
    def evolve_model(tend=0. | units.s, exactEnd=0):
        returns (walltime=0. | units.s)

    @remote_function
    def write_restart():
        returns()

class Dales(CommonCode):
    
    QT_FORCING_GLOBAL=0
    QT_FORCING_LOCAL=1
    QT_FORCING_VARIANCE=2

    def __init__(self,**options):
        self._input_file=options.get("input_file", "namoptions.001")

        CommonCode.__init__(self,  DalesInterface(**options), **options)
        self.stopping_conditions = StoppingConditions(self)

        # grid size
        self.itot = None
        self.jtot = None
        self.k = None
        self.xsize = None
        self.ysize = None
        self.dx = None
        self.dy = None

        if 'workdir' in options:
            # print('Dales.__init__() : setting workdir.')
            self.set_workdir(options['workdir'])            

        self._workdir=self.get_workdir()

        self.read_input_file()

    def read_input_file(self):
        inputfile=os.path.join(self._workdir,self.parameters.input_file)

        self.params = f90nml.read(inputfile)
        self.parameters.grid_size_xy = self.params["DOMAIN"]["itot"],self.params["DOMAIN"]["jtot"]
        self.parameters.grid_extent_xy = [self.params["DOMAIN"]["xsize"],self.params["DOMAIN"]["ysize"]] | units.m

        #print ("DalesInterface.__init__", options)
        
    def commit_parameters(self):
      
        if self.parameters.input_file != self._input_file:
          print "rereading namelist options file..."
          self.read_input_file()

        inputfile=os.path.join(self._workdir,self.parameters.input_file)
        dalesinputfile=os.path.join(self._workdir,"_"+self.parameters.input_file)
        if not os.path.exists(os.path.join(self._workdir,inputfile)):
            raise Exception("File %s not found"%inputfile)
        
        
        patch=dict()
        patch["RUN"]=dict(  lwarmstart=self.parameters.restart_flag,
                            startfile=self.parameters.restart_file, 
                            trestart=self.parameters.trestart.value_in(units.s) )
        if self.parameters.grid_size_xy:
            patch["DOMAIN"]=dict( itot=self.parameters.grid_size_xy[0],
                                  jtot=self.parameters.grid_size_xy[1] )
        f90nml.patch(inputfile,patch,dalesinputfile)
        self.set_input_file(dalesinputfile)

        # print "code options written to %s"%dalesinputfile

        # check this...
        dt = self.parameters.starttime
        if dt != 0:
            self.set_start_date(10000 * dt.year + 100 * dt.month + dt.day)
            self.set_start_time(10000 * dt.hour + 100 * dt.minute + dt.second)
        self.set_qt_forcing(self.parameters.qt_forcing)
        
        self.overridden().commit_parameters()
        self.get_params()

    def define_parameters(self, object):
        object.add_interface_parameter(
            "input_file",
            "the input file name",
            self._input_file,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "restart_flag",
            "warm start from restart file if True",
            False,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "restart_file",
            "restart file name",
            "initdlatestx000y000.001",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "trestart",
            "(simulation) time between writing restart files",
            3600 | units.s,
            "before_set_interface_parameter"
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
        object.add_boolean_parameter(
            "get_exact_end",
            "set_exact_end",
            "exactEndFlag",
            "parameter specifying whether evolve should end exactly at target time ",
            False
        )
        object.add_interface_parameter(
            "starttime",
            "absolute model start datetime",
            0,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "qt_forcing",
            "qt forcing type",
            0,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "grid_size_xy",
            "tuple of number grid cells in x and y direction (None means use namelist default)",
            None,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "grid_extent_xy",
            "tuple or vector of physical grid size in x and y direction (None means use namelist default)",
            None,
            "before_set_interface_parameter"
        )


    def define_properties(self, object):
        object.add_property('get_model_time', public_name = "model_time")
        object.add_property('get_timestep', public_name = "timestep")

    def commit_grid(self):        
        self.overridden().commit_grid()
        
    def define_state(self, object):
        object.set_initial_state("UNINITIALIZED")
        object.add_transition("UNINITIALIZED", "INITIALIZED", "initialize_code")
        object.add_method("!UNINITIALIZED", "before_get_parameter")
        object.add_transition("!UNINITIALIZED!STOPPED!INITIALIZED", "END", "cleanup_code")
        object.add_transition("END", "STOPPED", "stop", False)
        object.add_method("STOPPED", 'stop')
        object.add_method("UNINITIALIZED", 'stop')
        object.add_method("INITIALIZED", 'stop')

        object.add_transition("INITIALIZED","EDIT","commit_parameters")
        object.add_transition("EDIT","RUN","commit_grid")

        object.add_method("INITIALIZED", "set_input_file")
        object.add_method("INITIALIZED", "set_start_date_time")

        for state in ["EDIT","RUN","EVOLVED"]:
          object.add_method(state,"get_model_time")
          object.add_method(state,"get_timestep")
        for state in ["RUN","EVOLVED"]:
          object.add_method(state,"get_itot")

        object.add_transition("RUN", "EVOLVED", "evolve_model", False)
        object.add_method("EVOLVED", "evolve_model")

    # wrapping functions for hiding the dummy array passed to getter functions
    def get_profile_U(self):
        return self.get_profile_U_(numpy.arange(1,self.k+1))

    def get_profile_V(self):
        return self.get_profile_V_(numpy.arange(1,self.k+1))

    def get_profile_W(self):
        return self.get_profile_W_(numpy.arange(1,self.k+1))
        
    def get_profile_THL(self):
        return self.get_profile_THL_(numpy.arange(1,self.k+1))
        
    def get_profile_QT(self):
        return self.get_profile_QT_(numpy.arange(1,self.k+1))

    def get_profile_QL(self):
        return self.get_profile_QL_(numpy.arange(1,self.k+1))


    def get_profile_QL_ice(self):
        return self.get_profile_QL_ice_(numpy.arange(1,self.k+1))

    def get_profile_QR(self):
        return self.get_profile_QR_(numpy.arange(1,self.k+1))

    def get_profile_E12(self):
        return self.get_profile_E12_(numpy.arange(1,self.k+1))

    def get_profile_T(self):
        return self.get_profile_T_(numpy.arange(1,self.k+1))
    
    def get_zf(self):
        return self.get_zf_(numpy.arange(1,self.k+1))

    def get_zh(self):
        return self.get_zh_(numpy.arange(1,self.k+1))

    def get_presh(self):
        return self.get_presh_(numpy.arange(1,self.k+1))

    def get_presf(self):
        return self.get_presf_(numpy.arange(1,self.k+1))

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
        elif field == 'Qsat':
            field = self.get_field_Qsat(i,j,k)
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

        # set max indices from the size of a 
        try:
          imax = imin + a.shape[0]
          jmax = jmin + a.shape[1]
          kmax = kmin + a.shape[2]
          i,j,k = numpy.mgrid[imin:imax, jmin:jmax, kmin:kmax]
          i=i.flatten()
          j=j.flatten()
          k=k.flatten()
          a=a.flatten()
        except:
          i=imin
          j=jmin
          k=kmin

        #build index arrays
        
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
            profile = self.get_profile_U_(numpy.arange(1,self.k+1))
        elif field == 'V':
            profile = self.get_profile_V_(numpy.arange(1,self.k+1))
        elif field == 'W':
            profile = self.get_profile_W_(numpy.arange(1,self.k+1))
        elif field == 'THL':
            profile = self.get_profile_THL_(numpy.arange(1,self.k+1))
        elif field == 'QT':
            profile = self.get_profile_QT_(numpy.arange(1,self.k+1))
        elif field == 'QL':
            profile = self.get_profile_QL_(numpy.arange(1,self.k+1))
        elif field == 'E12':
            profile = self.get_profile_E12_(numpy.arange(1,self.k+1))
        elif field == 'T':
            profile = self.get_profile_T_(numpy.arange(1,self.k+1))
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
    def get_kmax(self):
        return self.k
    def get_xsize(self):
        return self.xsize
    def get_ysize(self):
        return self.ysize
    def get_dx(self):
        return self.dx
    def get_dy(self):
        return self.dy

    def get_grid_range(self):
        return (0,self.get_itot()-1,0,self.get_jtot()-1,0,self.get_kmax()-1)

    def get_grid_position(self,i,j,k):
        return i*self.dx,j*self.dy,self.zf[k]

    def get_surface_field_position(self,i,j):
        return i*self.dx,j*self.dy

    def get_profile_grid_range(self):
        return (0,self.get_kmax()-1)

    def get_surface_grid_range(self):
        return ()

    def get_surface_field_grid_range(self):
        return (0,self.get_itot()-1,0,self.get_jtot()-1)

    def get_profile_grid_position(self,k):
        return self.zf[k]
        
    def define_grids(self, object):
        object.define_grid('grid',axes_names = "xyz", grid_class=datamodel.RectilinearGrid,state_guard="before_new_set_instance")
        object.set_grid_range('grid', 'get_grid_range')
        object.add_getter('grid', 'get_grid_position', names="xyz")
        for x in ['U', 'V', 'W', 'THL', 'QT', 'QL', 'E12', 'T']:
            object.add_getter('grid', 'get_field_'+x,  names=[x])
        for x in ['U', 'V', 'W', 'THL', 'QT']:
            object.add_setter('grid', 'set_field_'+x,  names=[x])

        object.define_grid('profile_grid',axes_names = "z", grid_class=datamodel.RectilinearGrid,state_guard="before_new_set_instance")
        object.set_grid_range('profile_grid', 'get_profile_grid_range')
        object.add_getter('profile_grid', 'get_profile_grid_position', names="z")
        for x in ['U', 'V', 'W', 'THL', 'QT', 'QL', 'QL_ice', 'E12', 'T']:
            object.add_getter('profile_grid', 'get_profile_'+x+'_',  names=[x])

        object.define_grid('forcings',axes_names = "z", grid_class=datamodel.RectilinearGrid,state_guard="before_new_set_instance")
        object.set_grid_range('forcings', 'get_profile_grid_range')
        object.add_getter('forcings', 'get_profile_grid_position', names="z")
        for x in ['U', 'V', 'THL', 'QT']:
            object.add_getter('forcings', 'get_tendency_'+x+'_',  names=['tendency_'+x])
            object.add_setter('forcings', 'set_tendency_'+x+'_',  names=['tendency_'+x])

        object.define_grid('surface', grid_class=datamodel.RectilinearGrid,state_guard="before_new_set_instance")
        object.set_grid_range('surface', 'get_surface_grid_range')
        object.add_getter('surface', 'get_surface_pressure',  names=['pressure'])
        object.add_getter('surface', 'get_rain',  names=['rain'])

        object.define_grid('surface_field', grid_class=datamodel.RectilinearGrid,state_guard="before_new_set_instance")
        object.set_grid_range('surface_field', 'get_surface_field_grid_range')
        object.add_getter('surface_field', 'get_surface_field_position', names="xy")
        object.add_getter('surface_field', 'get_field_LWP',  names=['LWP'])
        object.add_getter('surface_field', 'get_field_TWP',  names=['TWP'])
        object.add_getter('surface_field', 'get_field_RWP',  names=['RWP'])
