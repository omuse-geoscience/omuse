import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option

class SwanInterface(CodeInterface, 
                      CommonCodeInterface,
                      LiteratureReferencesMixIn):
    """
    
    SWAN - 

    .. [#] swanmodel.sf.net
    
    """
    use_modules=['StoppingConditions','swan_interface']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'swan_worker'

    @remote_function
    def initialize_code(coordinates="cartesian", mode="stationary", grid_type="regular",input_grid_type="regular"):
        returns ()

    @remote_function
    def initialize_grids():
        returns ()

    @remote_function
    def commit_grids():
        returns ()


    @remote_function
    def get_wlev():
        returns (wlev=0. | units.m)
    @remote_function
    def set_wlev(wlev=0. | units.m):
        returns ()

    @remote_function
    def get_grav():
        returns (grav=9.81 | units.m/units.s**2)
    @remote_function
    def set_grav(grav=9.81 | units.m/units.s**2):
        returns ()

    @remote_function
    def get_rho():
        returns (rho=1025. | units.kg/units.m**3)
    @remote_function
    def set_rho(rho=1025. | units.kg/units.m**3):
        returns ()

    @remote_function
    def get_cdcap():
        returns (rho=99999. )
    @remote_function
    def set_cdcap(rho=99999. ):
        returns ()

    @remote_function
    def get_rearth():
        returns (earth_radius=6.371e6 | units.m)
    @remote_function
    def set_rearth(earth_radius=6.371e6 | units.m ):
        returns ()

    @remote_function
    def get_grid_xpc():
        returns (grid_xpc=0. | units.m)
    @remote_function
    def set_grid_xpc(grid_xpc=0. | units.m ):
        returns ()
    
    @remote_function
    def get_grid_ypc():
        returns (grid_ypc=0. | units.m)
    @remote_function
    def set_grid_ypc(grid_ypc=0. | units.m ):
        returns ()
    
    @remote_function
    def get_grid_alpc():
        returns (grid_alpc=0. | units.deg)
    @remote_function
    def set_grid_alpc(grid_alpc=0. | units.deg ):
        returns ()
    
    @remote_function
    def get_grid_xlenc():
        returns (grid_xlenc=0. | units.m)
    @remote_function
    def set_grid_xlenc(grid_xlenc=0. | units.m ):
        returns ()
    
    @remote_function
    def get_grid_ylenc():
        returns (grid_ylenc=0. | units.m)
    @remote_function
    def set_grid_ylenc(grid_ylenc=0. | units.m ):
        returns ()
    
    @remote_function
    def get_grid_mxc():
        returns (grid_mxc=0)
    @remote_function
    def set_grid_mxc(grid_mxc=0 ):
        returns ()

    @remote_function
    def get_grid_myc():
        returns (grid_myc=0)
    @remote_function
    def set_grid_myc(grid_myc=0 ):
        returns ()

    @remote_function
    def get_nfreq():
        returns (nfreq=0)
    @remote_function
    def set_nfreq(nfreq=0 ):
        returns ()

    @remote_function
    def get_ndir():
        returns (ndir=0)
    @remote_function
    def set_ndir(ndir=0 ):
        returns ()

    @remote_function
    def get_flow():
        returns (flow=0. | units.s**-1)
    @remote_function
    def set_flow(flow=0. | units.s**-1 ):
        returns ()

    @remote_function
    def get_fhigh():
        returns (fhigh=0. | units.s**-1)
    @remote_function
    def set_fhigh(fhigh=0. | units.s**-1 ):
        returns ()

    @remote_function
    def get_input_xp():
        returns (input_xp=0. | units.m)
    @remote_function
    def set_input_xp(input_xp=0. | units.m ):
        returns ()
    @remote_function
    def get_input_yp():
        returns (input_yp=0. | units.m)
    @remote_function
    def set_input_yp(input_yp=0. | units.m ):
        returns ()
    @remote_function
    def get_input_dx():
        returns (input_dx=0. | units.m)
    @remote_function
    def set_input_dx(input_dx=0. | units.m ):
        returns ()
    @remote_function
    def get_input_dy():
        returns (input_dy=0. | units.m)
    @remote_function
    def set_input_dy(input_dy=0. | units.m ):
        returns ()
    @remote_function
    def get_input_alp():
        returns (input_alp=0. | units.deg)
    @remote_function
    def set_input_alp(input_alp=0. | units.deg ):
        returns ()
    @remote_function
    def get_input_mx():
        returns (input_mx=0)
    @remote_function
    def set_input_mx(input_mx=0 ):
        returns ()
    @remote_function
    def get_input_my():
        returns (input_my=0)
    @remote_function
    def set_input_my(input_my=0 ):
        returns ()


        
    @remote_function
    def get_calc_mode():
        returns (mode="s")
    @remote_function
    def get_grid_type():
        returns (grid_type="s")
    @remote_function
    def get_input_grid_type():
        returns (input_grid_type="s")
    @remote_function
    def get_coord():
        returns (coordinates="s")
    @remote_function
    def get_proj_method():
        returns (projection_method="s")
    @remote_function
    def get_ndim():
        returns (number_of_dimensions=0)
