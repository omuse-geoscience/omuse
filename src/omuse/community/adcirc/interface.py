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

from write_grid import adcirc_grid_writer,adcirc_parameter_writer

from amuse.datamodel.staggeredgrid import StaggeredGrid

class AdcircInterface(CodeInterface, 
                      CommonCodeInterface,
                      StoppingConditionInterface,
                      LiteratureReferencesMixIn):
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
    def get_rootdir():
        returns (rootdir="s")
    @remote_function
    def set_rootdir(rootdir="."):
        returns ()

    @remote_function
    def commit_grid():
        pass

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
    def get_node_eta(index=0):
        returns (eta=0. | units.m)
    @remote_function(can_handle_array=True)
    def set_node_eta(index=0,eta=0. | units.m):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_vx(index=0):
        returns (vx=0.| units.m/units.s)
    @remote_function(can_handle_array=True)
    def set_node_vx(index=0,vx=0.| units.m/units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_vy(index=0):
        returns (vy=0.| units.m/units.s)
    @remote_function(can_handle_array=True)
    def set_node_vy(index=0,vy=0.| units.m/units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_deta_dt(index=0):
        returns (deta_dt=0. | units.m/units.s)
    @remote_function(can_handle_array=True)
    def set_node_deta_dt(index=0,deta_dt=0. | units.m/units.s):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_status(index=0):
        returns (status='s')
    @remote_function(can_handle_array=True)
    def set_node_status(index=0,status='wet'):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_coriolis_f(index=0):
        returns (coriolis_f=0. | units.s**-1)
    @remote_function(can_handle_array=True)
    def set_node_coriolis_f(index=0,coriolis_f=0. | units.s**-1):
        returns ()
        
    @remote_function(can_handle_array=True)
    def get_node_wind_stress(index=0):
        returns (tau_x=0. | units.Pa,tau_y=0. | units.Pa)
    @remote_function(can_handle_array=True)
    def set_node_wind_stress(index=0,tau_x=0.| units.Pa,tau_y=0. | units.Pa):
        returns ()        

    @remote_function(can_handle_array=True)
    def get_node_wave_stress(index=0):
        returns (tau_x=0. | units.Pa,tau_y=0. | units.Pa)
    @remote_function(can_handle_array=True)
    def set_node_wave_stress(index=0,tau_x=0.| units.Pa,tau_y=0. | units.Pa):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_atmospheric_pressure(index=0):
        returns (pressure=0. | units.mbar)
    @remote_function(can_handle_array=True)
    def set_node_atmospheric_pressure(index=0,pressure=0. | units.mbar):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_tidal_potential(index=0):
        returns (tidal_potential=0. | (units.m**2/units.s**2) )
    @remote_function(can_handle_array=True)
    def set_node_tidal_potential(index=0,tidal_potential=0. | (units.m**2/units.s**2) ):
        returns ()

    @remote_function(can_handle_array=True)
    def get_node_surface_heat_flux(index=0):
        returns (shf=0. | units.W/units.m**2)
    @remote_function(can_handle_array=True)
    def set_node_surface_heat_flux(index=0,shf=0. | units.W/units.m**2):
        returns ()


    @remote_function(can_handle_array=True)
    def get_node_position(index='i'):
        returns (x=0.| units.m,y=0. | units.m)
    @remote_function(can_handle_array=True)
    def get_node_coordinates(index='i'):
        returns (lon=0.| units.rad,lat=0. | units.rad)

    @remote_function(can_handle_array=True)
    def get_node_depth(index='i'):
        returns (depth=0.| units.m)

    @remote_function(can_handle_array=True)
    def get_node_sigma(index='i',zindex='i'):
        returns (sigma='d',z='d' | units.m)

    @remote_function(can_handle_array=True)
    def get_element_position(index='i'):
        returns (x=0.| units.m,y=0. | units.m)
    @remote_function(can_handle_array=True)
    def get_element_coordinates(index='i'):
        returns (lon=0.| units.rad,lat=0. | units.rad)

    @remote_function(can_handle_array=True)
    def get_node_velocities_3d(index='i',zindex='i'):
        returns (vx=0.| units.m/units.s, vy=0. | units.m/units.s, vz=0. | units.m/units.s)

    @remote_function(can_handle_array=True)
    def get_node_temperature_3d(index='i',zindex='i'):
        returns (temperature=0. | units.Celsius)

    @remote_function(can_handle_array=True)
    def get_element_nodes(index=0):
        returns (n1=0,n2=0,n3=0)

    @remote_function(can_handle_array=True)
    def get_element_status(index=0):
        returns (status='s')
    @remote_function(can_handle_array=True)
    def set_element_status(index=0,status='wet'):
        returns ()

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

    @remote_function
    def get_number_of_vertical_nodes():
        returns (n_znodes=0)

    @remote_function(can_handle_array=True)
    def get_elevation_boundary_node(index=0,index_of_segment=0):
        returns (node_index=0)

    @remote_function(can_handle_array=True)
    def set_elevation_boundary_eta(index=0,eta=0. | units.m,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_elevation_boundary_eta(index=0,index_of_segment=0):
        returns (eta=0. | units.m)

    @remote_function(can_handle_array=True)
    def set_boundary_lnm(index=0,lnm=0. | units.m,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_boundary_lnm(index=0,index_of_segment=0):
        returns (lnm=0. | units.m)

    @remote_function(can_handle_array=True)
    def set_boundary_salinity(index=0,zindex=0,salinity=0. | units.psu,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_boundary_salinity(index=0,zindex=0,index_of_segment=0):
        returns (salinity=0. | units.psu)

    @remote_function(can_handle_array=True)
    def set_boundary_temperature(index=0,zindex=0,temperature=0. | units.Celsius,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_boundary_temperature(index=0,zindex=0,index_of_segment=0):
        returns (temperature=0. | units.Celsius)

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

    @remote_function(can_handle_array=True)
    def set_flow_boundary_fluxx(index=0,flux_x=0. | units.m**2/units.s,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_flow_boundary_fluxx(index=0,index_of_segment=0):
        returns (flux_x=0. | units.m**2/units.s)

    @remote_function(can_handle_array=True)
    def set_flow_boundary_fluxy(index=0,flux_y=0. | units.m**2/units.s,index_of_segment=0):
        returns ()

    @remote_function(can_handle_array=True)
    def get_flow_boundary_fluxy(index=0,index_of_segment=0):
        returns (flux_y=0. | units.m**2/units.s)

    @remote_function
    def get_use_interface_elevation_boundary():
        returns (use_interface_elevation_boundary='b')
    @remote_function
    def set_use_interface_elevation_boundary(use_interface_elevation_boundary='b'):
        returns ()

    @remote_function
    def get_use_interface_flow_boundary():
        returns (use_interface_flow_boundary='b')
    @remote_function
    def set_use_interface_flow_boundary(use_interface_flow_boundary='b'):
        returns ()

    @remote_function
    def get_use_interface_met_forcing():
        returns (use_interface_met_forcing='b')
    @remote_function
    def set_use_interface_met_forcing(use_interface_met_forcing='b'):
        returns ()

    @remote_function
    def get_use_interface_wave_forcing():
        returns (use_interface_wave_forcing='b')
    @remote_function
    def set_use_interface_wave_forcing(use_interface_wave_forcing='b'):
        returns ()

    @remote_function
    def get_use_interface_tidal_forcing():
        returns (use_interface_tidal_forcing='b')
    @remote_function
    def set_use_interface_tidal_forcing(use_interface_tidal_forcing='b'):
        returns ()

    @remote_function
    def get_use_interface_surface_heat_forcing():
        returns (use_interface_surface_heat_forcing='b')
    @remote_function
    def set_use_interface_surface_heat_forcing(use_interface_surface_heat_forcing='b'):
        returns ()

    @remote_function
    def get_use_interface_lnm_boundary():
        returns (use_interface_lnm_boundary='b')
    @remote_function
    def set_use_interface_lnm_boundary(use_interface_lnm_boundary='b'):
        returns ()

    @remote_function
    def get_use_interface_salinity_boundary():
        returns (use_interface_salinity_boundary='b')
    @remote_function
    def set_use_interface_salinity_boundary(use_interface_salinity_boundary='b'):
        returns ()

    @remote_function
    def get_use_interface_temperature_boundary():
        returns (use_interface_temperature_boundary='b')
    @remote_function
    def set_use_interface_temperature_boundary(use_interface_temperature_boundary='b'):
        returns ()

    @remote_function
    def get_reference_pressure():
        returns (pressure=0. | units.mbar)

class Adcirc(CommonCode):
  
    MODE_2D = "2D"
    MODE_3D = "3D"

    bottom_friction_laws={"linear":0,"quadratic":1,"hybrid":2}     
    bottom_friction_laws_3d={"noslip":0,"linear":1,"loglayer":2,"quadratic":3}     
    baroclinic_density_forcings={"none":0, "sigmaT":1, "salinity":2,"temperature":3,"salinity+temperature":4}
    vertical_grid_types={"user":0, "uniform":1, "log":2, "loglinear":3, "doublelog":4, "pgrid":5, "sine":6}
    vertical_eddy_viscosity_formulations={"user":0,"constant":1}
    equations_of_state={"CushmanRoisin":1,"McDougall":2,"unesco":3}
    
    def __init__(self, mode=MODE_2D, coordinates="cartesian", **options):
        self.mode=mode
        self.coordinates=coordinates
        CommonCode.__init__(self,  AdcircInterface(**options), **options)
        self._nodes=None
        self._elements=None
        self._elev_boundary=None
        self._flow_boundary=None
        self._parameters=None
        self.stopping_conditions = StoppingConditions(self)
       
    def get_coordinates(self):
        return self.coordinates
       
    def assign_grid_and_boundary(self,nodes,elements,elev_boundary, flow_boundary):
        self._nodes=nodes
        self._elements=elements
        self._elev_boundary=elev_boundary
        self._flow_boundary=flow_boundary

    def commit_parameters(self):      
        if self.parameters.use_interface_parameters:
          param=adcirc_parameter_writer()
          if self._parameters is not None:
            param.parameters=self._parameters
          else:
            self._parameters=param.parameters

          if self.parameters.calculate_coriolis and not self.coordinates=="spherical":
            raise Exception("calculate_coriolis needs spherical coordinates")

          if self.parameters.bottom_friction_law not in self.bottom_friction_laws:
            raise Exception("invalid/ unimplemented bottom friction law: %s"%self.parameters.bottom_friction_law)
          if self.parameters.bottom_friction_law_3d not in self.bottom_friction_laws_3d:
            raise Exception("invalid/ unimplemented 3D bottom friction law: %s"%self.parameters.bottom_friction_law_3d)
          if self.parameters.baroclinic_density_forcing not in self.baroclinic_density_forcings:
            raise Exception("invalid/ unimplemented baroclinic forcing: %s"%self.parameters.baroclinic_density_forcing)
          if self.parameters.vertical_grid_type not in self.vertical_grid_types:
            raise Exception("invalid vertical grid type: %s"%self.parameters.vertical_grid_type)
          if self.parameters.vertical_eddy_viscosity_formulation not in self.vertical_eddy_viscosity_formulations:
            raise Exception("invalid vertical eddy viscosity, check! %s"%self.parameters.vertical_eddy_viscosity_formulation)
          if self.parameters.equation_of_state not in self.equations_of_state:
            raise Exception("invalid equation of state: %s"%self.paramters.equation_of_state)

          if self.mode==self.MODE_2D and self.parameters.baroclinic_flag==False:
            IM,IDEN=0,0
          elif self.mode==self.MODE_3D and self.parameters.baroclinic_flag==False:   
            IM,IDEN=1,0
          elif self.mode==self.MODE_2D and self.parameters.baroclinic_flag==True:
            IM=20
            IDEN=self.baroclinic_density_forcings[self.parameters.baroclinic_density_forcing]   
          elif self.mode==self.MODE_3D and self.parameters.baroclinic_flag==True:
            IM=21
            IDEN=self.baroclinic_density_forcings[self.parameters.baroclinic_density_forcing] 
          
          if (self.parameters.use_interface_lnm_boundary or
              self.parameters.use_interface_temperature_boundary or
              self.parameters.use_interface_salinity_boundary 
             ) and not (self.mode==self.MODE_3D and self.parameters.baroclinic_flag==True):
               raise Exception("interface lnm, salinity or temperature boundary needs 3D baroclinic")
            
          if self.parameters.bottom_friction_law_3d == "noslip":
            KP=0.
          elif self.parameters.bottom_friction_law_3d == "linear":
            KP=self.parameters.linear_3d_bottom_friction_coeff.value_in(units.m/units.s)
          else:
            KP=self.parameters.quadratic_3d_bottom_friction_coeff
          NFEN=len(self.parameters.sigma_levels) if self.parameters.vertical_grid_type=="user" else self.parameters.number_of_vertical_levels
          param.update( IM=IM,
                        IDEN=IDEN,
                        NCOR=1 if self.parameters.calculate_coriolis else 0,
                        ESLM=self.parameters.A_H.value_in(units.m**2/units.s),
                        SLAM0=trigo.in_deg(self.parameters.central_longitude),
                        SFEA0=trigo.in_deg(self.parameters.central_latitude),
                        ICS=1 if self.parameters.coordinates=="cartesian" else 2,
                        DTDP=(-1 if self.parameters.use_predictor_corrector else 1)*self.parameters.timestep.value_in(units.s),
                        NOLIBF=self.bottom_friction_laws[self.parameters.bottom_friction_law],
                        TAU0=self.parameters.GWCE_weighting_factor,
                        TAU=self.parameters.linear_bottom_friction_coeff.value_in(units.s**-1),
                        CF=self.parameters.quadratic_bottom_friction_coeff,
                        STATIM=self.parameters.begin_time.value_in(units.day),
                        CONVCR=self.parameters.convergence_criterion,
                        ITMAX=self.parameters.maximum_iterations,
                        ISLDIA=self.parameters.solver_verbosity,
# NOLIFA, NOLICA, NOLICAT should all either 0 or >0...
                        NOLIFA=self.parameters.finite_amplitude_term_parameter,
                        NOLICA=self.parameters.spatial_derivative_advective_term_parameter,
                        NOLICAT=self.parameters.time_derivative_advective_term_parameter,
                        ISLIP=self.bottom_friction_laws_3d[self.parameters.bottom_friction_law_3d],
                        KP=KP,
                        IGC=self.vertical_grid_types[self.parameters.vertical_grid_type],
                        NFEN=NFEN,
                        SIGMA=self.parameters.sigma_levels,
                        IEVC=self.vertical_eddy_viscosity_formulations[self.parameters.vertical_eddy_viscosity_formulation],
                        EVMIN=self.parameters.minimum_vertical_eddy_viscosity.value_in(units.m**2/units.s),
                        EVCON=self.parameters.vertical_eddy_viscosity.value_in(units.m**2/units.s),
                        EVTOT=self.parameters.vertical_eddy_viscosities.value_in(units.m**2/units.s),
                        NLSD=self.parameters.lateral_salinity_diffusion_coefficient.value_in(units.m**2/units.s),
                        NVSD=self.parameters.vertical_salinity_diffusion_coefficient.value_in(units.m**2/units.s),
                        NLTD=self.parameters.lateral_temperature_diffusion_coefficient.value_in(units.m**2/units.s),
                        NVTD=self.parameters.vertical_temperature_diffusion_coefficient.value_in(units.m**2/units.s),
                        EQNSTATE=self.equations_of_state[self.parameters.equation_of_state],
                        RES_BC_FLAG=0, # these can be zero  because corresponding boundary
                        BCFLAG_LNM=0,  # will be used depending on use_interface_.. parameters
                        BCFLAG_TEMP=-1 if self.parameters.use_interface_surface_heat_forcing else 0, 
                        ALP1=self.parameters.time_weight_coefficient_coriolis,
                        ALP2=self.parameters.time_weight_coefficient_bottom_friction,
                        ALP3=self.parameters.time_weight_coefficient_vertical_diffusion,
                        ALP4=self.parameters.time_weight_coefficient_transport,
                        H0=self.parameters.minimum_depth.value_in(units.m),
                        NRAMP=1 if self.parameters.use_ramping else 0,
                        DRAMP=self.parameters.ramping_time.value_in(units.day),
                        HBREAK=self.parameters.hybrid_bottom_friction_hbreak.value_in(units.m),
                        FTHETA=8, # fixed atm
                        FGAMMA=0.33333 # fixed atm
                        )
          NETA=None
          NFLUX=None
          if self.parameters.use_interface_grid:
            NETA=0
            for b in self._elev_boundary: NETA+=len(b)
            NFLUX=0
            for b in self._flow_boundary: 
              if b[0].type in [2,12,22,52]: NFLUX+=len(b)
          param.write(NETA=NETA,NFLUX=NFLUX)
        if self.parameters.use_interface_grid:
          agw=adcirc_grid_writer(coordinates=self.coordinates,nodes=self._nodes,
            elements=self._elements, elevation_boundaries=self._elev_boundary,
            flow_boundaries=self._flow_boundary)
          agw.write_grid()
          if self.mode==self.MODE_3D and self.parameters.baroclinic_flag==True:
              agw.write_density(density_forcing=self.parameters.baroclinic_density_forcing,
               number_of_vertical_levels=NFEN)
        self.overridden().commit_parameters()

    def get_firstlast_node(self):
        return 1,self.get_number_of_nodes()
    def get_firstlast_element(self):
        return 1,self.get_number_of_elements()
    def get_firstlast_node_of_elevation_boundary_segment(self,index_of_segment):
        return 1,self.get_number_of_nodes_in_elevation_boundary_segment(index_of_segment)
    def get_firstlast_node_of_flow_boundary_segment(self,index_of_segment):
        return 1,self.get_number_of_nodes_in_flow_boundary_segment(index_of_segment)
    def get_firstlast_vertical_index(self):
        return 1,self.get_number_of_vertical_nodes()
    def get_firstlast_grid3d(self):
        return 1,self.get_number_of_nodes(),1,self.get_number_of_vertical_nodes()

    def get_node_position_3d(self,index,zindex):
        #~ sigma,z=self.get_node_sigma(index,zindex)
        x,y=self.get_node_position(index)
        return x,y
    def get_node_coordinates_3d(self,index,zindex):
        x,y=self.get_node_coordinates(index)
        return x,y



    def define_parameters(self, object):
        object.add_default_form_parameter(
            "rootdir", 
            "set the root directory", 
            "."
        )

        object.add_default_form_parameter(
            "use_interface_elevation_boundary", 
            "toggle the use of interface boundary conditions for the elevation specified boundaries", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_flow_boundary", 
            "toggle the use of interface boundary conditions for the flow specified boundaries", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_met_forcing", 
            "toggle the use of interface meteorological forcings", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_wave_forcing", 
            "toggle the use of interface wave stress forcings", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_tidal_forcing", 
            "toggle the use of interface tidal forcing", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_surface_heat_forcing", 
            "toggle the use of interface surface heat flux forcing (baroclinic runs)", 
            False
        )
        object.add_interface_parameter(
            "A_H",
            "turbulent lateral friction coefficient",
            100.0 | units.m**2/units.s,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "timestep",
            "ADCIRC timestep",
            360.0 | units.s,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "use_predictor_corrector",
            "flag for use of predictor corrector integrator",
            True,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "calculate_coriolis",
            "flag to let adcirc calculate coriolis frequence (needs spherical coordinates)",
            False,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "use_interface_parameters",
            "flag for use of interface parameters (i.e. write fort.15)",
            False,
            "before_set_interface_parameter"
        ) 
        object.add_interface_parameter(
            "use_interface_grid",
            "flag for use of interface grid (i.e. write fort.14)",
            False,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "central_latitude",
            "central latitude used for CPP projection between x,y and lon,lat",
            0. | units.deg,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "central_longitude",
            "central longitude used for CPP projection between x,y and lon,lat",
            0. | units.deg,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "bottom_friction_law",
            "type of bottom friction law [linear,quadratic, hybrid]",
            "linear",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "GWCE_weighting_factor",
            "factor in GWCE Weighting primitive and wave contributions [-5,-4,-3,-2,-1, or >0 suggested: 0.005-0.1]",
            0.005,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "linear_bottom_friction_coeff",
            "linear bottom friction coefficient [1/s], only used for linear bottom friction law",
            0. | units.s**-1,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "quadratic_bottom_friction_coeff",
            "quadratic bottom friction coefficient [dimensionless], only used for quadratic bottom friction law",
            0.,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "linear_3d_bottom_friction_coeff",
            "linear 3D bottom friction coefficient [0. m/s], only used for linear bottom friction law",
            0. | units.m/units.s,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "quadratic_3d_bottom_friction_coeff",
            "quadratic 3D bottom friction coefficient [dimensionless], only used for quadratic bottom friction law",
            0.,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "bottom_friction_law_3d",
            "type of 3D bottom friction law [noslip, linear, quadratic, loglayer]",
            "noslip",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "begin_time",
            "begin time of the simulation",
            0. | units.day,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "convergence_criterion",
            "absolute convergence criteria (should be no smaller than 500 times the machine precision) [1.e-10] ",
            1.e-10,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "maximum_iterations",
            "maximum number of iterations for iterative solver [25] ",
            25,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "solver_verbosity",
            "verbosity of solver (special for debug) [0-5] ",
            0,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "finite_amplitude_term_parameter",
            "term selecting treatment of finite amplitude terms [0 (linearalize using bathym.) ,1 (actual depth) or 2 (1+wetting/drying), default = 1] ",
            1,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "spatial_derivative_advective_term_parameter",
            "term selecting treatment spatial advective terms [0 (not include) or 1 (include), default = 1] ",
            1,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "time_derivative_advective_term_parameter",
            "term selecting treatment time derivative advective terms [0 (not include) or 1 (include), default = 1] ",
            1,
            "before_set_interface_parameter"
        )
        object.add_method_parameter(
            "get_reference_pressure", 
            None,
            "atmospheric_reference_pressure", 
            "reference atmospheric pressure", 
            None
        )        
        object.add_method_parameter(
            "get_coordinates", 
            None,
            "coordinates", 
            "type of coordinates (spherical or cartesian)", 
            None
        )
        object.add_interface_parameter(
            "baroclinic_flag",
            "flag selecting barotropic or baroclinic model",
            False,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "baroclinic_density_forcing",
            "type of density forcing for baroclinic runs [sigmaT, salinity, temperature or  salinity+temperature]",
            "none",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "number_of_vertical_levels",
            "number of nodes in the vertical grid, only for 3D runs",
            21,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "vertical_grid_type",
            "type of vertical grid spacing: [user, uniform, log, loglinear, doublelog, pgrid, sine]",
            "uniform",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "sigma_levels",
            "dimensionless levels of the vertical grid, from -1 (bottom) to +1 (surface)",
            [],
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "vertical_eddy_viscosity_formulation",
            "form of eddy viscosity used [user, constant, ...]",
            "constant",
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "vertical_eddy_viscosity",
            "(reference) value of eddy viscosity used [m**2/s]",
            0.05 | units.m**2/units.s,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "minimum_vertical_eddy_viscosity",
            "minimal value of eddy viscosity [m**2/s]",
            1.e-6 | units.m**2/units.s,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "vertical_eddy_viscosities",
            "value of eddy viscosity of each depth level [m**2/s]",
            [] | units.m**2/units.s,
            "before_set_interface_parameter"
        )
        
        object.add_interface_parameter(
            "lateral_salinity_diffusion_coefficient",
            "value of the lateral salinity diffusion coefficient for baroclinic runs [m**2/s]",
            0. | units.m**2/units.s,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "vertical_salinity_diffusion_coefficient",
            "value of the vertical salinity diffusion coefficient for baroclinic runs  [m**2/s]",
            0. | units.m**2/units.s,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "lateral_temperature_diffusion_coefficient",
            "value of the lateral temperature diffusion coefficient for baroclinic runs  [m**2/s]",
            0. | units.m**2/units.s,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "vertical_temperature_diffusion_coefficient",
            "value of the vertical temperature diffusion coefficient for baroclinic runs  [m**2/s]",
            0. | units.m**2/units.s,
            "before_set_interface_parameter"
        )

        object.add_interface_parameter(
            "equation_of_state",
            "equation of state rho(T,sal) to use for baroclinic runs [CushmanRoisin, McDougall, unesco]",
            "CushmanRoisin",
            "before_set_interface_parameter"
        )
        
        object.add_interface_parameter(
            "time_weight_coefficient_coriolis",
            "time weighting for coriolis terms in 3D velocity eq. (0.=fully explicit, 1. = fully implicit)",
            0.5,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "time_weight_coefficient_bottom_friction",
            "time weighting for bottom friction terms in 3D velocity eq. (0.=fully explicit, 1. = fully implicit)",
            0.5,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "time_weight_coefficient_vertical_diffusion",
            "time weighting for vertical diffusion terms in 3D velocity eq. (0.=fully explicit, 1. = fully implicit)",
            0.5,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "time_weight_coefficient_transport",
            "time weighting for transport eq. terms in baroclinic sims (0.=fully explicit, 1. = fully implicit)",
            0.5,
            "before_set_interface_parameter"
        )
        object.add_default_form_parameter(
            "use_interface_lnm_boundary", 
            "toggle the use of interface boundary conditions for the level of no motion (3D baroclinic runs)", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_salinity_boundary", 
            "toggle the use of interface boundary conditions for salinity (3D baroclinic runs)", 
            False
        )
        object.add_default_form_parameter(
            "use_interface_temperature_boundary", 
            "toggle the use of interface boundary conditions for temperature (3D baroclinic runs)", 
            False
        )
        object.add_interface_parameter(
            "use_ramping",
            "whether to use ADCIRC internal (hyperbolic) ramping for elev and flux boundary, tidal potential, met and wave stress forcing.",
            False,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "ramping_time",
            "ramping time scale",
            12 | units.hour,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "minimum_depth",
            "minimum depth or threshold water level for wetting/drying",
            1. | units.m,
            "before_set_interface_parameter"
        )
        object.add_interface_parameter(
            "hybrid_bottom_friction_hbreak",
            "transition depth for hybrid bottom friction law ",
            50 | units.m,
            "before_set_interface_parameter"
        )
        
    def define_properties(self, object):
        object.add_property('get_model_time', public_name = "model_time")


    def get_grid(self):
        return StaggeredGrid(self.elements, self.nodes)


    def define_particle_sets(self, object):
        axes_names = ['x', 'y']
        if self.coordinates == "spherical":
            axes_names = ['lon','lat']
        object.define_grid('nodes', axes_names = axes_names, grid_class=datamodel.UnstructuredGrid)
        object.set_grid_range('nodes', 'get_firstlast_node')
        object.add_getter('nodes', 'get_node_position', names=('x','y'))
        object.add_getter('nodes', 'get_node_coordinates', names=('lon','lat'))          
        object.add_getter('nodes', 'get_node_depth', names=('depth',))
        object.add_getter('nodes', 'get_node_eta', names=('eta',))
        object.add_getter('nodes', 'get_node_vx', names=('vx',))
        object.add_getter('nodes', 'get_node_vy', names=('vy',))
        object.add_getter('nodes', 'get_node_deta_dt', names=('deta_dt',))
        object.add_getter('nodes', 'get_node_status', names=('status',))
        object.add_setter('nodes', 'set_node_eta', names=('eta',))
        object.add_setter('nodes', 'set_node_vx', names=('vx',))
        object.add_setter('nodes', 'set_node_vy', names=('vy',))
        object.add_setter('nodes', 'set_node_deta_dt', names=('deta_dt',))
        object.add_setter('nodes', 'set_node_status', names=('status',))

                
        if self.mode in [self.MODE_3D]:
            object.define_grid('grid3d',axes_names=['x','y','z'])
            object.set_grid_range('grid3d', 'get_firstlast_grid3d')
            object.add_getter('grid3d', 'get_node_sigma', names = ('sigma','z'))
            object.add_getter('grid3d', 'get_node_velocities_3d', names = ('wx','wy','wz'))
            object.add_getter('grid3d', 'get_node_temperature_3d', names = ('temperature',))
            object.add_getter('grid3d', 'get_node_position_3d', names = ('x','y'))
            object.add_getter('grid3d', 'get_node_coordinates_3d', names = ('lon','lat'))

            object.add_gridded_getter('nodes', 'get_node_sigma','get_firstlast_vertical_index', names = ('sigma','z'))
            object.add_gridded_getter('nodes', 'get_node_velocities_3d','get_firstlast_vertical_index', names = ('wx','wy','wz'))
            object.add_gridded_getter('nodes', 'get_node_temperature_3d','get_firstlast_vertical_index', names = ('temperature',))

        object.define_grid('forcings',axes_names = axes_names, grid_class=datamodel.UnstructuredGrid)
        object.set_grid_range('forcings', 'get_firstlast_node')
        object.add_getter('forcings', 'get_node_coriolis_f', names=('coriolis_f',))
        object.add_setter('forcings', 'set_node_coriolis_f', names=('coriolis_f',))
        object.add_getter('forcings', 'get_node_wind_stress', names=('tau_x','tau_y'))
        object.add_setter('forcings', 'set_node_wind_stress', names=('tau_x','tau_y'))
        object.add_getter('forcings', 'get_node_wave_stress', names=('wave_tau_x','wave_tau_y'))
        object.add_setter('forcings', 'set_node_wave_stress', names=('wave_tau_x','wave_tau_y'))
        object.add_getter('forcings', 'get_node_atmospheric_pressure', names=('pressure',))
        object.add_setter('forcings', 'set_node_atmospheric_pressure', names=('pressure',))
        object.add_getter('forcings', 'get_node_tidal_potential', names=('tidal_potential',))
        object.add_setter('forcings', 'set_node_tidal_potential', names=('tidal_potential',))
        object.add_getter('forcings', 'get_node_surface_heat_flux', names=('surface_heat_flux',))
        object.add_setter('forcings', 'set_node_surface_heat_flux', names=('surface_heat_flux',))
        object.add_getter('forcings', 'get_node_position', names=('x','y'))
        object.add_getter('forcings', 'get_node_coordinates', names=('lon','lat'))


        object.define_grid('elements',axes_names = axes_names, grid_class=datamodel.UnstructuredGrid)
        object.set_grid_range('elements', 'get_firstlast_element')    
        object.add_getter('elements', 'get_element_nodes', names=('n1','n2','n3'))
        object.add_getter('elements', 'get_element_position', names=('x','y'))
        object.add_getter('elements', 'get_element_coordinates', names=('lon','lat'))
        object.add_getter('elements', 'get_element_status', names=('status',))
        object.add_setter('elements', 'set_element_status', names=('status',))

          
    def elevation_boundaries(self):
        n=self.get_number_of_elevation_boundary_segments()
        for i in range(1,n+1):
            yield self._create_new_grid(self.specify_elevation_boundary_grid, index = i)

    def specify_elevation_boundary_grid(self, definition, index=1):
        definition.set_grid_range('get_firstlast_node_of_elevation_boundary_segment') 
        definition.add_getter('get_elevation_boundary_node', names=('node',))
        definition.add_getter('get_elevation_boundary_eta', names=('eta',))
        definition.add_setter('set_elevation_boundary_eta', names=('eta',))

        if self.parameters.use_interface_lnm_boundary:
            definition.add_getter('get_boundary_lnm', names=('lnm',))
            definition.add_setter('set_boundary_lnm', names=('lnm',))

        if self.parameters.use_interface_salinity_boundary:
            definition.add_gridded_getter('get_boundary_salinity', 'get_firstlast_vertical_index',names=('salinity',))
            definition.add_gridded_setter('set_boundary_salinity', 'get_firstlast_vertical_index',names=('salinity',))
        if self.parameters.use_interface_temperature_boundary:
            definition.add_gridded_getter('get_boundary_temperature', 'get_firstlast_vertical_index',names=('temperature',))
            definition.add_gridded_setter('set_boundary_temperature', 'get_firstlast_vertical_index',names=('temperature',))

        definition.define_extra_keywords({'index_of_segment':index})


    def flow_boundaries(self):
        n=self.get_number_of_flow_boundary_segments()
        for i in range(1,n+1):
            yield self._create_new_grid(self.specify_flow_boundary_grid, index = i)

    def specify_flow_boundary_grid(self, definition, index=1):
        definition.set_grid_range('get_firstlast_node_of_flow_boundary_segment') 
        definition.add_getter('get_flow_boundary_node', names=('node',))
        definition.add_getter('get_flow_boundary_type', names=('type',))

        definition.add_getter('get_flow_boundary_fluxx', names=('flux_x',))
        definition.add_setter('set_flow_boundary_fluxx', names=('flux_x',))
        definition.add_getter('get_flow_boundary_fluxy', names=('flux_y',))
        definition.add_setter('set_flow_boundary_fluxy', names=('flux_y',))

        definition.define_extra_keywords({'index_of_segment':index})

    def define_state(self, object):
        object.set_initial_state('UNINITIALIZED')
        object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        object.add_method('!UNINITIALIZED', 'before_get_parameter')
        object.add_method('!UNINITIALIZED', 'before_set_parameter')
        object.add_method('END', 'before_get_parameter')
        object.add_transition('!UNINITIALIZED!STOPPED', 'END', 'cleanup_code')
        object.add_transition('END', 'STOPPED', 'stop', False)
        object.add_method('STOPPED', 'stop')

        object.add_transition('INITIALIZED','EDIT','commit_parameters')

        #~ object.set_initial_state('UNINITIALIZED')
        #~ object.add_transition('!STOPPED', 'END', 'cleanup_code')
        #~ object.add_transition('UNINITIALIZED', 'INITIALIZED', 'initialize_code')
        #~ object.add_transition('END', 'STOPPED', 'stop', False)
        #~ object.add_method('STOPPED', 'stop')        
 
        object.add_method('INITIALIZED', 'assign_grid_and_boundary')
        object.add_method('INITIALIZED', 'before_set_interface_parameter')
        object.add_method('INITIALIZED', 'set_rootdir')

        for state in ["RUN","EDIT","EVOLVED"]:
          object.add_method(state,"get_model_time")
          object.add_method(state,"get_timestep")          
          object.add_method(state,"get_number_of_nodes")
          object.add_method(state,"get_number_of_elements")
          object.add_method(state,"get_number_of_nodes_in_elevation_boundary_segment")
          object.add_method(state,"get_number_of_nodes_flow_boundary_segment")
          object.add_method(state,"get_number_of_vertical_nodes")
          object.add_method(state,"get_grid")

        object.add_transition('EDIT', 'RUN', 'commit_grid')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_method('EVOLVED', 'evolve_model')
