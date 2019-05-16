import numpy

from omuse.units import units
from amuse.units import trigo

from amuse.datamodel import Grid

default_parameters={
 "RUNDES": dict(dtype=str, value="AMUSE", description="32 CHARACTER ALPHANUMERIC RUN DESCRIPTION"),
 "RUNID": dict(dtype=str,value="AMUSE test",description="24 CHARACTER ALPANUMERIC RUN IDENTIFICATION"),
 "NFOVER": dict(dtype=int , value=0, description="NONFATAL ERROR OVERRIDE OPTION"),
 "NABOUT": dict(dtype=int , value=0, description="ABREVIATED OUTPUT OPTION PARAMETER"),
 "NSCREEN": dict(dtype=int , value=1, description="UNIT 6 OUTPUT OPTION PARAMETER"),
 "IHOT": dict(dtype=int , value=0, description="HOT START PARAMETER"),
 "ICS": dict(dtype=int , value=1, description="COORDINATE SYSTEM SELECTION PARAMETER"),
 "IM": dict(dtype=int , value=0, description="MODEL SELECTION PARAMETER "),
 "NOLIBF": dict(dtype=int , value=1, description="BOTTOM FRICTION TERM SELECTION PARAMETER"),
 "NOLIFA": dict(dtype=int , value=1, description="FINITE AMPLITUDE TERM SELECTION PARAMETER"),
 "NOLICA": dict(dtype=int , value=1, description="SPATIAL DERIVATIVE CONVECTIVE SELECTION PARAMETER"),
 "NOLICAT": dict(dtype=int , value=1, description="TIME DERIVATIVE CONVECTIVE TERM SELECTION PARAMETER"),
 "NWP": dict(dtype=int , value=0, description="VARIABLE BOTTOM FRICTION AND LATERAL VISCOSITY OPTION PARAMETER"),
 "NCOR": dict(dtype=int , value=0, description="VARIABLE CORIOLIS IN SPACE OPTION PARAMETER"),
 "NTIP": dict(dtype=int , value=0, description="TIDAL POTENTIAL OPTION PARAMETER"),
 "NWS": dict(dtype=int , value=0, description="WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER"),
 "NRAMP": dict(dtype=int , value=0, description="RAMP FUNCTION OPTION"),
 "G": dict(dtype=float , value=9.81, description="ACCELERATION DUE TO GRAVITY: DETERMINES UNITS"),
 "TAU0": dict(dtype=float , value=0.005, description="WEIGHTING FACTOR IN GWCE"),
 "DTDP": dict(dtype=float , value=-360., description="TIME STEP (IN SECONDS)"),
 "STATIM": dict(dtype=float , value=0.00, description="STARTING TIME (IN DAYS)"),
 "REFTIM": dict(dtype=float , value=0.00, description="REFERENCE TIME (IN DAYS)"),
 "RNDAY": dict(dtype=float , value=100.0, description="TOTAL LENGTH OF SIMULATION (IN DAYS)"),
 "DRAMP": dict(dtype=float , value=2.0, description="DURATION OF RAMP FUNCTION (IN DAYS)"),
 "A00": dict(dtype=float, value=0.35, description="TIME WEIGHTING FACTOR FOR THE GWCE EQUATION"),
 "B00": dict(dtype=float, value=0.30, description="TIME WEIGHTING FACTOR FOR THE GWCE EQUATION"),
 "C00": dict(dtype=float, value=0.35, description="TIME WEIGHTING FACTOR FOR THE GWCE EQUATION"),
 "H0": dict(dtype=float , value=1.0, description="MIMINMUM CUTOFF DEPTH"),
 "SLAM0": dict(dtype=float , value=0.0, description="CENTER OF CPP PROJECTION (NOT USED IF ICS=1, NTIP=0, NCOR=0)"),
 "SFEA0": dict(dtype=float , value=0.0, description="CENTER OF CPP PROJECTION (NOT USED IF ICS=1, NTIP=0, NCOR=0)"),
 "TAU": dict(dtype=float , value=0.0, description="HOMOGENEOUS LINEAR  BOTTOM FRICTION COEFFICIENT"),
 "CF": dict(dtype=float , value=0.00, description="HOMOGENEOUS QUADRATIC BOTTOM FRICTION COEFFICIENT"),
 "ESLM": dict(dtype=float , value=1000.0, description="LATERAL EDDY VISCOSITY COEFFICIENT; IGNORED IF NWP =1"),
 "CORI": dict(dtype=float , value=0.0, description="CORIOLIS PARAMETER: IGNORED IF NCOR = 1"),
 "NTIF": dict(dtype=int , value=0, description="TOTAL NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED"),
 "NBFR": dict(dtype=int , value=-1, description="TOTAL NUMBER OF FORCING FREQUENCIES ON OPEN BOUNDARIES"),
 "NFFR": dict(dtype=int , value=-2, description="TOTAL NUMBER OF FORCING FREQUENCIES ON FLOW BOUNDARIES"),
 "ANGINN": dict(dtype=float , value=110.0, description="INNER ANGLE THRESHOLD"),
 "ITITER": dict(dtype=int, value=1, description="ALGEBRAIC SOLUTION PARAMETER: iterative or lumped"),
 "ISLDIA": dict(dtype=int, value=0, description="ALGEBRAIC SOLUTION PARAMETER: solver verbosity"),
 "CONVCR": dict(dtype=int, value=1.e-10, description="ALGEBRAIC SOLUTION PARAMETER: convergence criterion"),
 "ITMAX": dict(dtype=int, value=25, description="ALGEBRAIC SOLUTION PARAMETER: max iterations"),
 "ISLIP": dict(dtype=int, value=1, description="slip code"),
 "KP": dict(dtype=float, value=0.01, description="slip coefficient"),
 "Z0S": dict(dtype=float, value=0.01, description="free surface roughness (const horiz)"),
 "Z0B": dict(dtype=float, value=0.1, description="bottom roughness (const horiz)"),
 "ALP1": dict(dtype=float, value=0.5, description="timestepping coefficient (coriolis)"),
 "ALP2": dict(dtype=float, value=0.5, description="timestepping coefficient (bottom fric.)"),
 "ALP3": dict(dtype=float, value=0.5, description="timestepping coefficient (vert. diffusion)"),
 "IGC": dict(dtype=int, value=1, description="type of vertical grid"),
 "NFEN": dict(dtype=int, value=21, description="number of nodes in vertical grid"),
 "IEVC": dict(dtype=int, value=1, description="vertical eddy viscosity code"),
 "EVMIN": dict(dtype=float, value=1.e-6, description="minimal vertical eddy viscosity"),
 "EVCON": dict(dtype=float, value=0.05, description="vertical eddy viscosity constant"),
 "HBREAK": dict(dtype=float, value=1., description="hybrid bottom friction break depth"),
 "FTHETA": dict(dtype=float, value=10., description="dimensionless shallow/deep water bottom friction law transition parameter"),
 "FGAMMA": dict(dtype=float, value=0.33333333, description="friction factor parameter"),
 "RES_BC_FLAG" : dict(dtype=int, value=0, descrption="Controls the type of boundary conditions used in the 3D baroclinic simulations (=IDEN)"),
 "BCFLAG_LNM" : dict(dtype=int, value=0, description="Flag to control the levels of no motion parameterization used in the 3D ADCIRC run."),
 "BCFLAG_TEMP" : dict(dtype=int, value=0, description="surface heat flux parameterization (must be 0 or 1)"),
 "RBCTIMEINC" : dict(dtype=float, value=-1, description=" Time interval between data sets for the level of no motion boundary condition data, in seconds"),
 "SBCTIMEINC" : dict(dtype=float, value=-1, description=" Time interval between data sets for the level of salinity condition data, in seconds"),
 "TBCTIMEINC" : dict(dtype=float, value=-1, description=" Time interval between data sets for the level of temperature condition data, in seconds"),
 "TTBCTIMEINC" : dict(dtype=float, value=-1, description=" Time interval between data sets for the surface heat flux condition data, in seconds"),
 "BCSTATIM" : dict(dtype=float, value=0., desciption="Starting time (in seconds since ADCIRC cold start) for boundary condition data for the level of no motion boundary condition."),
 "SBCSTATIM" : dict(dtype=float, value=0., desciption="Starting time (in seconds since ADCIRC cold start) for boundary condition data for the salinity boundary condition."),
 "TBCSTATIM" : dict(dtype=float, value=0., desciption="Starting time (in seconds since ADCIRC cold start) for boundary condition data for the temperature boundary condition."),
 "TTBCSTATIM" : dict(dtype=float, value=0., desciption="Starting time (in seconds since ADCIRC cold start) for boundary condition data for the surface heat flux boundary condition."),
 "SPONGEDIST" : dict(dtype=float, value=0., description="spatial ramp (in m) on the wind and advection terms"),
 "EQNSTATE" : dict(dtype=int, value=1, description="equation of state to use (1,2 or 3)"),
 "NLSD" : dict(dtype=float, value=100., description="lateral salinity diffusion coefficient"),
 "NLTD" : dict(dtype=float, value=100., description="lateral temperature diffusion coefficient"),
 "NVSD" : dict(dtype=float, value=1.e-4, description="vertical salinity diffusion coefficient"),
 "NVTD" : dict(dtype=float, value=1.e-4, description="vertical temperature diffusion coefficient"),
 "ALP4" : dict(dtype=float, value=0.5, description="timestepping coefficient for transport eq. terms"),
 "NTF" : dict(dtype=int, value=0, description="temperature boundary condition file type"),
 "IDEN" : dict(dtype=int, value=0, description="type of density forcing for baroclinic runs"),
 }
# THETA1, THETA2

def get_default_parameter_set():
  return dict([(key,val["value"]) for key,val in default_parameters.items()])

outputblock1="""\
 0 0.0 5.0  3                        ! NOUTE,TOUTSE,TOUTFE,NSPOOLE:ELEV STATION OUTPUT INFO (UNIT  61)
 0                                   ! TOTAL NUMBER OF ELEVATION RECORDING STATIONS
 0 0.0 5.0  3                        ! NOUTV,TOUTSV,TOUTFV,NSPOOLV:VEL STATION OUTPUT INFO (UNIT  62)
 0                                   ! TOTAL NUMBER OF VELOCITY RECORDING STATIONS
 0 0.0 5.0 3                         ! NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : GLOBAL ELEVATION OUTPUT INFO (UNIT  63)
 0 0.0 5.0 3                         ! NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV : GLOBAL VELOCITY  OUTPUT INFO (UNIT  64)    
 0                                   ! NHARFR - NUMBER OF CONSTITUENTS TO BE INCLUDED IN THE HARMONIC ANALYSIS
 4.00  5.00  1  0.0                  ! THAS,THAF,NHAINC,FMV - HARMONIC ANALYSIS PARAMETERS         
 0 0 0 0                             ! NHASE,NHASV,NHAGE,NHAGV - CONTROL HARMONIC ANALYSIS AND OUTPUT TO UNITS 51,52,53,54
 0 1236                               ! NHSTAR,NHSINC - HOT START FILE GENERATION PARAMETERS                  
"""

outputblock2="""\
 0  0.0  5.0  3                      ! DTS station output
 0
 0  0.0  5.0  3                      ! velocity station output
 0
 0  0.0  5.0  3                      ! turbulence station output
 0
 0   0.0  5.0     3                  ! DTS global output
 0   0.0  5.0     3                  ! velocity global output
 0   0.0  5.0     3                  ! turbulence global output
"""

class adcirc_file_writer(file):
  def write_var(self,*var):
    self.write(' '.join([str(x) for x in var])+"\n")
  def write_var_rows(self,*var):
    for v in zip(var):
      self.write_var(*v)


class adcirc_parameter_writer(object):
  def __init__(self,filename="fort.15"):
    self.filename=filename
    self.parameters=get_default_parameter_set()
  def write(self, NETA=0, NFLUX=None):
    if NFLUX is None: NFLUX=0
    if NETA is None: NETA=0
    with adcirc_file_writer(self.filename,'w') as f:
      param=self.parameters
      f.write_var(param["RUNDES"])    
      f.write_var(param["RUNID"])    
      f.write_var(param["NFOVER"])
      f.write_var(param["NABOUT"])
      f.write_var(param["NSCREEN"])
      f.write_var(param["IHOT"])
      f.write_var(param["ICS"])
      f.write_var(param["IM"])
      if param["IM"] in [21,31]:
        f.write_var(param["IDEN"])
      else:
        pass
      f.write_var(param["NOLIBF"])
      f.write_var(param["NOLIFA"])
      f.write_var(param["NOLICA"])
      f.write_var(param["NOLICAT"])
      f.write_var(param["NWP"])
      if param["NWP"]>0:
        for x in param["AttrName"]:
          f.write_var(x)
      f.write_var(param["NCOR"])
      f.write_var(param["NTIP"])
      f.write_var(param["NWS"])
      f.write_var(param["NRAMP"])
      f.write_var(param["G"])    
      f.write_var(param["TAU0"])
      if param["TAU0"]==-5.0:
        f.write_var(param["Tau0FullDomainMin"])
        f.write_var(param["Tau0FullDomainMax"])
      else:
        pass
      f.write_var(param["DTDP"])
      f.write_var(param["STATIM"])
      f.write_var(param["REFTIM"])
      if param["NWS"] not in [0]:
        raise Exception("NWS not equal zero, check..")
  #~ WTIMINC  Supplemental Meteorological/Wave/Ice Parameters Line
      f.write_var(param["RNDAY"])
      if param["NRAMP"] in [0,1]:
        f.write_var(param["DRAMP"])
      elif param["NRAMP"] in [2,3,4,5,6,7,8]:
        raise Exception("DRAMP for NRAMP>1 tbd")
      else:
        pass
      f.write_var(param["A00"],param["B00"],param["C00"])
      if param["NOLIFA"] in [0,1]: 
        f.write_var(param["H0"])
      elif param["NOLIFA"] in [2,3]:
        f.write_var(param["H0"],0,0,param["VELMIN"])
      f.write_var(param["SLAM0"],param["SFEA0"])
      if param["NOLIBF"]==0:
        f.write_var(param["TAU"])
      elif param["NOLIBF"]==1:
        f.write_var(param["CF"])
      elif param["NOLIBF"]==2:
        f.write_var(param["CF"],param["HBREAK"],param["FTHETA"],param["FGAMMA"])
      if param["IM"]!=10:
        f.write_var(param["ESLM"])
      elif param["IM"]==10:
        f.write_var(param["ESLM"],param["ESLC"])
      f.write_var(param["CORI"])
      f.write_var(param["NTIF"])
      if param["NTIF"]>0:
        for tag,tpk,amigt,etrf,fft,facet in zip(param["TIPOTAG"],param["TPK"],param["AMIGT"],param["ETRF"],param["FFT"],param["FACET"]):
          f.write_var(tag)
          f.write_var(tpk,amigt,etrf,fft,facet)
      f.write_var(param["NBFR"])
      if param["NBFR"]>0:
        for tag,amig,ff,face in zip(param["BOUNTAG"],param["AMIG"],param["FF"],param["FACE"]):
          f.write_var(tag)
          f.write_var(amig,ff,face)
        for tag,data in param["BOUNDARY_FORCING_DATA"].items():
            f.write_var(tag)
            for emo,efa in zip(*data):
                f.write_var(emo,efa)      
      f.write_var(param["ANGINN"])
      if NFLUX>0:
        f.write_var(param["NFFR"])
      if param["NFFR"]>0:
        for tag,amig,ff,face in zip(param["FBOUNTAG"],param["FAMIGT"],param["FFF"],param["FFACE"]):
          f.write_var(tag)
          f.write_var(amig,ff,face)
        for tag,data in param["FLUX_BOUNDARY_FORCING_DATA"].items():
            f.write_var(tag)
            for qnam,qnph in zip(*data):
                f.write_var(qnam,qnph)
            #~ for qnam,qnph,enam,enph in zip(*data):
                #~ f.write_var(qnam,qnph,enam,enph)
      f.write(outputblock1)
      f.write_var(param["ITITER"],param["ISLDIA"],
                  param["CONVCR"],param["ITMAX"])
      if param['IM'] in [1,11,21,31,2]:
        # continue writing 3D info
        f.write_var(param["IDEN"])
        f.write_var(param['ISLIP'],param['KP'])
        f.write_var(param['Z0S'],param['Z0B'])
        f.write_var(param['ALP1'],param['ALP2'],param['ALP3'])
        f.write_var(param['IGC'],param['NFEN'])
        if  param['IGC']==0:
          assert  len(param['SIGMA'])==param['NFEN']
          f.write_var_rows(param['SIGMA'])
        f.write_var(param['IEVC'],param['EVMIN'],param['EVCON'])
        if param['IEVC'] in [50,51]:
          f.write_var(param['THETA1'],param['THETA2'])
        if  param['IEVC']==0:
          assert  len(param['EVTOT'])==param['NFEN']
          f.write_var_rows(param['EVTOT'])
        f.write(outputblock2)
        if param['IDEN'] != 0: # param['IM'] in [21,31]:
          if param['RES_BC_FLAG'] not in [0,param['IDEN']]:
            print param['RES_BC_FLAG'], param['IDEN']
            raise Exception("RES_BC_FLAG should be equal to IDEN or 0")
          f.write_var(param['RES_BC_FLAG'],param['BCFLAG_LNM'],param['BCFLAG_TEMP'])
          if param['RES_BC_FLAG'] not in range(-4,5):
            raise Exception("unexpected RES_BC_FLAG value")
          if param['RES_BC_FLAG']<0:
            if abs(param['RES_BC_FLAG'])>=1 and NETA>0:
              f.write_var(param['RBCTIMEINC'])
              f.write_var(param['BCSTATIM'])
          elif param['RES_BC_FLAG']>0:
            if abs(param['RES_BC_FLAG'])==1 and NETA>0:
              f.write_var(param['RBCTIMEINC'])
              f.write_var(param['BCSTATIM'])
            elif abs(param['RES_BC_FLAG'])==2 and NETA>0:
              f.write_var(param['RBCTIMEINC'],param['SBCTIMEINC'])
              f.write_var(param['BCSTATIM'],param['SBCSTATIM'])
            elif abs(param['RES_BC_FLAG'])==3 and NETA>0:
              f.write_var(param['RBCTIMEINC'],param['TBCTIMEINC'])
              f.write_var(param['BCSTATIM'],param['TBCSTATIM'])
              if param['BCFLAG_TEMP']!=0:
                  f.write_var(param['TTBCTIMEINC'],param['TTBCSTATIM'])
            elif abs(param['RES_BC_FLAG'])==4 and NETA>0:
              f.write_var(param['RBCTIMEINC'],param['SBCTIMEINC'],param['TBCTIMEINC'])
              f.write_var(param['BCSTATIM'],param['SBCSTATIM'],param['TBCSTATIM'])
              if param['BCFLAG_TEMP']!=0:
                  f.write_var(param['TTBCTIMEINC'],param['TTBCSTATIM'])
          f.write_var(param['SPONGEDIST'])
          f.write_var(param['EQNSTATE'])
        if param['IDEN']>0:
          f.write_var(param['NLSD'],param['NVSD'])
          f.write_var(param['NLTD'],param['NVTD'])
          f.write_var(param['ALP4'])
        #~ if param['IDEN'] in [3,4]:
          #~ f.write_var(param['NTF'])



      
  def update(self, **kwargs):
    self.parameters.update(kwargs)

class adcirc_grid_writer(object):
  
  def __init__(self,filename="fort.14",coordinates="cartesian", density_filename="fort.11",
                nodes=None,elements=None, elevation_boundaries=None, flow_boundaries=None):
    self.nodes=nodes
    self.elements=elements
    self.elevation_boundaries=elevation_boundaries
    self.flow_boundaries=flow_boundaries
    self.filename=filename
    self.density_filename=density_filename
    self.coordinates=coordinates
    if coordinates not in ["cartesian","spherical"]:
      raise Exception("coordinates must be cartesian or spherical")  

  def write(self,x,y,depth,elements,elev_boundary,flow_boundary):
    f=adcirc_file_writer(self.filename,'w')
    f.write_var("AMUSE grid")
    f.write_var(len(elements),len(x))
    for i,x_,y_,d_ in zip(range(1,len(x)+1),x,y,depth):
      f.write_var(i,x_,y_,d_)
    for i,n in enumerate(elements):
      f.write_var(i+1,3,n[0],n[1],n[2])

    f.write_var(len(elev_boundary))
    NETA=sum(map(lambda x:len(x[0]), elev_boundary))
    f.write_var(NETA)
    for b,t in elev_boundary:
      f.write_var(len(b), t)
      for node in b: f.write_var(node)    
    
    f.write_var(len(flow_boundary))
    NVEL=sum(map(lambda x:len(x[0]), flow_boundary))
    f.write_var(NVEL)
    for b,t in flow_boundary:
      f.write_var(len(b),t)
      for node in b: f.write_var(node)    
      
    f.close()

  def write_density(self, nodes=None, density_forcing=None, number_of_vertical_levels=None,
      Tconst=10. | units.Celsius):
    if nodes is None: nodes=self.nodes
    if density_forcing!="temperature":
      raise Exception("not implemented yet")
    else:
      print "assumes constant temperature:", Tconst

    NP=len(nodes)
    NFEN=number_of_vertical_levels
    f=adcirc_file_writer(self.density_filename,'w')
    f.write_var("header 1")
    f.write_var("header 2")
    f.write_var(NFEN,NP)
    for i in range(NP):
      for j in range(NFEN):
        f.write_var(i+1,j+1,Tconst.value_in(units.Celsius))
      
  def write_grid(self, nodes=None, elements=None, elevation_boundaries=None,flow_boundaries=None):
    if nodes is None: nodes=self.nodes
    if elements is None: elements=self.elements
    if elevation_boundaries is None: elevation_boundaries=self.elevation_boundaries
    if flow_boundaries is None: flow_boundaries=self.flow_boundaries
    
    if self.coordinates=="cartesian":
      x=nodes.x.value_in(units.m)
      y=nodes.y.value_in(units.m)
    else:
      x=nodes.lon.value_in(units.deg)
      y=nodes.lat.value_in(units.deg)      
    depth=nodes.depth.value_in(units.m)
    try:
        element_nodes=elements.nodes+1
    except:
        element_nodes=numpy.column_stack((elements.n1,elements.n2,elements.n3))+1
    elev_boundary_nodes=[(b.nodes+1,0) for b in elevation_boundaries] # type = always zero atm
    flow_boundary_nodes=[(b.nodes+1,b[0].type) for b in flow_boundaries]
    self.write(x,y,depth,element_nodes,elev_boundary_nodes,flow_boundary_nodes)

if __name__=="__main__":
  from simple_triangulations import square_domain_sets
  
  nodes,elements,boundary=square_domain_sets(N=10)
  nodes.depth=4. | units.km

  adcirc_grid_writer().write_grid(nodes,elements,[(boundary,0)])
  
  adcirc_parameter_writer("test").write()
