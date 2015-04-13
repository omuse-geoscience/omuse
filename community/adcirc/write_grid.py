import numpy

from amuse.units import units

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
 "NRAMP": dict(dtype=int , value=1, description="RAMP FUNCTION OPTION"),
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
 "TAU": dict(dtype=float , value=0.00, description="HOMOGENEOUS LINEAR  BOTTOM FRICTION COEFFICIENT"),
 "CF": dict(dtype=float , value=0.00, description="HOMOGENEOUS QUADRATIC BOTTOM FRICTION COEFFICIENT"),
 "ESLM": dict(dtype=float , value=1000.0, description="LATERAL EDDY VISCOSITY COEFFICIENT; IGNORED IF NWP =1"),
 "CORI": dict(dtype=float , value=0.0, description="CORIOLIS PARAMETER: IGNORED IF NCOR = 1"),
 "NTIF": dict(dtype=int , value=0, description="TOTAL NUMBER OF TIDAL POTENTIAL CONSTITUENTS BEING FORCED"),
 "NBFR": dict(dtype=int , value=0, description="TOTAL NUMBER OF FORCING FREQUENCIES ON OPEN BOUNDARIES"),
 "ANGINN": dict(dtype=float , value=110.0, description="INNER ANGLE THRESHOLD"),
 "ITITER": dict(dtype=int, value=1, description="ALGEBRAIC SOLUTION PARAMETER: iterative or lumped"),
 "ISLDIA": dict(dtype=int, value=0, description="ALGEBRAIC SOLUTION PARAMETER: solver verbosity"),
 "CONVCR": dict(dtype=int, value=1.e-10, description="ALGEBRAIC SOLUTION PARAMETER: convergence criterion"),
 "ITMAX": dict(dtype=int, value=25, description="ALGEBRAIC SOLUTION PARAMETER: max iterations"),
 }

footer="""\
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

class adcirc_file_writer(file):
  def write_var(self,*var):
    self.write(' '.join([str(x) for x in var])+"\n")

class adcirc_parameter_writer(object):
  def __init__(self,filename="fort.15"):
    self.filename=filename
    self.parameters=dict([(key,val["value"]) for key,val in default_parameters.items()])
  def write(self):
    f=adcirc_file_writer(self.filename,'w')
    param=self.parameters
    f.write_var(param["RUNDES"])    
    f.write_var(param["RUNID"])    
    f.write_var(param["NFOVER"])
    f.write_var(param["NABOUT"])
    f.write_var(param["NSCREEN"])
    f.write_var(param["IHOT"])
    f.write_var(param["ICS"])
    f.write_var(param["IM"])
    if param["IM"] in [20,30]:
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
#~ WTIMINC  Supplemental Meteorological/Wave/Ice Parameters Line
    f.write_var(param["RNDAY"])
    if param["NRAMP"] in [0,1]:
      f.write_var(param["DRAMP"])
    elif param["NRAMP"] in [2,3,4,5,6,7,8]:
      raise Exception("tbd")
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
    if param["IM"] in [0,1,2]:
      f.write_var(param["ESLM"])
    elif param["IM"]==10:
      f.write_var(param["ESLM"],param["ESLC"])
    f.write_var(param["CORI"])
    f.write_var(param["NTIF"])
    f.write_var(param["NBFR"])
    f.write_var(param["ANGINN"])
    f.write(footer)
    f.write_var(param["ITITER"],param["ISLDIA"],
                param["CONVCR"],param["ITMAX"])
    f.close()
  def set_non_default(self,A_H=100. | units.m**2/units.s, timestep=360. | units.s, 
         use_predictor_corrector=True):
    self.parameters["ESLM"]=A_H.value_in(units.m**2/units.s)
    if use_predictor_corrector:
      self.parameters["DTDP"]=-timestep.value_in(units.s)
    else:
      self.parameters["DTDP"]=timestep.value_in(units.s)
      

class adcirc_grid_writer(object):
  
  def __init__(self,filename="fort.14"):
    self.filename=filename
    self.unit_length = units.m
    self.unit_position = units.m

  def write(self,x,y,depth,elements,boundary):
    f=adcirc_file_writer(self.filename,'w')
    f.write_var("AMUSE grid")
    f.write_var(len(elements),len(x))
    for i,x_,y_,d_ in zip(range(1,len(x)+1),x,y,depth):
      f.write_var(i,x_,y_,d_)
    for i,n in enumerate(elements):
      f.write_var(i+1,3,n[0],n[1],n[2])
    f.write_var(0)
    f.write_var(0)
    f.write_var(len(boundary))
    NVEL=sum(map(lambda x:len(x[0]), boundary))
    f.write_var(NVEL)
    for b,t in boundary:
      f.write_var(len(b),t)
      for node in b: f.write_var(node)    
    f.close()

  def write_grid(self, nodes, elements, boundary):
    x=nodes.x.value_in(self.unit_length)
    y=nodes.y.value_in(self.unit_length)
    depth=nodes.depth.value_in(self.unit_length)
    element_nodes=elements.nodes+1
    boundary_nodes=[(b.nodes+1,t) for b,t in boundary] 
    self.write(x,y,depth,element_nodes,boundary_nodes)

if __name__=="__main__":
  from simple_triangulations import square_domain_sets
  
  nodes,elements,boundary=square_domain_sets(N=10)
  nodes.depth=4. | units.km

  adcirc_grid_writer().write_grid(nodes,elements,[(boundary,0)])
  
  adcirc_parameter_writer("test").write()
