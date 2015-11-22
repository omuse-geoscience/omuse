import numpy

from amuse.units import units

from amuse.datamodel import Grid

class adcirc_file_reader(file):
  def read_string(self,n=80):
    return self.readline()[:n]
  def read_int(self,n=1):
    result=[int(l) for l in self.readline().split()[:n]]
    return result[0] if len(result)==1 else result    
  def read_value(self, *types):
    if len(types)==0: types=(float,)
    result=[x(l) for x,l in zip(types,self.readline().split())]    
    return result[0] if len(result)==1 else result    
  def read_value_rows(self,n,*types):
    result=[]
    for i in range(n):
      result.append(self.read_value(*types))
    return result
  def read_single_type_attributes(self,n,m,dtype=int):
    result=numpy.zeros((n,m),dtype)
    _n=0
    while _n<n:
      line_=self.readline().split()
      index=int(line_[0])-1
      result[index,:]=[dtype(line_[i]) for i in range(1,m+1)]
      _n+=1    
    return result
  def read_boundary_segment(self,n):
    seg=[]
    i=0
    while i<n:
      line=self.readline()
      line_=line.split()
      seg.append( int(line_[0]))
      i+=1
    return numpy.array(seg)    
  def read_boundary_segments(self,nseg,nnodes,default_type=None):
    i=0
    _nnodes=0
    result=[]
    while i<nseg:
      line=self.readline()
      line_=line.split()
      n=int(line_[0])
      if default_type is not None:
        _type=default_type
      else:
        _type=int(line_[1])
      _nnodes+=n
      seg=self.read_boundary_segment(n)
      i+=1
      result.append((_type,seg))
    assert nnodes==_nnodes
    return result

class adcirc_parameter_reader(object):  
  def __init__(self,filename="fort.15"):
    self.filename=filename
  
  def read_parameters(self,NETA=None):
    with adcirc_file_reader(self.filename,'r') as f:
      param=dict()
      param["RUNDES"]=f.read_string(32)    
      param["RUNID"]=f.read_string(24)    
      param["NFOVER"]=f.read_int()
      param["NABOUT"]=f.read_int()
      param["NSCREEN"]=f.read_int()
      param["IHOT"]=f.read_int()
      param["ICS"]=f.read_int()
      param["IM"]=f.read_int()
      if param["IM"] in [20,30,21,31]:
        param["IDEN"]=f.read_int()
      else:
        param["IDEN"]=None
      param["NOLIBF"]=f.read_int()
      param["NOLIFA"]=f.read_int()
      param["NOLICA"]=f.read_int()
      param["NOLICAT"]=f.read_int()
      param["NWP"]=f.read_int()
      param["AttrName"]=[]
      for i in range(param["NWP"]):
        param["AttrName"].append(f.read_string(80))
      param["NCOR"]=f.read_int()
      param["NTIP"]=f.read_int()
      param["NWS"]=f.read_int()
      param["NRAMP"]=f.read_int()
      param["G"]=f.read_value()    
      param["TAU0"]=f.read_value()
      if param["TAU0"]==-5.0:
        _x,_y=f.read_value(float,float)
        param["Tau0FullDomainMin"]=_x
        param["Tau0FullDomainMax"]=_y
      else:
        param["Tau0FullDomainMin"]=None
        param["Tau0FullDomainMax"]=None
      param["DTDP"]=f.read_value()
      param["STATIM"]=f.read_value()
      param["REFTIM"]=f.read_value()
      if param["NWS"] in [0,1,11]:
        pass
      elif param["NWS"] in [4,5,6,7,10,12,15]:
        param["WTIMINC"]=f.read_value()
        
  #~ WTIMINC  Supplemental Meteorological/Wave/Ice Parameters Line
      param["RNDAY"]=f.read_value()
      if param["NRAMP"] in [0,1]:
        param["DRAMP"]=f.read_value()
      elif param["NRAMP"] in [2,3,4,5,6,7,8]:
        raise Exception("tbd")
      else:
        param["DRAMP"]=None
      param["A00"],param["B00"],param["C00"]=f.read_value(float,float,float)
      if param["NOLIFA"] in [0,1]:
        param["H0"]=f.read_value()
      elif param["NOLIFA"] in [2,3]:
        param["H0"], dummy1, dummy2, param["VELMIN"]=f.read_value(float,int, int, float)
      param["SLAM0"],param["SFEA0"]=f.read_value(float,float)

      param["TAU"],param["CF"],param["HBREAK"],param["FTHETA"],param["FGAMMA"]=0.,0.,0.,0.,0.

      if param["NOLIBF"]==0:
        param["TAU"]=f.read_value()
      elif param["NOLIBF"]==1:
        param["CF"]=f.read_value()
      elif param["NOLIBF"]==2:
        param["CF"],param["HBREAK"],param["FTHETA"],param["FGAMMA"]=f.read_value(float,float,float,float)
      if param["IM"] in [0,1,2]:
        param["ESLM"]=f.read_value()
      elif param["IM"]==10:
        param["ESLM"],param["ESLC"]=f.read_value(float,float)
      param["CORI"]=f.read_value()
      param["NTIF"]=f.read_int()
      param["NBFR"]=f.read_int()
      param["BOUNTAG"]=[]
      param["AMIG"]=[]
      param["FF"]=[]
      param["FACE"]=[]
      param["BOUNDARY_FORCING_DATA"]=dict()
      for i in range(param["NBFR"]):
        param["BOUNTAG"].append(f.read_string(10))
        amig,ff,face=f.read_value(float,float,float)
        param["AMIG"].append(amig)
        param["FF"].append(ff)
        param["FACE"].append(face)
      for i in range(param["NBFR"]):
          if NETA is not None:
            tag=f.read_string(10)
            EMO=[]
            EFA=[]
            for i in range(NETA):
              emo,efa=f.read_value(float,float)
              EMO.append(emo)
              EFA.append(efa)
            param["BOUNDARY_FORCING_DATA"][tag]=(EMO,EFA)
          else:
            f.close()
            raise Exception("expect NETA to be provided")
      param["ANGINN"]=f.read_value()
# dummy read
      f.read_value(int,float,float,int)
      NSTAE=f.read_int()
      for i in range(NSTAE): f.readline()
      f.read_value(int,float,float,int)
      NSTAV=f.read_int()
      for i in range(NSTAV): f.readline()
      if param['IM']==10:
        f.read_value(int,float,float,int)
        NSTAC=f.read_int()
        for i in range(NSTAC): f.readline()
      if param['NWS']!=0:
        f.read_value(int,float,float,int)
        NSTAM=f.read_int()
        for i in range(NSTAM): f.readline()
      f.read_value(int,float,float,int)
      f.read_value(int,float,float,int)
      if param['IM']==10: f.read_int(4)
      if param['NWS']!=0: f.read_int(4)
      NHARF=f.read_int()
      for i in range(NHARF): 
        f.readline()
        f.readline()
      f.readline()
      f.readline()
      f.readline()
      param["ITITER"],param["ISLDIA"],param["CONVCR"],param["ITMAX"]=f.read_value(int,int,float,int)
      if param['IM'] in [1,11,21,31,2]:
        # continue reading 3D info
        param["IDEN"]=f.read_int()
        param['ISLIP'],param['KP']=f.read_value(int,float)
        param['Z0S'],param['Z0B']=f.read_value(float,float)
        param['ALP1'],param['ALP2'],param['ALP3']=f.read_value(float,float,float)
        param['IGC'],param['NFEN']=f.read_int(2)
        param["SIGMA"]=None
        if  param['IGC']==0:
          param['SIGMA']=f.read_value_rows(param['NFEN'],float)
        param['IEVC'],param['EVMIN'],param['EVCON']=f.read_value(int,float,float)
        if param['IEVC'] in [50,51]:
          param['THETA1'],param['THETA2']=f.read_value(float,float)
        param['EVTOT']=None
        if param['IEVC']==0: 
          param[EVTOT]=f.read_value_rows(param['NFEN'],float)
        for i in range(9):
          f.readline()
        if param['IM'] in [21,31]:
          param['RES_BC_FLAG'],param['BCFLAG_LNM'],param['BCFLAG_TEMP']=f.read_int(3)
          raise Exception("tbd: 3D baroclinic input")

#~ H0  include this line if NOLIFA =0, 1
#~ H0, INTEGER, INTEGER, VELMIN  include this line if NOLIFA =2, 3
#~ SLAM0, SFEA0
#~ TAU include this line only if NOLIBF = 0
#~ CF  include this line only if NOLIBF =1
#~ CF, HBREAK, FTHETA, FGAMMA  include this line only if NOLIBF =2
#~ ESLM  include this line only if IM =0, 1, 2
#~ ESLM, ESLC  include this line only if IM =10
#~ CORI
#~ NTIF

      

    self.parameters=param
    
class adcirc_grid_reader(object):
  
  def __init__(self,filename="fort.14",coordinates="cartesian"):
    self.filename=filename
    self.coordinates=coordinates
    if coordinates not in ["cartesian","spherical"]:
      raise Exception("coordinates must be cartesian or spherical")  
      
  def read_grid(self):
    f=adcirc_file_reader(self.filename,'r')
    param=dict()
    param["AGRID"]=f.read_string()

    NE,NP=f.read_int(2)
    
    self.p=f.read_single_type_attributes(NP,3,float)
    self.t=f.read_single_type_attributes(NE,4,int)

    assert numpy.all(self.t[:,0]==3)
            
    NOPE=f.read_int(1)
    NETA=f.read_int(1)

    self.elev_spec_boundary_seg=f.read_boundary_segments(NOPE,NETA,0)

    NBOU=f.read_int(1)
    NVEL=f.read_int(1)
    self.flow_spec_boundary_seg=f.read_boundary_segments(NBOU,NVEL)
    
    param["NE"]=NE
    param["NP"]=NP
    param["NOPE"]=NOPE
    param["NBOU"]=NBOU
    param["NETA"]=NETA
    param["NVEL"]=NVEL
    self.parameters=param
        
    f.close()

  def get_sets(self):
    nodes=Grid(self.parameters["NP"])
    if self.coordinates=="cartesian":
      nodes.x=self.p[:,0] | units.m
      nodes.y=self.p[:,1] | units.m
    else:
      nodes.lon=self.p[:,0] | units.deg
      nodes.lat=self.p[:,1] | units.deg
    nodes.depth=self.p[:,2] | units.m
        
    elements=Grid(self.parameters["NE"])
    elements.nodes=[(x[1]-1,x[2]-1,x[3]-1) for x in self.t]
    
    elev_boundary=[]
    for type_,seg in self.elev_spec_boundary_seg:
      indices=seg-1
      b=Grid(len(indices))
      b.nodes=indices
      b.type=type_
      elev_boundary.append(b)
    
    flow_boundary=[]
    for type_,seg in self.flow_spec_boundary_seg:
      indices=seg-1
      b=Grid(len(indices))
      b.nodes=indices
      b.type=type_
      flow_boundary.append(b)
    
    return nodes,elements,elev_boundary,flow_boundary
  
def assign_neighbours(nodes,elements):  
    for n in nodes:
      n.neighbours=set()
    
    for e in elements:
      p1=e.nodes[0]
      p2=e.nodes[1]
      p3=e.nodes[2]
      nodes[p1].neighbours.add(p2)
      nodes[p1].neighbours.add(p3)
      nodes[p2].neighbours.add(p1)
      nodes[p2].neighbours.add(p3)
      nodes[p3].neighbours.add(p1)
      nodes[p3].neighbours.add(p2)
     
def get_edges(elements):
    edges=set()
    for e in elements:
      for i,j in [(0,1),(1,2),(2,0)]:
        edges.add( frozenset([e.nodes[i],e.nodes[j]]) )
    return edges

if __name__=="__main__":
    a=adcirc_grid_reader()
    a.read_grid()
    nodes,elements,boundary=a.get_sets()
    assign_neighbours(nodes,elements)
    edges=get_edges(elements)
    
    print nodes[1].neighbours
    print 
    print elements
    print
    print boundary
    
    a=adcirc_parameter_reader()
    a.read_parameters()
    print a.parameters
    
