import numpy

from amuse.units import units

from amuse.datamodel import Grid

from omuse.community.adcirc.read_grid import adcirc_file_reader

from matplotlib import pyplot,tri

class read_triangle_mesh(object):
    def __init__(self,filebase):
        self.filebase=filebase
        self.coordinates="cartesian"
    def read_grid(self):
      param=dict()

      f=adcirc_file_reader(self.filebase+'.node','r')
  
      NP,NDIM,NATTR,NBM=f.read_int(4)

      assert NDIM==2
      
      p=f.read_value_rows(*( (NP,)+(int,)+NDIM*(float,)+NATTR*(float,)+NBM*(int,) ) )
  
      index=map(lambda x: x[0],p)
      pos=map(lambda x: x[1:3],p)
      vmark=map(lambda x: x[3],p)
      index=numpy.array(index)-1
      self.p=numpy.array(pos)[index]
      self.vmark=numpy.array(vmark)[index]
  
      f.close()
      
      f=adcirc_file_reader(self.filebase+'.ele','r')
  
      NE,NPT,NEATTR=f.read_int(3)

      assert NPT==3
            
      t=f.read_value_rows(*( (NE,)+(int,)+NPT*(int,)+NEATTR*(float,)) )
      index=map(lambda x: x[0],t)
      ele=map(lambda x: x[1:4],t)
      index=numpy.array(index)-1
      self.t=numpy.array(ele)[index]-1
  
      f.close()

      #~ NOPE=f.read_int(1)
      #~ NETA=f.read_int(1)
  
      #~ self.elev_spec_boundary_seg=f.read_boundary_segments(NOPE,NETA,0)
  #~ 
      #~ NBOU=f.read_int(1)
      #~ NVEL=f.read_int(1)
      #~ self.flow_spec_boundary_seg=f.read_boundary_segments(NBOU,NVEL)
      
      param["NE"]=NE
      param["NP"]=NP
      param["NDIM"]=NDIM
      param["NATTR"]=NATTR
      param["NEATTR"]=NEATTR
      param["NPT"]=NPT
      param["NBM"]=NBM
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
        
        nodes.vmark=self.vmark
            
        elements=Grid(self.parameters["NE"])
        elements.n1=[(x[0]) for x in self.t]
        elements.n2=[(x[1]) for x in self.t]
        elements.n3=[(x[2]) for x in self.t]
            
        return nodes,elements#,elev_boundary,flow_boundary


if __name__=="__main__":
    rt=read_triangle_mesh("f32hari")
    rt.read_grid()
    nodes,elements=rt.get_sets()
    print nodes,elements

    x=nodes.x.number
    y=nodes.y.number
    vmark=nodes.vmark
    n1=elements.n1[:]
    n2=elements.n2[:]
    n3=elements.n3[:]

    elements=numpy.column_stack((n1,n2,n3))
    print elements.min(),elements.max()
    print len(nodes)
    triangulation=tri.Triangulation(x,y,elements)
        
    pyplot.figure(figsize=(14,8))
    f1=pyplot.subplot(121)
    f1.set_aspect('equal')
    pyplot.triplot(triangulation)
    f2=pyplot.subplot(122)
    f2.set_aspect('equal')
    pyplot.tripcolor(triangulation,vmark,shading='gouraud')
    pyplot.show()
  
  
