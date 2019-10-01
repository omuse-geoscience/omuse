import numpy

from omuse.units import units
from amuse.datamodel import UnstructuredGrid, Particles,Particle

def square_domain(L=1.,N=10):
    dx=1./N
    
    x,y= numpy.mgrid[0:L+dx/2:dx,0:L+dx/2:dx]
    x1,y1= numpy.mgrid[dx/2:L+dx/2:dx,dx/2:L+dx/2:dx]
    
    s=x.shape
    s1=x1.shape
    
    x_=x.flatten()
    y_=y.flatten()
    x1_=x1.flatten()
    y1_=y1.flatten()
    
    index=numpy.arange(len(x_)+len(x1_))
    
    l1=numpy.product(x.shape)
    
    i=index[:l1].reshape(s)
    i1=index[l1:].reshape(s1)
    
    e1=numpy.zeros((N**2,3),dtype='i')
    e2=numpy.zeros((N**2,3),dtype='i')
    e3=numpy.zeros((N**2,3),dtype='i')
    e4=numpy.zeros((N**2,3),dtype='i')
    
    e1[:,0]=i[:-1,:-1].flatten()
    e1[:,1]=i[1:,:-1].flatten()
    e1[:,2]=i1[:,:].flatten()
    
    e2[:,0]=i[1:,:-1].flatten()
    e2[:,1]=i[1:,1:].flatten()
    e2[:,2]=i1[:,:].flatten()
    
    e3[:,0]=i[1:,1:].flatten()
    e3[:,1]=i[:-1,1:].flatten()
    e3[:,2]=i1[:,:].flatten()
    
    e4[:,0]=i[:-1,:-1].flatten()
    e4[:,1]=i1[:,:].flatten()
    e4[:,2]=i[:-1,1:].flatten()
    
    x=numpy.concatenate([x_,x1_])
    y=numpy.concatenate([y_,y1_])
    
    elements=numpy.zeros((4*N**2,3),dtype='i')
    elements[0::4,:]=e1
    elements[1::4,:]=e2
    elements[2::4,:]=e3
    elements[3::4,:]=e4
  
    boundaries=[xx.flatten() for xx in [i[:,0],i[-1,:],i[::-1,-1],i[0,::-1]] ]
  
    return x,y,elements,boundaries


#~ def unstructured_square_domain(L=1000 | units.km,N=10):
    #~ x,y,triangles,edges=square_domain(N=N)
    #~ nodes=UnstructuredGrid(len(x))
    #~ elements=UnstructuredGrid(len(triangles))
    #~ nodes.x=x*L
    #~ nodes.y=y*L
    #~ elements.nodes=triangles
    #~ boundaries=Particles()
    #~ for e in edges:
      #~ ee=nodes[e]
      #~ boundary=Particle(nodes=ee)
      #~ boundaries.add_particle(boundary)
    #~ return nodes,elements,boundaries

def unstructured_square_domain_sets(L=1000 | units.km,N=10):
    x,y,triangles,boundaries=square_domain(N=N)
    nodes=UnstructuredGrid(len(x))
    elements=UnstructuredGrid(len(triangles))
    nb=[len(x) for x in boundaries]
    bnodes=numpy.zeros(numpy.sum(nb)-3, dtype='i')
    bnodes[:nb[0]]=boundaries[0][:]
    bnodes[nb[0]:nb[0]+nb[1]-1]=boundaries[1][1:]
    bnodes[nb[0]+nb[1]-1:nb[0]+nb[1]+nb[2]-2]=boundaries[2][1:]
    bnodes[nb[0]+nb[1]+nb[2]-2:]=boundaries[3][1:]
    boundary=UnstructuredGrid(numpy.sum(nb)-3)
    nodes.x=x*L
    nodes.y=y*L
    elements.nodes=triangles
    boundary.nodes=bnodes
    boundary.type=0
    return nodes,elements,[], [boundary]

def unstructured_square_domain(L=1000 | units.km,N=10):
    x,y,triangles,edges=square_domain(N=N)
    nodes=UnstructuredGrid(len(x))
    elements=UnstructuredGrid(len(triangles))
    nodes.x=x*L
    nodes.y=y*L
    for i in range(len(triangles)):
        elements[i].nodes=nodes[triangles[i,:]]
    boundaries=Particles()
    for e in edges:
      boundary=Particle(nodes=nodes[e])
      boundaries.add_particle(boundary)
    return nodes,elements,boundaries

if __name__=="__main__":
    nodes,elements,boundaries=unstructured_square_domain(L=1. ,N=10)
    
    for i,b in enumerate(boundaries):
      b.nodes.vmark=i+1

    print(elements[0].nodes.indices())
    print(boundaries.nodes)
