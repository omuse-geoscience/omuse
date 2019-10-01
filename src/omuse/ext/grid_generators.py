from omuse.units import units
from amuse.units.trigo import *
from amuse.datamodel import Grid
import numpy

zero=units.zero
pi=units.pi

colat=lambda fi: (0.5*pi-fi)
mlog=lambda x: numpy.log(x+(x==0.)) # m.. means masked

class displaced_pole(object):
    """
    analytic displaced pole grid of Roberts et al., Ocean Modelling 12 (2006) 16-31

    """
    def __init__(self,fi_JN=10. | units.deg, fi_pole=75.| units.deg, l_pole=320. | units.deg,
                 grid_center_function="parabolic"):
        self.fi_JN=fi_JN
        self.fi_pole=fi_pole
        self.l_pole=l_pole

        a=tan(colat(fi_pole))/tan(colat(fi_JN))

        if grid_center_function=="linear":
          df=lambda r: -a
          f= lambda r:a*(1-r)
          F= lambda r: -a*numpy.log(r)
        elif grid_center_function=="parabolic":
          df=lambda r: -2*a*(1-r)
          f= lambda r:a*(1-r)**2
          F= lambda r: 2*a*r-2*a*numpy.log(r)
        elif grid_center_function=="cubic":
          df=lambda r: -3*a*(1-r)**2
          f= lambda r: a*(1-r)**3
          F= lambda r: 6*a*r-3*a*r**2/2.-3*a*numpy.log(r)
        
        self.a=a
        self.df=df
        self.f=f
        self.F=F  

    def _forward_transform(self,fi,l):
        l=l%(2*pi)
        #polar stereographic projection
        r1=tan(colat(fi))/tan(colat(self.fi_JN))
        theta1=(l>=self.l_pole)*(l-self.l_pole)+(l<self.l_pole)*(2*pi-(self.l_pole-l))
        
        #integration constant c
        c=(theta1%pi != zero)*(-mlog(abs(tan(theta1/2)))+self.F(1.))

        #transformed longitude
        theta2=(theta1==zero)*zero +  \
               (zero<theta1)*(theta1<pi)*(2*arctan(numpy.exp(self.F(r1)-c))) + \
               (theta1==pi)*pi + \
               (pi<theta1)*(theta1<2*pi)*(2*pi-2*arctan(numpy.exp(self.F(r1)-c)))

        #orthogonal grid projection
        x=self.f(r1)+r1*cos(theta2)
        y=r1*sin(theta2)
        
        #polar coordinates
        r3=(x**2+y**2)**0.5
        theta3=arctan2(y,x)
        
        #inverse stereographic projection            
        l=theta3+self.l_pole
        l=(l<2*pi)*l + (l>=2*pi)*(l-2*pi)
        fi=colat(arctan(r3*tan(colat(self.fi_JN))))
      
        return fi,l

    def forward_transform(self,fi,l):
        mask=(fi>self.fi_JN)
        fi2=fi*mask+self.fi_JN*(1-mask)
        fi2,l2=self._forward_transform(fi2,l)
        return mask*fi2+(1-mask)*fi,mask*l2+(1-mask)*l

    def backward_iterative_solver(self,x,y,r0,t0):
        i=0
        while True:
            i+=1
            f=self.f(r0)
            r1=(y**2+(x-f)**2)**0.5
            if numpy.all(abs(r1-r0)< 1.e-15):
                break
            r0=r1
        t1=arctan2(y,x-self.f(r1))
        t1=(t1>=zero)*(t1)+(t1<zero)*(2*pi+t1)
        return r1,t1

    def _backward_transform(self,fi,l):
        l=l%(2*pi)
        
        #stereographic projection
        r3=tan(colat(fi))/tan(colat(self.fi_JN))    
        theta3=(l>=self.l_pole)*(l-self.l_pole)+(l<self.l_pole)*(2*pi-(self.l_pole-l))      

        #inverse grid projection
        x=r3*cos(theta3)
        y=r3*sin(theta3)
        r1,theta2=self.backward_iterative_solver(x,y,r3,theta3)

        #integration constant c
        c=(theta2%pi != zero)*(-mlog(abs(tan(theta2/2)))+self.F(r1))

        #transform back theta
        theta1=(theta2==zero)*zero +  \
               (zero<theta2)*(theta2<pi)*(2*arctan(numpy.exp(self.F(1.)-c))) + \
               (theta2==pi)*pi + \
               (pi<theta2)*(theta2<2*pi)*(2*pi-2*arctan(numpy.exp(self.F(1.)-c)))

        #inverse stereographic projection
        l=theta1+self.l_pole
        l=(l<2*pi)*l + (l>=2*pi)*(l-2*pi)
        fi=colat(arctan(r1*tan(colat(self.fi_JN))))

        return fi,l

    def backward_transform(self,fi,l):
        mask=(fi>self.fi_JN)
        fi2=fi*mask+self.fi_JN*(1-mask)
        fi2,l2=self._backward_transform(fi2,l)
        return mask*fi2+(1-mask)*fi,mask*l2+(1-mask)*l

    def generate_grid(self,n=10,m=10, FI0=-80. | units.deg, FI1=80. | units.deg, 
                           L0=0. | units.deg, L1=360. | units.deg,
                           dFI=None, dL=None):
        if dFI is None:
          dFI=(FI1-FI0)/n
        if dL is None:
          dL=(L1-L0)/m
        FI,L=numpy.mgrid[in_deg(FI0):in_deg(FI1):in_deg(dFI),
                           in_deg(L0):in_deg(L1):in_deg(dL)]
        FI=FI|units.deg
        L=L|units.deg
        n,m=FI.shape
        fi,l=self.forward_transform(FI,L)
        grid=Grid((n,m))
        grid.FI=FI
        grid.L=L
        grid.lat=to_deg(fi)
        grid.lon=to_deg(l)
        return grid

def test_transform():
    gg=displaced_pole()
    N=0.5
    extent = (0,360.,0.05,89.95)
    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    xs=xs.flatten() | units.deg
    ys=ys.flatten() | units.deg
    x0=xs
    y0=ys
    ys,xs=gg.forward_transform(ys,xs)
    xs1=xs
    ys1=ys
    ys,xs=gg.backward_transform(ys,xs)
    a1=abs(xs-x0)
    a2=abs(xs-2*pi-x0)
    dth=(a1<a2)*a1+(a1>=a2)*a2
    assert abs(ys-y0).max().number < 1.e-12
    assert in_rad(dth).max() < 1.e-12
    print("test ok")


if __name__=="__main__":
    from matplotlib import pyplot
    gg=displaced_pole()

    grid=gg.generate_grid(FI1=90|units.deg,dFI=5|units.deg,dL=5| units.deg)
#~ 
    lat=grid.lat.value_in(units.deg)
    lon=grid.lon.value_in(units.deg)

    pyplot.plot(lon,lat,'r.')
    pyplot.show()    

    #~ test_transform()
