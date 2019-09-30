import numpy

from matplotlib import pyplot

data=numpy.load("data.npz")

print data.files

ac2=data['ac2']

mxc,myc,mdc,msc=ac2.shape
fhigh=1.
flow=0.0512


thetas=numpy.arange(mdc+1)*2*numpy.pi/mdc
fac=(fhigh/flow)**(1./(msc-1))
fs=flow*fac**numpy.arange(msc)

print fs
print thetas

f,theta=numpy.meshgrid(fs,thetas)

print theta.shape,f.shape,ac2[0,0,:,:].shape

x=f*numpy.cos(theta)
y=f*numpy.sin(theta)

pyplot.ion()
pyplot.figure(figsize=(10,10))
pyplot.show()
iy=50
ix=20
for iy in range(myc):
  print ix,iy
  pyplot.clf()
  pyplot.pcolormesh(x,y,ac2[ix,iy,:,:])
  pyplot.draw()

raw_input()

#~ pyplot.imshow(ac2[0,0,:,:], origin="lower")
#~ 
#~ pyplot.show()
