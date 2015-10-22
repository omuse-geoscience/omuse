import numpy

from omuse.community.swan.interface import SwanInterface

from matplotlib import pyplot

from amuse.test.amusetest import TestWithMPI

def read_bot_data(filename="f31hari.bot",n=88,m=117):
  f=open(filename,"r")
  lines=f.readlines()
  f.close()
  
  dat=[]
  for line in lines:
    for s in line.split():
      dat.append(float(s))
  
  return numpy.transpose(numpy.array(dat).reshape((m,n)))

class TestHaringvliet(TestWithMPI):
    def test1(self):
      bathymetry=read_bot_data()


if __name__=="__main__":
  bathymetry=read_bot_data()
  bathymetry=numpy.array(bathymetry,dtype="float32")

  pyplot.ion()
  f=pyplot.figure(figsize=(8,6))
  pyplot.show()
  pyplot.imshow(numpy.transpose(bathymetry),origin="lower")
  pyplot.draw()
  
  s=SwanInterface(redirection="none")
  
  print s.get_coordinates()
  print s.get_projection_method()
  print s.get_grid_type()
  print s.get_input_grid_type()
  print s.get_calc_mode()
  print s.get_number_dimensions()
  
  print s.initialize_code()
  
  s.set_wlev(0.3)
  wlev,err=s.get_wlev()
  
  s.set_grid_xpc(6960.2)
  s.set_grid_ypc(0.)
  s.set_grid_alpc(0.)
  s.set_grid_xlenc(14789.8)
  s.set_grid_ylenc(22000.)
  mxc=98
  myc=88

  s.set_grid_mxc(mxc)
  s.set_grid_myc(myc)


  msc=32
  mdc=36

  s.set_msc(msc)
  s.set_slow(2*numpy.pi*0.0521)
  s.set_shig(2*numpy.pi*1.)
  s.set_mdc(mdc)

  print s.get_msc()
  print s.get_slow()
  print s.get_shig()
  print s.get_mdc()

  s.set_input_xp(0.)
  s.set_input_yp(0.)
  s.set_input_alp(0.)
  s.set_input_dx(250.)
  s.set_input_dy(250.)
  s.set_input_mx(87)
  s.set_input_my(116)

  print s.initialize_grids()
  
  print s.set_u10(12.)
  print s.set_wdip(8.8)

  s.set_west_boundary_spec_file("f31har01.bnd")

  s.set_use_gen3(True)
  s.set_use_breaking(True)
  s.set_use_triads(True)
  s.set_use_friction(True)
  s.set_use_uniform_wind(True)

  print s.commit_parameters()
 
  exc,err=s.get_exc_value(1)
  
  bathymetry[bathymetry==-99.]=exc
  print (bathymetry==-1.e20).sum(),exc
  input_shape=bathymetry.shape
  ii,jj=numpy.mgrid[1:input_shape[0]+1,1:input_shape[1]+1]
  print s.set_depth_regular(ii.flatten(),jj.flatten(),bathymetry.flatten())

  print s.commit_grids()

  print s.initialize_boundary()

  i,j=numpy.mgrid[1:mxc+1,1:myc+1]
  x,y,err=s.get_grid_position_regular(i.flatten(),j.flatten())

  print s.evolve_model(0.)

  i,j,k,l=numpy.mgrid[1:mxc+1,1:myc+1,1:mdc+1,1:msc+1]
  ac2,err=s.get_ac2_regular(i.flatten(),j.flatten(),k.flatten(),l.flatten())
  print err.max(),err.min()

  ac2=ac2.reshape(mxc,myc,mdc,msc)


  numpy.savez("data.npz", ac2=ac2,x=x,y=y)
  print "done"
  raw_input()
