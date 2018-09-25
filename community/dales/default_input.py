import numpy
import StringIO
from omuse.units import units
from omuse.units.quantities import to_quantity

from amuse.datamodel import new_rectilinear_grid

large_scale_forcing_input="""
# #large scale forcing
# height     ug         vg         wfls       dqtdxls    dqtqyls    dqtdtls    dthlrad       
     12.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
     37.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
     62.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
     87.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    112.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    137.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    162.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    187.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    212.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    237.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    262.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    287.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    312.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    337.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    362.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    387.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    412.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    437.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    462.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    487.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    512.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    537.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    562.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    587.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    612.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    637.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    662.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    687.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    712.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    737.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    762.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    787.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    812.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    837.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    862.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    887.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    912.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    937.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    962.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
    987.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1012.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1037.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1062.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1087.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1112.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1137.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1162.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1187.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1212.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1237.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1262.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1287.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1312.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1337.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1362.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1387.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1412.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1437.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1462.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1487.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1512.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1537.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1562.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1587.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1612.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1637.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1662.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1687.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1712.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1737.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1762.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1787.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1812.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1837.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1862.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1887.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1912.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1937.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1962.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   1987.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2012.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2037.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2062.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2087.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2112.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2137.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2162.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2187.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2212.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2237.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2262.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2287.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2312.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2337.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2362.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
   2387.500      2.500      0.000      0.000      0.000      0.000      0.000      0.000
"""

initial_profile_input="""
# #profile
# zf         thl        qt         u          v          tke             
     12.500    292.500      0.005      2.500      0.000      1.000
     37.500    292.500      0.005      2.500      0.000      1.000
     62.500    292.500      0.005      2.500      0.000      1.000
     87.500    292.500      0.005      2.500      0.000      1.000
    112.500    292.500      0.005      2.500      0.000      1.000
    137.500    292.500      0.005      2.500      0.000      1.000
    162.500    292.500      0.005      2.500      0.000      1.000
    187.500    292.500      0.005      2.500      0.000      1.000
    212.500    292.500      0.005      2.500      0.000      1.000
    237.500    292.500      0.005      2.500      0.000      1.000
    262.500    292.500      0.005      2.500      0.000      1.000
    287.500    292.500      0.005      2.500      0.000      1.000
    312.500    292.500      0.005      2.500      0.000      1.000
    337.500    292.500      0.005      2.500      0.000      1.000
    362.500    292.500      0.005      2.500      0.000      1.000
    387.500    292.500      0.005      2.500      0.000      1.000
    412.500    292.500      0.005      2.500      0.000      1.000
    437.500    292.500      0.005      2.500      0.000      1.000
    462.500    292.500      0.005      2.500      0.000      1.000
    487.500    292.500      0.005      2.500      0.000      1.000
    512.500    292.500      0.005      2.500      0.000      1.000
    537.500    292.500      0.005      2.500      0.000      1.000
    562.500    292.500      0.005      2.500      0.000      1.000
    587.500    292.500      0.005      2.500      0.000      1.000
    612.500    292.500      0.005      2.500      0.000      1.000
    637.500    292.500      0.005      2.500      0.000      1.000
    662.500    292.500      0.005      2.500      0.000      1.000
    687.500    292.500      0.005      2.500      0.000      1.000
    712.500    292.500      0.005      2.500      0.000      1.000
    737.500    292.500      0.005      2.500      0.000      1.000
    762.500    292.500      0.005      2.500      0.000      1.000
    787.500    292.500      0.005      2.500      0.000      1.000
    812.500    295.075      0.005      2.500      0.000      1.000
    837.500    295.225      0.005      2.500      0.000      1.000
    862.500    295.375      0.005      2.500      0.000      1.000
    887.500    295.525      0.005      2.500      0.000      1.000
    912.500    295.675      0.005      2.500      0.000      1.000
    937.500    295.825      0.005      2.500      0.000      1.000
    962.500    295.975      0.005      2.500      0.000      1.000
    987.500    296.125      0.005      2.500      0.000      1.000
   1012.500    296.275      0.005      2.500      0.000      1.000
   1037.500    296.425      0.005      2.500      0.000      1.000
   1062.500    296.575      0.005      2.500      0.000      1.000
   1087.500    296.725      0.005      2.500      0.000      1.000
   1112.500    296.875      0.005      2.500      0.000      1.000
   1137.500    297.025      0.005      2.500      0.000      1.000
   1162.500    297.175      0.005      2.500      0.000      1.000
   1187.500    297.325      0.005      2.500      0.000      1.000
   1212.500    297.475      0.005      2.500      0.000      1.000
   1237.500    297.625      0.005      2.500      0.000      1.000
   1262.500    297.775      0.005      2.500      0.000      1.000
   1287.500    297.925      0.005      2.500      0.000      1.000
   1312.500    298.075      0.005      2.500      0.000      1.000
   1337.500    298.225      0.005      2.500      0.000      1.000
   1362.500    298.375      0.005      2.500      0.000      1.000
   1387.500    298.525      0.005      2.500      0.000      1.000
   1412.500    298.675      0.005      2.500      0.000      1.000
   1437.500    298.825      0.005      2.500      0.000      1.000
   1462.500    298.975      0.005      2.500      0.000      1.000
   1487.500    299.125      0.005      2.500      0.000      1.000
   1512.500    299.275      0.005      2.500      0.000      1.000
   1537.500    299.425      0.005      2.500      0.000      1.000
   1562.500    299.575      0.005      2.500      0.000      1.000
   1587.500    299.725      0.005      2.500      0.000      1.000
   1612.500    299.875      0.005      2.500      0.000      1.000
   1637.500    300.025      0.005      2.500      0.000      1.000
   1662.500    300.175      0.005      2.500      0.000      1.000
   1687.500    300.325      0.005      2.500      0.000      1.000
   1712.500    300.475      0.005      2.500      0.000      1.000
   1737.500    300.625      0.005      2.500      0.000      1.000
   1762.500    300.775      0.005      2.500      0.000      1.000
   1787.500    300.925      0.005      2.500      0.000      1.000
   1812.500    301.075      0.005      2.500      0.000      1.000
   1837.500    301.225      0.005      2.500      0.000      1.000
   1862.500    301.375      0.005      2.500      0.000      1.000
   1887.500    301.525      0.005      2.500      0.000      1.000
   1912.500    301.675      0.005      2.500      0.000      1.000
   1937.500    301.825      0.005      2.500      0.000      1.000
   1962.500    301.975      0.005      2.500      0.000      1.000
   1987.500    302.125      0.005      2.500      0.000      1.000
   2012.500    302.275      0.005      2.500      0.000      1.000
   2037.500    302.425      0.005      2.500      0.000      1.000
   2062.500    302.575      0.005      2.500      0.000      1.000
   2087.500    302.725      0.005      2.500      0.000      1.000
   2112.500    302.875      0.005      2.500      0.000      1.000
   2137.500    303.025      0.005      2.500      0.000      1.000
   2162.500    303.175      0.005      2.500      0.000      1.000
   2187.500    303.325      0.005      2.500      0.000      1.000
   2212.500    303.475      0.005      2.500      0.000      1.000
   2237.500    303.625      0.005      2.500      0.000      1.000
   2262.500    303.775      0.005      2.500      0.000      1.000
   2287.500    303.925      0.005      2.500      0.000      1.000
   2312.500    304.075      0.005      2.500      0.000      1.000
   2337.500    304.225      0.005      2.500      0.000      1.000
   2362.500    304.375      0.005      2.500      0.000      1.000
   2387.500    304.525      0.005      2.500      0.000      1.000
"""

def resample_array(k,x,**kwargs):
    k_=1.*numpy.arange(len(x))
    k_=k_/k_[-1]

    knew=1.*numpy.arange(k)
    knew=knew/knew[-1]
    return numpy.interp(knew,k_,x.number, **kwargs) | x.unit

def interp(x,xp,fp,**kwargs):
    x=to_quantity(x)
    xp=to_quantity(xp)
    fp=to_quantity(fp)
    return numpy.interp(x.value_in(x.unit), xp.value_in(x.unit), fp.number, **kwargs) | fp.unit

def read_initial_profile(filename=None):
    data=numpy.loadtxt(filename or StringIO.StringIO(initial_profile_input))
    zf=data[:,0] | units.m
    grid=new_rectilinear_grid( (len(zf),), cell_centers= [zf], axes_names=["z"])
    grid.QT=data[:,2] | units.shu
    grid.THL=data[:,1] | units.K
    grid.U=data[:,3] | units.m/units.s
    grid.V=data[:,4] | units.m/units.s
    grid.E12=data[:,5] | units.m/units.s
    return grid

def read_large_scale_forcings(filename=None):
    data=numpy.loadtxt(filename or StringIO.StringIO(large_scale_forcing_input))
    zf=data[:,0] | units.m
    grid=new_rectilinear_grid( (len(zf),), cell_centers= [zf], axes_names=["z"])
    grid.ug=data[:,1] | units.m/units.s
    grid.vg=data[:,2] | units.m/units.s
    grid.wfls=data[:,3]
    grid.dqtdxls=data[:,4]
    grid.dqtdyls=data[:,5]
    grid.dqtdtls=data[:,6]
    grid.dthlrad=data[:,7]
    return grid

def resample_grid(grid, kmax):
    zf=grid.positions
    assert len(zf.shape)==1
    zf=resample_array(kmax, grid.zf)
    grid=new_rectilinear_grid( (len(zf),), cell_centers= [zf], axes_names=["z"])
    grid.QT=interp(zf, zf, grid.QT)
    grid.THL=interp(zf, zf, grid.THL)
    grid.U=interp(zf, zf, grid.U)
    grid.V=interp(zf, zf, grid.V)
    grid.E12=interp(zf, zf, grid.E12)
    return grid

def write_initial_profile_file(grid, filename="prof.inp.001", kmax=None):
    if kmax and kmax!=len(grid.z):
        grid=resample_grid(grid,kmax)
    header=""" profile
 zf         thl        qt         u          v          tke             """
    numpy.savetxt( filename,numpy.column_stack((
      grid.z.value_in(units.m),
      grid.THL.value_in(units.K),
      grid.QT.value_in(units.shu),
      grid.U.value_in(units.m/units.s),
      grid.V.value_in(units.m/units.s),
      grid.E12.value_in(units.m/units.s))), header=header)

def write_large_scale_forcing_file(grid, filename="lscale.inp.001", kmax=None):
    if kmax and kmax!=len(grid.z):
        grid=resample_grid(grid,kmax)
    header=""" large scale forcing
 height     ug         vg         wfls       dqtdxls    dqtqyls    dqtdtls    dthlrad"""
    numpy.savetxt(filename,numpy.column_stack((
      grid.z.value_in(units.m),
      grid.ug.value_in(units.m/units.s),
      grid.vg.value_in(units.m/units.s),
      grid.wfls, # missing units!!
      grid.dqtdxls,
      grid.dqtdyls,
      grid.dqtdtls,
      grid.dthlrad,
      )), header=header)


if __name__=="__main__":
    pass
