import time

import matplotlib
matplotlib.use("TkAgg") # select appropiate backend
from matplotlib import pyplot

from omuse.units import units

from omuse.community.era5.interface import ERA5

e=ERA5(variables=["2m_temperature", "total_precipitation"], 
       nwse_boundingbox=[70, -15, 40, 15]| units.deg)

print("starting date:", e.start_datetime)
print(e.grid) # note grid has prepended the names with _ (because long names are not valid python var names)

f=pyplot.figure(figsize=(8,8))
pyplot.ion()
pyplot.show()

dt=1. | units.hour
tend=7. | units.day

while e.model_time < tend:
    print("tnow:",e.model_time)
    e.evolve_model(e.model_time+dt)
    val=e.grid._2m_temperature.value_in(units.K)
    val=0.25*(val[:-1,:-1]+val[1:,:-1]+val[:-1,1:]+val[1:,1:])
    x=e.grid.lon.value_in(units.deg)
    y=e.grid.lat.value_in(units.deg)
        
    pyplot.clf()
    pyplot.pcolormesh(x,y,val,vmin=260,vmax=288)
    pyplot.xlabel("lon (deg)")
    pyplot.ylabel("lat (deg)")
    f.canvas.flush_events()

    time.sleep(1/20.)
