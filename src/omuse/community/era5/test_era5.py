import time

from omuse.units import units

from omuse.community.era5.interface import era5cached

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot

e=era5cached(variables=["2m_temperature", "total_precipitation"], nwse_boundingbox=[70, -15, 40, 15]| units.deg)

print(e.grid) # note grid has prepended the names with _ (because long names are not valid python var names)

f=pyplot.figure()
pyplot.ion()
pyplot.show()

#~ pyplot.imshow(e.grid.land_sea_mask)

dt=1. | units.hour
tend=28. | units.day


while e.model_time < tend:
    print("tnow:",e.model_time)
    e.evolve_model(e.model_time+dt)
    t2=e.grid._2m_temperature
        
    pyplot.clf()
    pyplot.imshow(t2.value_in(units.K),vmin=260,vmax=288)
    f.canvas.flush_events()

    #~ input()
    time.sleep(1/20.)
