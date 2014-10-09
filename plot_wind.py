import numpy

from amuse.units import units
from interface import single_gyre_wind_model,jans_wind_model,dijkstra_wind_model

from matplotlib import pyplot

L=1000. | units.km

a=numpy.arange(1001)/1000.*L

tau=1.

w0=single_gyre_wind_model(0.*a,a,L,tau)
w1=jans_wind_model(0.*a,a,L,tau)
w2=dijkstra_wind_model(0.*a,a,L,tau,sigma=1)
w3=dijkstra_wind_model(0.*a,a,L,tau,sigma=0.5)

pyplot.plot(w0[0],a/L,'r')
pyplot.plot(w1[0],a/L,'g')
pyplot.plot(w2[0],a/L,'b')
pyplot.plot(w3[0],a/L,'y')
pyplot.show()

