from amuse.units.units import *
from amuse.units.core import named_unit
from amuse.units.quantities import zero

Sv=named_unit("Sverdrup","Sv",1.e6*m**3/s)
dyn=named_unit("Dyne","dyn", 1.e-5*N)
bar=named_unit("Bar", "bar", 1.e5*Pa)
dbar=deci(bar)
mbar=milli(bar)

salt=named("absolute reference salinity","Sr", g/kg)
psu=named("practical salinity unit","psu", (35.16504/35.) * salt)

nautical_mile=named_unit("Nautical Mile", "nMile", 1852*m)
knot=named_unit("Knot","knot",nautical_mile/hour)

Celsius=named_unit("Celsius","celsius", K)
