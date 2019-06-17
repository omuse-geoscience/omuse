from amuse.units.units import *
from amuse.units.core import named_unit
from amuse.units.quantities import zero

Rearth = named('Earth radius', 'Rearth', 6371.0088 * km)
Sv=named_unit("Sverdrup","Sv",1.e6*m**3/s)
dyn=named_unit("Dyne","dyn", 1.e-5*N)
bar=named_unit("Bar", "bar", 1.e5*Pa)
dbar=deci(bar)
mbar=milli(bar)

salt=named("absolute reference salinity","Sr", g/kg)
psu=named("practical salinity unit","psu", (35.16504/35.) * salt)

nautical_mile=named_unit("Nautical Mile", "nMile", 1852.*m)
knot=named_unit("Knot","knot",nautical_mile/hour)

Celsius=named_unit("Celsius","celsius", K)

shu=named("specific humidity unit","shu", kg/kg)
ahu=named("absolute humidity unit","ahu", kg/m**3)
rhu=named("relative humidity unit","rhu", percent)

mfu=named("mass fraction unit","mfu", kg/kg)
ccu=named("cloud coverage unit","ccu", m**2/m**2)

ppb=named("parts per billion", "ppb", 1.e-9 * none)

vsmc=named("volumetric soil moisture content", "vsmc", m**3/m**3)
