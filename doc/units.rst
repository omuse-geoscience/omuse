.. _units:

Units in OMUSE
==============

OMUSE code interfaces use quantities with units.
This has the advantages that the units of any
quantitiy is explicitly specified, and that
automatic conversion can be done between different units.

Unit example::

  from omuse.units import units

  # units are attached to a number with the | operator:
  velocity = 3 | units.m / units.s
  mass = 0.002 | units.kg

  print('velocity', velocity)

  # a quantity can be separated into a number and a unit
  print('number:', velocity.number, 'unit:', velocity.unit)
  print()

  print('mass', mass)

  # get value in a different unit
  print(mass.value_in(units.g), 'g')
  print()

  # arithmetic on quantities:
  print('momentum:', mass * velocity)

output::
  
  velocity 3 m / s
  number: 3 unit: m / s

  mass 0.002 kg
  2.0 g

  momentum: 0.006 m * kg * s**-1



