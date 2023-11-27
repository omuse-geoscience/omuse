Data interface
==============

Using the OMUSE data interface you can expose any dataset in a 
code-like manner. The ``Data`` class has a model interface 
with a ``evolve_model`` method which handles the loading and caching of 
local data in a transparent fashion.

Variables
---------

Variables are available on standard OMUSE ``grid`` object.

Prerequisites
-------------


Example
-------

example::

    from omuse.units import units
    from omuse.community.data.interface import Data

    e=Data(path="path/to/local/data/file.nc", variables=["2m_temperature", "total_precipitation"], 
            nwse_boundingbox=[70, -15, 40, 15]| units.deg)

    print("starting date:", e.start_datetime)
    print(e.grid) # note grid has prepended the names with _ (because ERA5 variable names are not always valid python var names)

    e.evolve_model(128 | units.hour)
    val=e.grid._2m_temperature.value_in(units.K)

.. autoclass:: omuse.community.data.interface.Data
   :members:

.. _NetCDF: hhttps://docs.unidata.ucar.edu/nug/current/
