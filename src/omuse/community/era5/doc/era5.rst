ERA5 interface
==============

Using the OMUSE interface to ERA5 you can access the CDS ERA5 dataset in a 
code-like manner. The ``ERA5`` class has a model interface 
with a ``evolve_model`` method which handles the download and caching of 
ERA5 data in a transparent fashion.

Variables
---------

Currently the ERA5 interface exposes single level variables from ERA5_.
Variables are available on standard OMUSE ``grid`` object.

Prerequisites
-------------

``ERA5`` needs the CDSAPI_ installed as well as a valid 
access key in `` ~/.cdsapirc``

Example
-------

example::

    from omuse.units import units
    from omuse.community.era5.interface import ERA5

    e=ERA5(variables=["2m_temperature", "total_precipitation"], 
            nwse_boundingbox=[70, -15, 40, 15]| units.deg)

    print("starting date:", e.start_datetime)
    print(e.grid) # note grid has prepended the names with _ (because ERA5 variable names are not always valid python var names)

    e.evolve_model(128 | units.hour)
    val=e.grid._2m_temperature.value_in(units.K)

.. autoclass:: omuse.community.era5.interface.ERA5
   :members:

.. _ERA5: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
.. _CDSAPI: https://pypi.org/project/cdsapi/
