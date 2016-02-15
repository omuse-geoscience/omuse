# OMUSE #

OMUSE stands for Oceanographic MUltipurpose Software Environment. It is a package to conduct numerical
experiments in Oceanography and Climate Science.

### What is this repository for? ###

This repository contains the source branch for OMUSE. 

### How do I get set up? ###

To install first install [AMUSE](www.amusecode.org) from source. Then in the AMUSE trunk, clone the OMUSE
source into /src, such that this directory contains src/amuse and src/omuse subdirectories. The community
codes of OMUSE are not yet automatically build by the framework, so you need to build them manually by going into each of:
* src/omuse/community/adcirc
* src/omuse/community/swan
* src/omuse/community/pop
* src/omuse/community/qgmodel
and typing: `make`
This has been tested on OSX and linux machines, with ifort and gfortran compilers.

In addition to the AMUSE dependencies, OMUSE needs/ can use the following packages:
* matplotlib basemap
* netCDF and netCDF for fortran and the python bindings

### Contribution guidelines ###

At the moment, this repository is private, since it contains the source of ADCIRC, which is not 
freely distributable.