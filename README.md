# OMUSE #

OMUSE stands for Oceanographic MUltipurpose Software Environment. It is a package to conduct numerical
experiments in Oceanography and Climate Science.

### What is this repository for? ###

This repository contains the source tree for OMUSE. 

### How do I get set up? ###

To install first install [AMUSE](http://www.amusecode.org) from source. Then in the AMUSE trunk, clone the OMUSE
source into /src, such that this directory contains src/amuse and src/omuse subdirectories. The community
codes of OMUSE can be build manually by going into each of:

 + src/omuse/community/adcirc
 + src/omuse/community/swan
 + src/omuse/community/pop
 + src/omuse/community/qgmodel

and typing: (first `make download` for some) `make`
This has been tested on OSX and linux machines, with ifort and gfortran compilers, on desktop machines and the Carthesius
supercomputer.

In addition to the AMUSE dependencies, OMUSE needs/ can use the following packages:

 + matplotlib basemap
 + netCDF and netCDF for fortran and the python bindings

### Reporting issues ###

Issues can be reported at the OMUSE issue tracker; for framework issues, report them at the 
AMUSE [repository](https://github.com/amusecode/amuse).

### Contribution guidelines ###

Contributions are welcome. Note that most framework development happens at the AMUSE [repository](https://github.com/amusecode/amuse)
A primer for writing code interfaces and other documentation can be found on the amuse website
(www.amusecode.org).
