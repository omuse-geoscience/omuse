# OMUSE #

OMUSE stands for Oceanographic MUltipurpose Software Environment. It is a 
package to conduct numerical experiments in Oceanography and Climate 
Science.

### What is this repository for? ###

This repository contains the source tree for OMUSE. 

### How do I get set up? ###

Easiest way is to use a pip developer install:

- (optional) first do a develop install of [AMUSE](http://www.amusecode.org), 
- Set the environment variable `DOWNLOAD_CODES` to `latest`.
- then do `pip install -e .` from the root of the package
- Then type `python setyp.py build_codes --inplace` to build the codes. 
- The file `build.log` will report errors in the build process.

The community codes of OMUSE can be build manually by going into ea:

 + src/omuse/community/adcirc
 + src/omuse/community/swan
 + etc

and typing: first `make download` (for some) and then `make`

OMUSE has been tested on OSX and linux machines, with ifort and gfortran 
compilers, on desktop machines and on the Carthesius supercomputer.

In addition to the AMUSE dependencies, OMUSE needs/ can use the following 
packages:

 + matplotlib basemap
 + netCDF and netCDF for fortran and the python bindings
 + GRIB_API

### Reporting issues ###

Issues can be reported at the OMUSE issue tracker; for framework issues, 
report them at the AMUSE [repository](https://github.com/amusecode/amuse).

### Contribution guidelines ###

Contributions are welcome. Note that most framework development happens at 
the AMUSE [repository](https://github.com/amusecode/amuse) A primer for 
writing code interfaces and other documentation can be found on the amuse 
website (www.amusecode.org).
