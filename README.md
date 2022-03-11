# OMUSE #

OMUSE stands for Oceanographic MUltipurpose Software Environment. It is a 
package to conduct numerical experiments in oceanography and other Earth 
sciences. Example OMUSE applications can be found in the examples 
[repository](https://github.com/omuse-geoscience/omuse-examples).

### For whom is OMUSE? ###

OMUSE aims to be useable by any researcher or student with a basic knowledge of 
PYTHON.

### What is this repository for? ###

This repository contains the source tree for OMUSE, including OMUSE specific framework
components and community codes.

### How do I get set up? ###

While there are some packages available on [pipy](www.pypi.org), we recommend at the moment 
to do a pip developer install:

- setup a python environment, e.g. using virtualenv, and activate it.
- in a suitable working directory clone the [AMUSE](https://github.com/amusecode/amuse) repository: `git clone https://github.com/amusecode/amuse`
- go into the created directory: `cd amuse`
- do the developer install from here: `pip install -e . [MPI]` The MPI is optional. 
- Going back to the working directory (`cd ..`) also clone the OMUSE repository: `git clone https://github.com/omuse-geoscience/omuse`,
- go into the source directory `cd omuse` and set the environment variable `DOWNLOAD_CODES`, e.g. `export DOWNLOAD_CODES=latest`.
- now, do `pip install -e .` from the root of the package
- type `python setup.py build_codes --inplace` to build the codes. 
- the file `build.log` will report any errors in the build process.

This installs amuse-devel and omuse-devel. The community codes of OMUSE can 
be build manually by going into each directory:

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

### Documentation ###

Documentation can be found [here](https://omuse.readthedocs.io). In addition the base  [AMUSE documentation](https://amuse.readthedocs.io) can be consulted.

### Reporting issues ###

Issues can be reported at the OMUSE issue tracker; for framework issues, 
report them at the AMUSE [repository](https://github.com/amusecode/amuse).

### Contribution guidelines ###

Contributions are welcome. Note that most framework development happens at 
the AMUSE [repository](https://github.com/amusecode/amuse) A primer for 
writing code interfaces and other documentation can be found on the amuse 
[website](www.amusecode.org).
