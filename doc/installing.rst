Installing OMUSE
================

Install required programs and libraries:

* make
* cmake
* python 3
* gcc
* gfortran
* mpi
* netCDF including Fortran bindings
* git
* mercurial
  
Set up and activate a virtual Python environment::
  
    python3 -m venv omuse_env
    source omuse_env/bin/activate

Get the OMUSE source code, install its Python dependencies and set up a development build of OMUSE, where the codes will be built in place::
  
    hg clone https://bitbucket.org/omuse/omuse/
    cd omuse/
    pip install -e .
    export DOWNLOAD_CODES=all

Build codes, select the ones needed::
  
    python setup.py build_code --code-name dales   --inplace
    python setup.py build_code --code-name cdo     --inplace
    python setup.py build_code --code-name qgcm    --inplace
    python setup.py build_code --code-name qgmodel --inplace
    python setup.py build_code --code-name swan    --inplace
    python setup.py build_code --code-name pop     --inplace
   



Alternatively, see :ref:`Singularity-section` for instructions for setting up and using a Singularity container with
OMUSE and Jupyter.
