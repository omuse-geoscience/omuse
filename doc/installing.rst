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
  
    git clone https://github.com/omuse-geoscience/omuse/
    cd omuse/
    pip install -e .
    export DOWNLOAD_CODES=1

Build codes, select the ones needed::
  
    python setup.py build_code --code-name dales   --inplace
    python setup.py build_code --code-name cdo     --inplace
    python setup.py build_code --code-name qgcm    --inplace
    python setup.py build_code --code-name qgmodel --inplace
    python setup.py build_code --code-name swan    --inplace
    python setup.py build_code --code-name pop     --inplace
   
or try to build all of them::
  
    python setup.py develop_build

Install Jupyter in the virtual environment, and make the virtual environment's Python available as a kernel in Jupyter 

    python -m pip install ipykernel matplotlib
    python -m ipykernel install --user --name=omuse-env


Alternatively, see :ref:`Singularity-section` for instructions for setting up and using a Singularity container with
OMUSE and Jupyter.


Code versions
-------------

For DALES, there are additional options controlling which version of the code is used:
setting `DOWNLOAD_CODES=1` performs a shallow checkout of a single tag, while `DOWNLOAD_CODES="all"`
clones the whole DALES git repository, which is useful for development.
The environment variable `DALES_GIT_TAG` can be used to control which
branch or version tag to check out.
By default the variable points to a version tag in the DALES repository, which is tested to
work with the current OMUSE. 
