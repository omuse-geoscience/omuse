# Singularity recipe for OMUSE and Jupyter notebook server
# 29.7.2019

Bootstrap: docker
From: centos:7

%setup
mkdir ${SINGULARITY_ROOTFS}/opt/splib

%files

%post
yum -y update
yum install -y epel-release
yum groupinstall -y "Development Tools"
yum install -y git mercurial gcc-gfortran cmake python-devel python-pip python-wheel wget openmpi-devel mpi4py-openmpi netcdf-devel netcdf-fortran-devel fftw-devel gmp-devel mpfr-devel gsl-devel atls atlas-devel blas-devel lapack-devel perl-Digest-MD5 perl-Time-Piece perl-IO-Compress

pip install --upgrade --ignore-installed pip setuptools
pip install moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil cython

pip install jupyter
mkdir /opt/notebooks

# load MPI module
export MODULEPATH=/etc/modulefiles
eval `/usr/bin/modulecmd sh load mpi/openmpi-x86_64`


export DOWNLOAD_CODES=latest

# work-around for DALES for finding netCDF, for Fedora and CentOS 
export FCFLAGS="-I/usr/include -I/usr/lib64/gfortran/modules"

# more agressive DALES optimization options - targets the architecture where the image is built
# export SYST=gnu-fast  


cd /opt
hg clone https://bitbucket.org/omuse/omuse/
cd omuse/

pip install -e .

# compile all codes - halts at oifs asking for git user name
# python setup.py develop_build

# build only specific codes, the ones that currently work in this Singularity recipe :
python setup.py build_code --code-name dales   --inplace
python setup.py build_code --code-name cdo     --inplace
python setup.py build_code --code-name qgcm    --inplace
python setup.py build_code --code-name qgmodel --inplace
python setup.py build_code --code-name swan    --inplace


%environment
#PYTHONPATH=/opt/amuse-framework/src/

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
#  LD_LIBRARY_PATH=/.singularity.d/libs:/usr/local/lib/
MODULEPATH=/etc/modulefiles
export LD_LIBRARY_PATH MODULEPATH

%runscript
     echo "Starting notebook..."
     # echo "Open browser to localhost:8888"
     exec /usr/bin/jupyter notebook --notebook-dir=/opt/notebooks --port=8888 --no-browser


# OMUSE codes working:
# * cdo
# * dales
# * qgcm (hogg2006,ocean_only,hogg2006_20km)
# * qgmodel
# * swan

# codes not working
# * adcirc    "make sure the (unpatched) ADCIRC source is in ../adcirc/adcirc_v51_52_07"
# * oifs      requires password, also grib_api missing in this recipe
# * pop       doesn't find netCDF

