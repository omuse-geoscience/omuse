# Singularity recipe for OMUSE and Jupyter notebook server
#
# This recipe uses CentOS 8.
# AMUSE uses Python3 since version 13 (Nov 2019).
#
# To build the image:
# sudo singularity build omuse.img Singularity
#
# For debugging. Generate a driectory tree named "sandbox" instead of an image file,
# and don't remove it even if the build fails.
# sudo singularity build --no-cleanup -s sandbox Singularity
#
# To start a Jupyter notebook:
#
# mkdir run
# singularity run -B src/omuse/community/dales/example/:/opt/notebooks/,run:/run/user/ ./omuse.img
#
# for troubleshooting, try preventing singularity from accessing the user's home directory
# on the host system with the option -H
#
# Tested 11.2.2020

Bootstrap: docker
From: centos:8

%setup
mkdir ${SINGULARITY_ROOTFS}/opt/splib

%files

%post
yum -y update
yum install -y epel-release dnf-plugins-core
# EPEL packages assume that the 'PowerTools' repository is enabled.
dnf config-manager --set-enabled PowerTools

yum groupinstall -y "Development Tools"
yum install -y git mercurial gcc-gfortran cmake python3-devel python3-pip python3-wheel python3-tkinter wget netcdf-devel netcdf-fortran-devel fftw-devel gmp-devel mpfr-devel gsl-devel atlas atlas-devel   perl-Digest-MD5 perl-Time-Piece perl-IO-Compress  blas-devel lapack-devel environment-modules

# yum install python3-mpi4py-openmpi - installs but not sure it works, installing with pip below instead.

# yum install openmpi-devel
# currently broken - install from source below.
# bug: https://github.com/open-mpi/ompi/issues/6671
# yum installed UCX 1.4.0, OpenMPI 4.0.1  (tested 10.2.2020)
#
## load MPI module
#export MODULEPATH=/etc/modulefiles
#eval `/usr/bin/modulecmd sh load mpi/openmpi-x86_64`

# install OpenMPI from source, as long as the packaged version is broken
# note: prefix=/usr does not work due to https://github.com/open-mpi/ompi/issues/5613
cd /opt
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.gz
tar -xzf openmpi-4.0.2.tar.gz
cd openmpi-4.0.2
./configure --prefix=/usr/local
make -j 6
make install

# make python3 the default version, makes "python" point to python3
alternatives --set python /usr/bin/python3

python3 -m pip install --upgrade --ignore-installed pip setuptools
python3 -m pip install moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil cython
python3 -m pip install jupyter ipykernel
mkdir /opt/notebooks

export DOWNLOAD_CODES=latest

# work-around for DALES for finding netCDF, for Fedora and CentOS
# where nf-config is broken and does not report the module path correctly
export FCFLAGS="-I/usr/include -I/usr/lib64/gfortran/modules -lnetcdff"

# more agressive DALES optimization options - targets the architecture where the image is built
# export SYST=gnu-fast  

cd /opt
hg clone https://bitbucket.org/omuse/omuse/
cd omuse/

python3 -m pip install -e .

# compile all codes - halts at oifs asking for git user name
# python setup.py develop_build

# build only specific codes, the ones that currently work in this Singularity recipe :
python3 setup.py build_code --code-name dales   --inplace
python3 setup.py build_code --code-name qgcm    --inplace
python3 setup.py build_code --code-name qgmodel --inplace

# codes which are currently broken, Python3 conversion reasons ?
# python3 setup.py build_code --code-name swan    --inplace
# python3 setup.py build_code --code-name cdo     --inplace

# also broken ?
# python3 setup.py build_code --code-name pop     --inplace

# requires user name / password for access to OpenIFS source code
# python3 setup.py build_code --code-name oifs     --inplace


%environment
#PYTHONPATH=/opt/amuse-framework/src/

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
#  LD_LIBRARY_PATH=/.singularity.d/libs:/usr/local/lib/
MODULEPATH=/etc/modulefiles
export LD_LIBRARY_PATH MODULEPATH

%runscript
     echo "Starting notebook..."
     # echo "Open browser to localhost:8888"
     exec /usr/local/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip=127.0.0.1 --port=8888 --no-browser

# --ip is not always needed, but was necessary on Arch Linux (with extra docker-generated network interfaces ?)


# OMUSE codes working:
# * dales
# * qgcm (hogg2006,ocean_only,hogg2006_20km)
# * qgmodel

# codes not working
# * adcirc    "make sure the (unpatched) ADCIRC source is in ../adcirc/adcirc_v51_52_07"
# * oifs      requires password, also grib_api missing in this recipe
# * pop       doesn't find netCDF
# * cdo       python 3
# * swan      "NameError: name 'getset' is not defined"
