# Singularity recipe for OMUSE and Jupyter notebook server
# see https://omuse.readthedocs.io/en/latest/singularity.html
# 
# This recipe uses CentOS 8 and supports Python 3.
# AMUSE uses Python3 since version 13 (Nov 2019).
#
# To build the image:
# sudo singularity build omuse.img Singularity
#
# For debugging: generate a directory tree named "sandbox" instead of an image file,
# and don't remove it even if the build fails.
# sudo singularity build --no-cleanup -s sandbox Singularity
# The build log from compiling DALES or other OMUSE codes can then be seen at
# sandbox/opt/omuse/build.log . If the image building fails, the sandbox directory
# may be named rootfs-xxxxx.
#
#
# To start a Jupyter notebook:
#  mkdir run
#  singularity run --contain -B examples:/opt/notebooks,run:/run/user omuse.img
#
# To start a shell in the container:
#  mkdir run
#  singularity shell --contain omuse.img -B run:/run/user
#
# --contain prevents the outside environment and home directory from being imported,
# which can lead to problems e.g. if OMUSE is also installed on the host system
#
# -B run:/run/user is needed for OpenMPI to work.
#
# Tested 18.8.2020

Bootstrap: docker
From: centos:8

%setup
mkdir ${SINGULARITY_ROOTFS}/opt/splib
mkdir ${SINGULARITY_ROOTFS}/opt/omuse

%files
# Files to copy into the image. 
# The following is sufficient to install OMUSE in the container.
# By default, we instead do an independent git clone of omuse in the container
# since build results present in the source tree (on the host) can interfere 
# with building in the container
#
# LICENSE /opt/omuse/LICENSE
# MANIFEST.in /opt/omuse/MANIFEST.in
# NOTICE /opt/omuse/NOTICE
# omuse_logo.png /opt/omuse/omuse_logo.png
# pyproject.toml /opt/omuse/pyproject.toml
# README.md /opt/omuse/README.md
# setup.py /opt/omuse/setup.py
# src /opt/omuse/src
# support /opt/omuse/support
# test /opt/omuse/test
# examples /opt/omuse/examples
# doc /opt/omuse/doc

%post
yum -y update
yum install -y epel-release dnf-plugins-core
# EPEL packages assume that the 'PowerTools' repository is enabled.
dnf config-manager --set-enabled PowerTools

yum groupinstall -y "Development Tools"
yum install -y git mercurial gcc-gfortran cmake python3-devel python3-pip python3-wheel python3-tkinter wget netcdf-devel netcdf-fortran-devel fftw-devel gmp-devel mpfr-devel gsl-devel atlas atlas-devel   perl-Digest-MD5 perl-Time-Piece perl-IO-Compress  blas-devel lapack-devel environment-modules

yum install -y openmpi-devel

## load MPI module
export MODULEPATH=/etc/modulefiles
eval `/usr/bin/modulecmd sh load mpi/openmpi-x86_64`

yum install -y python3-mpi4py-openmpi 
# mpi4py can be installed with pip, automatically as an omuse dependency,
# but currently that doesn't work (MPI spawn fails) 

# The shipped OpenMPI was broken on early versions of centos8 - then installed from source below.
# bug: https://github.com/open-mpi/ompi/issues/6671
# yum installed UCX 1.4.0, OpenMPI 4.0.1                                (tested 10.2.2020, was broken)
#               ucx-1.6.1-1.el8.x86_64, openmpi-devel-4.0.2-2.el8.x86_64.rpm (tested 17.8.2020, works) 
#
# install OpenMPI from source, as long as the packaged version is broken
# note: prefix=/usr does not work due to https://github.com/open-mpi/ompi/issues/5613
# cd /opt
# wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.gz
# tar -xzf openmpi-4.0.2.tar.gz
# cd openmpi-4.0.2
# ./configure --prefix=/usr/local
# make -j 6
# make install

# make python3 the default version, makes "python" point to python3
alternatives --set python /usr/bin/python3

python3 -m pip install --upgrade --ignore-installed pip setuptools
python3 -m pip install moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil cython
python3 -m pip install jupyter ipykernel
mkdir /opt/notebooks

export DOWNLOAD_CODES=1

# work-around for DALES for finding netCDF, for Fedora and CentOS
# where nf-config is broken and does not report the module path correctly
# bug: https://bugzilla.redhat.com/show_bug.cgi?id=1738541
# The bug is supposed to be fixed in centos8, but compilation still fails without this line. 
# Probably -I/usr/lib64/gfortran/modules is missing, obtainable with nf-config --fflags. 
export FCFLAGS="-I/usr/include -I/usr/lib64/gfortran/modules -lnetcdff"

# optionally enable more agressive DALES optimization options
# - targets the architecture where the image is built, not portable.
# export SYST=gnu-fast  

# clone OMUSE from github. See the %files section for an alternative
cd /opt
git clone https://github.com/omuse-geoscience/omuse
cd /opt/omuse

# install omuse dependencies with pip
python3 -m pip install -e .

# compile all codes - halts at oifs asking for git user name
# python setup.py develop_build

# build only specific codes, the ones that currently work in this Singularity recipe :
python3 setup.py build_code --code-name dales   --inplace
python3 setup.py build_code --code-name qgcm    --inplace
python3 setup.py build_code --code-name qgmodel --inplace
python3 setup.py build_code --code-name swan    --inplace
python3 setup.py build_code --code-name cdo     --inplace

# Codes currently not working in the image :

# OpenIFS: requires user name / password for access to OpenIFS source code
# python3 setup.py build_code --code-name oifs     --inplace

# pop: doesn't find netCDF
# python3 setup.py build_code --code-name pop     --inplace



%environment
#PYTHONPATH=/opt/amuse-framework/src/

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
#  LD_LIBRARY_PATH=/.singularity.d/libs:/usr/local/lib/
MODULEPATH=/etc/modulefiles
export LD_LIBRARY_PATH MODULEPATH


%runscript
     echo "Starting omuse singularity container..."
     echo "Open browser to localhost:8888"
     exec /usr/local/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip=127.0.0.1 --port=8888 --no-browser

# --ip is not always needed, but was necessary on Arch Linux (with
    extra docker-generated network interfaces ?)


# OMUSE codes working:
# * dales
# * qgcm (hogg2006,ocean_only,hogg2006_20km)
# * qgmodel
# * swan
# * cdo

# codes not working
# * adcirc    "make sure the (unpatched) ADCIRC source is in ../adcirc/adcirc_v51_52_07"
# * oifs      requires password, also grib_api missing in this recipe
# * pop       doesn't find netCDF
