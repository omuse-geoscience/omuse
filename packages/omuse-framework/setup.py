import sys
import os

from setuptools import setup, find_packages
import support
support.use("system")
support.set_package_name("omuse")
from support.setup_codes import setup_commands
from support.version import version
from support.classifiers import classifiers

name = 'omuse-framework'
author = 'The Amuse/ Omuse Team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'https://github.com/omuse-geoscience/omuse'
install_requires = [
    'wheel>=0.32',
    'docutils>=0.6',
    'numpy>=1.2.2',
    'nose>=0.11.1',
    'mpi4py>=1.1.0',
    'h5py>=1.1.0',
    'amuse-framework>=13.1.0',
    'netCDF4>=1.4.0',
    'f90nml>=1.0.0'
]
description = 'The Oceanographic Multi-purpose Software Environment: base framework '
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []

packages = find_packages('src', exclude=["omuse.community.*"])

package_data = {
}

mapping_from_command_name_to_command_class=setup_commands()

setup(
    name=name,
    version=version,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    #~ extras_require = {
        #~ "MPI" : ["mpi4py>=1.1.0"]
    #~ },
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir = {'': 'src'},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
)
