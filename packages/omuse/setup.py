from setuptools import setup
from support.version import version, main_version
from support.classifiers import classifiers


name = 'omuse'
author = 'The AMUSE/ OMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'https://github.com/omuse-geoscience/omuse'
install_requires = [
    'omuse-framework>=%s' % main_version,
    'omuse-qgmodel>=%s' % main_version
]
description = 'The Oceanographic Multi-purpose Software Environment: a package for multi-physics and multi-scale earth science simulations. '
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

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
    python_requires=">=3.5",
)
