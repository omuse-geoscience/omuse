from setuptools import setup
from support.classifiers import classifiers
from setuptools_scm import get_version

version_full = get_version(
    root='../..',
    relative_to=__file__,
)

version = '.'.join(version_full.split('.')[:2]

name = 'omuse'
author = 'The AMUSE/ OMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'https://github.com/omuse-geoscience/omuse'
install_requires = [
    'omuse-framework>=%s' % version,
    'omuse-qgmodel>=%s' % version
]
description = 'The Oceanographic Multi-purpose Software Environment: a package for multi-physics and multi-scale earth science simulations. '
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

try:
    from src.omuse.version import version
    use_scm_version = False
    setup_requires = []
except ImportError:
    version = False
    setup_requires = ['setuptools_scm']
    use_scm_version = {
        "root": "../..",
        "relative_to": __file__,
    }

setup(
    name=name,
    use_scm_version=use_scm_version,
    setup_requires=setup_requires,
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
