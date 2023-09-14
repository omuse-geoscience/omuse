#!/usr/bin/env python3
import sys

name = sys.argv[1]
setupstring = '''#!/usr/bin/env python3
from support.classifiers import classifiers
from setuptools import setup
import support
support.use("system")
support.set_package_name("omuse")
from support.setup_codes import setup_commands
name = 'omuse-{name_lowercase}'
author = 'The AMUSE/OMUSE team'
author_email = 'info@amusecode.org'
license_ = "Apache License 2.0"
url = 'https://github.com/omuse-geoscience/omuse'
install_requires = [
    'omuse-framework',
]
description = 'The Oceanographic Multi-purpose Software Environment - {name}'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []

packages = [
    'omuse.community.{name_lowercase}',
]
package_data = {{
}}

mapping_from_command_name_to_command_class = setup_commands()

setup(
    name=name,
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
    cmdclass=mapping_from_command_name_to_command_class,
    ext_modules=extensions,
    package_dir={{
        'omuse.community.{name_lowercase}': 'src/omuse/community/{name_lowercase}',
    }},
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
)'''
print(setupstring.format(name=name, name_lowercase=name.lower()))
