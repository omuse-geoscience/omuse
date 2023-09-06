#!/usr/bin/env python3
import sys

name = sys.argv[1]
setupstring = '''[build-system]

requires = [ "setuptools>=65.0.0", "wheel", "omuse-framework", "setuptools_scm", "numpy" ]

[project]
name = "omuse-{name_lowercase}"
dynamic = ["version"]

[tool.setuptools_scm]
write_to = "src/omuse/version.py"
root = "../.."
'''
print(setupstring.format(name=name, name_lowercase=name.lower()))
