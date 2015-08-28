#!/usr/bin/env python
from distutils.core import setup, Extension
from setup import *

CdiLib_module = Extension('_CdiLib',
                  sources=['cdilib_wrap.c'],
                  extra_compile_args = INCFLAGS,
                  library_dirs = LDFLAGS,
                  extra_objects = ['../../src/cdilib.o'],
                  extra_link_args = LIBS,
                  )

setup (name = 'CdiLib',
       version = '0.1',
       author      = "Ralf Mueller",
       description = """pyhton bindings to CDI function library""",
       ext_modules = [CdiLib_module],
       py_modules = ["CdiLib"],
       )
