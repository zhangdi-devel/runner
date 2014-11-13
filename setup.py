# $File: setup.py $
# $LastChangedDate:  $
# $Rev:  $
# copied from gaow's setup file

import sys, imp, os
for item in ['faulthandler']:
    try:
        imp.find_module(item)
    except ImportError:
        sys.exit('Cannot build package: missing module "{}"!'.format(item))

from src import VERSION
from distutils.core import setup
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py


NAME = "Runner"

setup(name = NAME,
    version = VERSION,
    description = "A general pipeline to run jobs in cluster",
    author = "Di Zhang",
    packages = [NAME],
    scripts = ['src/run.py'],
    package_dir = {NAME:'src'},
    cmdclass = {'build_py': build_py }
)
