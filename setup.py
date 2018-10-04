#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 12:08:41 2018

@author: Charlie
"""

#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, find_packages, Extension
#from Cython.Build import cythonize

#extensions = [
#    Extension("cfuncs", ["cfuncs.pyx"], 
#              include_dirs = [np.get_include()])]
#
#setup(
#    ext_modules = cythonize("cfuncs.pyx", include_path = [np.get_include()])
#)
ext_modules = [
    Extension('cfuncs', 
              ['cfuncs.pyx'])]
import numpy as np
setup(name='hogLab',
      packages=find_packages(),
      include_dirs = [np.get_include()],
      ext_modules = ext_modules,
     )