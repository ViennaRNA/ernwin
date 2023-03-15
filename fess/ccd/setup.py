from __future__ import absolute_import
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(name="cytvec", 
                sources = ["rot_mat_c.c", "cytvec.pyx"])]

setup(
  name = 'Cython Vector Functions',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
