from distutils.core import setup
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import Extension
import os
import numpy as np
from Cython.Build import cythonize
import cython_gsl

__version__ = '1.0'

gsl_funcs_ext = Extension('CBASS_Simulator_Modules.Tools.gsl_funcs',
                          ['CBASS_Simulator_Modules/Tools/gsl_funcs.pyx'],
                          libraries=cython_gsl.get_libraries(),
                          library_dirs=[cython_gsl.get_library_dir()],
                          include_dirs=[cython_gsl.get_include()],
                          extra_compile_args=['-fopenmp'],
                          extra_link_args=['-fopenmp']
                      )

pysla = Extension(name = 'CBASS_Simulator_Modules.Tools.pysla', 
                        sources = ['CBASS_Simulator_Modules/Tools/pysla.f90',
                                   'CBASS_Simulator_Modules/Tools/sla.f'])


config = {'name':'CBASS_Simulator_Modules',
          'version':__version__,
          'packages':['CBASS_Simulator_Modules.Simulators',
                      'CBASS_Simulator_Modules.Tools'],
          'ext_modules':cythonize([pysla,gsl_funcs_ext])}

setup(**config)
