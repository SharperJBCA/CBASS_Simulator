from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy import get_include
from Cython.Build import cythonize
import os

# SET SLALIB PATH TO LOCATION OF YOUR LOCAL LIBRARIES
# Run > python setup.py build_ext --inplace
try:
    slalib_path = os.environ['SLALIB_LIBS']
except KeyError:
    slalib_path = '/star/lib' # default path of Manchester machines


binFuncs = Extension(name='binFuncs',
                     include_dirs=[get_include()],
                     sources=['binFuncs.pyx'])

pysla = Extension(name = 'pysla', 
                  sources = ['pysla.f90'],
                  libraries=['sla'],
                  library_dirs =['{}'.format(slalib_path)],
                  f2py_options = [],
                  extra_f90_compile_args=['-L{}'.format(slalib_path)])

extensions = [binFuncs,pysla]
setup(ext_modules=cythonize(extensions))
