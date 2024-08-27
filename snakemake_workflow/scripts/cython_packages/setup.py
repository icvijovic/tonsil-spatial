from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy
sourcefiles = ["_hierarchy_single_uint8.pyx"]
include_dirs = numpy.get_include()
extensions = Extension("MST_SINGLE_UINT8",sourcefiles,include_dirs = [include_dirs])
setup(
    name='MST_SINGLE_UINT8',
    ext_modules=cythonize(extensions),
    zip_safe=False
)
