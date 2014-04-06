import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_module = Extension(name="gait2de",
                       sources=["gait2de.pyx",
                                "gait2de_al.c"],
                       include_dirs=[numpy.get_include()])

setup(name="gait2de",
      cmdclass = {'build_ext': build_ext},
      ext_modules = [ext_module])
