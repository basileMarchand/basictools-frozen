from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy,os

enable_MKL = True
compile_args = ['-O3', '-std=c++11']#

if enable_MKL:
    #compile_args.append("-DEIGEN_USE_MKL_ALL")
    compile_args.append("-DMKL_DIRECT_CALL")
    compile_args.append("-I/softs/python/miniconda2/include/")

modules = cythonize("FE/*.pyx", gdb_debug=True,annotate=True, include_path = [numpy.get_include(),os.environ['EIGEN_INC'] ]  )

for m in modules:
    m.include_dirs = [ numpy.get_include(),os.environ['EIGEN_INC'],"." ]
    m.extra_compile_args=compile_args
    if enable_MKL:
        m.extra_link_args.append("-L/softs/python/miniconda2/pkgs/mkl-2017.0.3-0/lib/")
        m.extra_link_args.append("-lmkl_core")
        m.extra_link_args.append("-lmkl_avx")
        m.extra_link_args.append("-lmkl_intel_lp64")
        m.extra_link_args.append("-lmkl_sequential")
        m.extra_link_args.append("-lmkl_def")


setup( ext_modules = modules)
