from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy, os

enable_MKL = True
debug = False
force = True # to force recompilation
annotate = debug # to generate annotation (html files )

if debug:
    compile_args = ['-g','-O', '-std=c++11']
else:
    compile_args = ['-O2', '-std=c++11']

if enable_MKL:
    #compile_args.append("-DEIGEN_USE_MKL_ALL")
    compile_args.append("-DMKL_DIRECT_CALL")

include_path = [numpy.get_include(), os.environ['EIGEN_INC']]
modules = cythonize("src/BasicTools/FE/*.pyx", gdb_debug=debug, annotate=annotate, include_path=include_path, force=force)

for m in modules:
    m.include_dirs = include_path + ["src/BasicTools"]
    m.extra_compile_args = compile_args
    if enable_MKL:
        m.extra_link_args.append("-lmkl_core")
        m.extra_link_args.append("-lmkl_avx")
        m.extra_link_args.append("-lmkl_intel_lp64")
        m.extra_link_args.append("-lmkl_sequential")
        m.extra_link_args.append("-lmkl_def")

setup(ext_modules=modules, package_dir={'': 'src'})
