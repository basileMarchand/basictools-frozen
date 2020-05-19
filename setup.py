# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import os, sys
from setuptools.command.build_ext import build_ext
from setuptools import setup


# Compilation options
enable_MKL = True
debug = False
force = True # to force recompilation
annotate = debug # to generate annotation (HTML files)
useOpenmp = True


# Cython modules
cython_src = (
        "Linalg/EigenSolver.pyx",
        "FE/Integrators/NativeIntegration.pyx",
        "FE/WeakForms/NativeNumericalWeakForm.pyx")

try:
    from Cython.Build import cythonize

    # Compiler-dependent configuration
    # See https://stackoverflow.com/questions/30985862
    BUILD_ARGS = {}
    for compiler in ('msvc', 'gcc', 'icc'):
        if debug:
            compile_args = ['-g','-O', '-std=c++11']
        else:
            compile_args = ['-O3', '-std=c++11']

        if useOpenmp:
            compile_args.append("-fopenmp")
            if compiler == "icc":
                compile_args.append("-inline-forceinline")

        if enable_MKL:
            compile_args.append("-DMKL_DIRECT_CALL")
            compile_args.append("-DEIGEN_USE_MKL_VML")

        BUILD_ARGS[compiler] = compile_args

    import numpy
    include_path = [numpy.get_include(), os.environ['EIGEN_INC']]
    if "CONDA_PREFIX" in os.environ:
        include_path.append(os.environ["CONDA_PREFIX"] + "/include/")
    modules = cythonize(["src/BasicTools/" + src for src in cython_src], gdb_debug=debug, annotate=annotate, include_path=include_path, force=force)
    for m in modules:
        m.include_dirs = include_path + ["src/BasicTools"]
        m.extra_compile_args = compile_args
        if enable_MKL:
            m.extra_link_args.append("-lmkl_core")
            m.extra_link_args.append("-lmkl_avx")
            m.extra_link_args.append("-lmkl_intel_lp64")
            m.extra_link_args.append("-lmkl_sequential")
            m.extra_link_args.append("-lmkl_def")

except ImportError as e:
    print(f"Compilation disabled since {e.name} package is missing")
    modules = []


class build_ext_compiler_check(build_ext):
    def build_extensions(self):
        compiler = os.path.basename(self.compiler.compiler[0])
        compiler_type = self.compiler.compiler_type
        print(f"Using compiler {compiler} of type {compiler_type}")
        args = BUILD_ARGS[compiler]
        for ext in self.extensions:
            ext.extra_compile_args = args
        build_ext.build_extensions(self)

setup(
      ext_modules=modules,
      cmdclass={'build_ext': build_ext_compiler_check})
