# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import os, sys
from setuptools.command.build_ext import build_ext
from setuptools import setup, Extension
from setuptools.config import read_configuration
import configparser 

# Compilation options
enable_MKL = "BASICTOOLS_DISABLE_MKL" not in os.environ
annotate = False # to generate annotation (HTML files)
useOpenmp = "BASICTOOLS_DISABLE_OPENMP" not in os.environ
__config = configparser.ConfigParser()
__config.read('setup.cfg')
debug = True if __config["build_ext"]["debug"].lower()  in ["1","true"] else False 
force = True if __config["build_ext"]["force"].lower()  in ["1","true"] else False 

#Cpp sources (relative to the cpp_src folder)
cpp_src = ("LinAlg/EigenTools.cpp",
           "ImplicitGeometry/ImplicitGeometryBase.cpp",
           "Containers/ElementFilter.cpp",
           "Containers/UnstructuredMesh.cpp",
           "Containers/UnstructuredMeshTools.cpp",
           "Containers/Tags.cpp",
           "FE/NativeIntegration.cpp",
           "FE/DofNumbering.cpp",
           "FE/Space.cpp",
           )

#Cpp sources generated files
cpp_src += ("Containers/ElementNames.cpp", # <-- this files is generated from ElementNames.py
               )

# Cython modules
cython_src = (
    "Linalg/NativeEigenSolver.pyx",
    "FE/Integrators/NativeIntegration.pyx",
    "FE/WeakForms/NativeNumericalWeakForm.pyx",
    "Containers/NativeUnstructuredMesh.pyx",
    "FE/Numberings/NativeDofNumbering.pyx",
    "FE/Spaces/NativeSpace.pyx",
    "Containers/NativeFilters.pyx",)

try:
    from Cython.Build import cythonize
    import eigency
    from Cython.Compiler import Options

    from BasicTools.Containers.ElementNames import GeneratertElementNamesCpp

    GeneratertElementNamesCpp()

    Options.fast_fail = True
    Options.embed = True

    modules = []

    cythonextension = []

    basictools_cpp_src_path = os.path.join("cpp_src")
    cpp_src_with_path = [os.path.join(basictools_cpp_src_path, src) for src in cpp_src]
    cppextensions = [Extension("libCppBasicTools", cpp_src_with_path)]
    cythonextension.extend(cppextensions)

    basictools_src_path = os.path.join("src", "BasicTools")
    cython_src_with_path  = [os.path.join(basictools_src_path, src) for src in cython_src]
    
    for n,m in zip(cython_src,cython_src_with_path):
        cythonextension.append(Extension("BasicTools."+n.split(".pyx")[0].replace("/","."), [m],libraries=["CppBasicTools"], include_dirs=["./cpp_src/"]))
        
    modules.extend(cythonize(cythonextension, gdb_debug=debug, annotate=annotate, force=force))

except ImportError as e:
    print(f"Compilation disabled since {e.name} package is missing")
    modules = []

# Compiler-dependent configuration
# See https://stackoverflow.com/questions/30985862
class build_ext_compiler_check(build_ext):
    def get_ext_filename(self, ext_name):
        #strip sufix for the libCppBasicTools.xxxxx.so
        filename = super().get_ext_filename(ext_name)
        if filename.find("libCppBasicTools") == 0 :
            return filename.split(".")[0]+"."+filename.split(".")[-1]
        return filename 

    def finalize_options(self):
        super().finalize_options()
        self.library_dirs.append(self.build_temp)
        self.library_dirs.append(self.build_temp.replace("temp.","lib."))

    def build_extensions(self):
        if os.name == 'nt':
            compiler = os.path.basename(self.compiler.compiler_type)
        else:
            compiler = os.path.basename(self.compiler.compiler[0])
        compiler_type = self.compiler.compiler_type
        print(f"Using compiler {compiler} of type {compiler_type}")
        compile_args = self._compile_args(compiler)
        link_args = self._link_args(compiler)
        include_dirs = self._include_dirs(compiler)
        for ext in self.extensions:
            ext.extra_compile_args.extend(compile_args)
            ext.extra_link_args.extend(link_args)
            ext.include_dirs.extend(include_dirs)
        build_ext.build_extensions(self)

    def _compile_args(self, compiler):
        if debug:
            if compiler =="msvc":
                compile_args = ['/Od']
            else:
                compile_args = ['-g', '-O', '-std=c++11']
        else:
            if compiler =="msvc":
                compile_args = ['/O2']
            else:
                compile_args = ['-O3', '-std=c++11']
        if useOpenmp:
            if compiler == "icc":
                compile_args.append("-fopenmp")
                compile_args.append("-inline-forceinline")
            elif compiler == "msvc":
                compile_args.append("/openmp")
            else:
                compile_args.append("-fopenmp")
        if enable_MKL:
            compile_args.append("-DMKL_DIRECT_CALL")
            compile_args.append("-DEIGEN_USE_MKL_VML")
        return compile_args

    def _link_args(self, compiler):
        link_args = []
        if compiler != "msvc":
            if enable_MKL:
                link_args.append("-lmkl_core")
                link_args.append("-lmkl_avx")
                link_args.append("-lmkl_intel_lp64")
                link_args.append("-lmkl_sequential")
                link_args.append("-lmkl_def")
            if useOpenmp:
                link_args.append("-lgomp")
        return link_args

    def _include_dirs(self, _):
        import numpy
        include_dirs =[numpy.get_include(),"cpp_src" ,"."]
        include_dirs.extend(eigency.get_includes(include_eigen=False) )
        if "EIGEN_INC" in os.environ:
            include_dirs.append(os.environ.get('EIGEN_INC'))
        if "CONDA_PREFIX" in os.environ:
            conda_prefix = os.environ["CONDA_PREFIX"]
            include_dirs.append(os.path.join(conda_prefix, "include"))
            include_dirs.append(os.path.join(conda_prefix, "include", "eigen3"))
            include_dirs.append(os.path.join(conda_prefix, "Library", "include"))
        include_dirs.append(os.path.join("src", "BasicTools"))
        return include_dirs
        
if __name__ == '__main__':
    setup(
        ext_modules=modules,
        cmdclass={'build_ext': build_ext_compiler_check})
