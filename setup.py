# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import os, sys
from setuptools.command.build_ext import build_ext
from setuptools.command.build_clib import build_clib
from setuptools import setup, Extension
import configparser

# Compilation options
enable_MKL = "BASICTOOLS_DISABLE_MKL" not in os.environ
annotate = False # to generate annotation (HTML files)
useEigencyEigen = "BASICTOOLS_USE_EIGENCYEIGEN" in os.environ
__config = configparser.ConfigParser()
__config.read('setup.cfg')
debug = True if __config["build_ext"]["debug"].lower()  in ["1","true"] else False
force = True if __config["build_ext"]["force"].lower()  in ["1","true"] else False

#Cpp sources (relative to the cpp_src folder)
cpp_src = ("LinAlg/BasicOperations.cpp",
           "LinAlg/EigenTools.cpp",
           "ImplicitGeometry/ImplicitGeometryBase.cpp",
           "Containers/ElementFilter.cpp",
           "Containers/UnstructuredMesh.cpp",
           "Containers/UnstructuredMeshTools.cpp",
           "Containers/Tags.cpp",
           "FE/NativeIntegration.cpp",
           "FE/NativeNumericalWeakForm.cpp",
           "FE/DofNumbering.cpp",
           "FE/Space.cpp",
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

cpp_generators = ["cpp_generators/IntegrationRuleGenerator.py",
                  "cpp_generators/ElementNameGenerator.py",
                  "cpp_generators/SpaceGenerator.py",
                  ]
def GetBasicToolsIncludeDirs():
    try:
        import numpy
        include_dirs =[numpy.get_include(),"cpp_src" ,"."]


        import eigency
        include_dirs.extend(eigency.get_includes(include_eigen=useEigencyEigen) )
        if not useEigencyEigen:
            if "EIGEN_INC" in os.environ:
                include_dirs.append(os.environ.get('EIGEN_INC'))
        if "CONDA_PREFIX" in os.environ:
            conda_prefix = os.environ["CONDA_PREFIX"]
            include_dirs.append(os.path.join(conda_prefix, "include"))
            include_dirs.append(os.path.join(conda_prefix, "include", "eigen3"))
            include_dirs.append(os.path.join(conda_prefix, "Library", "include"))
            include_dirs.append(os.path.join(conda_prefix, "Library", "include", "eigen3"))

        if "PREFIX" in os.environ:
            conda_prefix = os.environ["PREFIX"]
            include_dirs.append(os.path.join(conda_prefix, "include"))
            include_dirs.append(os.path.join(conda_prefix, "include", "eigen3"))
            include_dirs.append(os.path.join(conda_prefix, "Library", "include"))
            include_dirs.append(os.path.join(conda_prefix, "Library", "include", "eigen3"))
        return include_dirs
    except :
        return include_dirs
class add_path():
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        sys.path.insert(0, self.path)

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            sys.path.remove(self.path)
        except ValueError:
            pass

try:
    from Cython.Build import cythonize
    from Cython.Compiler import Options
    import eigency
    with add_path('./src/'):
        for generator in cpp_generators:
            code = compile(open(generator).read(),generator,"exec")
            res = {}
            exec(code,res)
            generated_file  = res["GetGeneratedFiles"]("")
            print("generation of files : \n"+ "\n".join("cpp_src/"+str(gf)for gf in generated_file))
            res["Generate"]("cpp_src/")
            cpp_src = generated_file + cpp_src

    Options.fast_fail = True
    Options.embed = True

    modules = []

    cythonextension = []

    basictools_cpp_src_path = os.path.join("cpp_src")
    cpp_src_with_path = [os.path.join(basictools_cpp_src_path, src) for src in cpp_src]

    define_macros = []
    libraries = ["libCppBasicTools"]
    if enable_MKL or True:
        define_macros.append(("MKL_DIRECT_CALL",""))
        define_macros.append(("EIGEN_USE_MKL_VML",""))

    if debug:
        cpp_src_with_path = sorted(cpp_src_with_path,key=os.path.getmtime,reverse=True )

    ext_libraries = [['libCppBasicTools', {
               'sources': cpp_src_with_path,
               'include_dirs':GetBasicToolsIncludeDirs(),
               'macros': define_macros,
               }
    ]]

    basictools_src_path = os.path.join("src", "BasicTools")
    cython_src_with_path  = [os.path.join(basictools_src_path, src) for src in cython_src]

    for n,m in zip(cython_src,cython_src_with_path):
        cythonextension.append(Extension("BasicTools."+n.split(".pyx")[0].replace("/","."), [m],
        libraries=libraries,
        include_dirs=GetBasicToolsIncludeDirs(),
        define_macros=define_macros,
        language="c++"  ))
    modules.extend(cythonextension)

except ImportError as e:
    print(f"Compilation disabled since {e.name} package is missing")
    modules = []
    ext_libraries = []
except Exception as e:
    print("Error during generation of cpp sources ")
    print(e)

if __name__ == '__main__':
    setup(
        ext_modules=modules,
        libraries=ext_libraries,
)
