# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

import os, sys
from setuptools.command.build_ext import build_ext
from setuptools.command.build_clib import build_clib
from setuptools.command.install import install
from setuptools import setup, Extension, Command
from distutils.command.build import build
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
    include_dirs = []
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
    import numpy
    include_dirs.extend([numpy.get_include(),"cpp_src" ,"."])
    import eigency
    include_dirs.extend(eigency.get_includes(include_eigen=useEigencyEigen) )
    if not useEigencyEigen:
        if "EIGEN_INC" in os.environ:
            include_dirs.append(os.environ.get('EIGEN_INC'))
    return list(set(include_dirs))

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
    with add_path('./src/'):
        for generator in cpp_generators:
            code = compile(open(generator).read(),generator,"exec")
            res = {}
            exec(code,res)
            generated_file  = res["GetGeneratedFiles"]("")
            cpp_src = generated_file + cpp_src

    modules = []
    cpp_src_with_path = [os.path.join("cpp_src", src) for src in cpp_src]


    if debug:
        cpp_src_with_path = sorted(cpp_src_with_path,key=os.path.getmtime,reverse=True )

    ext_libraries = [['libCppBasicTools', {
               'sources': cpp_src_with_path,
               }
    ]]

    cython_src_with_path  = [os.path.join("src", "BasicTools",src) for src in cython_src]

    for n,m in zip(cython_src,cython_src_with_path):
        modules.append(Extension("BasicTools."+n.split(".pyx")[0].replace("/","."), [m],
        libraries=["libCppBasicTools"],
        language="c++",
        ))

except ImportError as e:
    print(f"Compilation disabled since {e.name} package is missing")
    modules = []
    ext_libraries = []
except Exception as e:
    print("Error during generation of cpp sources ")
    print(e)

extra_compile_args = {
            'unix': ['-fopenmp'],
            'msvc': ['/openmp']
    }
extra_link_args = {
            'unix': ['-fopenmp'],
            'msvc': []
    }

class GenerateCommand(Command):
    description = "custom generate command that generate the c++ sources from python "
    user_options = []
    def initialize_options(self):
        self.cwd = None
    def finalize_options(self):
        self.cwd = os.getcwd()
    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd
        with add_path('./src/'):
            for generator in cpp_generators:
                code = compile(open(generator).read(),generator,"exec")
                res = {}
                exec(code,res)
                generated_file  = res["GetGeneratedFiles"]("")
                print("generation of files : \n"+ "\n".join("cpp_src/"+str(gf)for gf in generated_file))
                res["Generate"]("cpp_src/")
                #cpp_src = generated_file + cpp_src

class my_build(build):
    def run(self):
        self.run_command("generate")
        build.run(self)

class my_build_ext(build_ext):
    def build_extensions(self):
        for ext in self.extensions:
            ext.include_dirs.extend(GetBasicToolsIncludeDirs())
            ctype = self.compiler.compiler_type
            ext.extra_compile_args = extra_compile_args.get(ctype, [])
            ext.extra_link_args = extra_link_args.get(ctype, [])

        build_ext.build_extensions(self)

class my_build_clib(build_clib):

    def build_libraries(self,libraries):
        self.run_command("generate")
        define_macros = []
        if enable_MKL:
            define_macros.append(("MKL_DIRECT_CALL",None))
            define_macros.append(("EIGEN_USE_MKL_VML",None))

        include_dirs_BasicTools = GetBasicToolsIncludeDirs()
        for (lib_name, build_info) in libraries:
            include_dirs = build_info.get("include_dirs",[])
            include_dirs.extend(include_dirs_BasicTools)
            build_info["include_dirs"] = include_dirs

            macros = build_info.get('macros',[])
            macros.extend(define_macros)
            build_info["macros"] = macros

            ctype = self.compiler.compiler_type
            cflags = build_info.get('cflags',[])
            cflags.extend(extra_compile_args.get(ctype, []))
            build_info["cflags"] = cflags


        build_clib.build_libraries(self,libraries)



if __name__ == '__main__':
    setup(
        ext_modules=modules,
        libraries=ext_libraries,
        cmdclass={ 'build': build,'build_ext':my_build_ext,'build_clib': my_build_clib,'generate': GenerateCommand},
        data_files=[("ParaViewPlugins",["extras/BasicToolsParaViewBridge.py"])]
)
