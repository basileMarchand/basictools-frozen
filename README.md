
What is BasicTools
==================

BasicTools was primary designed as a basic set of tools to work on meshes in the context of finite element computation.
The main functionalities of the library are:

* IO support: A set of classes to read and write meshes (and solutions fields) from/to a large variety of file formats. BasicTools does not have a proper file format mainly because existent formats provide most, if not all, of the functionalities needed.
* Mesh manipulation: Routines to filter, define, extract and manipulate meshes in many ways.
* Fields manipulation: Finite element fields can be defined using different kinds of interpolation (P0/P1/P2), in the full mesh or only in restricted zones, and also at integration points. This classes have overloaded operators to make computation of quantities of interest an easy task.
* Integration: Routines for the integration of weak formulations (tangent matrices, right hand terms, integral over only a part of a mesh).
* Field transfer: Basic routine to transfer field from one mesh to another.
* Finite element solver: Using all the previous tools, some basic finite element solvers are available to solve generic partial differential equations on unstructured meshes.


Important URLs
==============

- Documentation: https://basictools.readthedocs.io/en/latest/
- Conda-forge Package: https://anaconda.org/conda-forge/basictools
- Sources: https://gitlab.com/drti/basic-tools
- Conda-forge feedstock: https://github.com/conda-forge/basictools-feedstock


Installing BasicTools
=====================

Conda
-----

If you use conda, you can install BasicTools from the conda-forge channel:

Best practice, use an environment rather than install in the base env

    conda create -n my-env
    conda activate my-env

The actual install command

    conda install -c conda-forge basictools

PIP
---

The pip installation requires a local compilation, so you need to have a C++ (C++17 compatible) compiler installed locally on your system.
To compile and install BasicTools (version 1.9.12 in this case) with pip:

    wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.zip
    unzip boost_1_82_0.zip
    set BASICTOOLS_USE_EIGENCYEIGEN=1
    set BASICTOOLS_EXTERNAL_BOOST_DIR=%cd%\boost_1_82_0
    pip install eigency mkl numpy sympy mkl-include cython wheel
    pip install BasicTools@https://gitlab.com/drti/basic-tools/-/archive/1.9.11/basic-tools-1.9.12.tar.bz2



or for the latest master version:

    wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.zip
    unzip boost_1_82_0.zip
    set BASICTOOLS_USE_EIGENCYEIGEN=1
    set BASICTOOLS_EXTERNAL_BOOST_DIR=%cd%\boost_1_82_0
    pip install eigency mkl numpy sympy mkl-include cython wheel
    pip install BasicTools@git+https://gitlab.com/drti/basic-tools.git

>Note
On linux/OsX you must:
>- Change the `set` to `export` or `setenv` depending on your os/shell
>- Change the `%cd%` to `$PWD`  depending on your os/shell

The user can set the environment variable `PREFIX` to point to external libraries (like mkl and eigen header). For advance configuration please read the setup.py file on the git repository.

It is also good practice to use a virtual environment when using pip.


>Note
We can not guarantee that all combinations of OS, Python Versions, packaging systems works.


For more complex installation (from sources) for developers please read the documentation.


Asking for help
===============

All questions can be addressed using the Issues system of Gitlab https://gitlab.com/drti/basic-tools/-/issues.


Dependencies
============

    python minimal version: 3.8

    PYTHON OPEN-SOURCE DEPENDENCIES

    * numpy >= 1.20
    * scipy
    * sympy
    * cython
    * eigency >=2
    * mkl
    * mkl-include
    * psutil

    Optionals Python packages (some functionalities may not be available without this packages):

    * scikit-sparse
    * vtk
    * matplotlib
    * pyamg
    * h5py
    * meshio
    * sphinx
    * sphinx-rtd-theme
    * setuptools-scm
    * pyvista
    * networkx >=3
    * mpi4py

    C++ OPEN-SOURCE DEPENDENCIES:

    * Eigen (http://eigen.tuxfamily.org)
      (the pypi eigency package has the Eigen library already inside the package)
      ( a conda-forge package is available for eigen)

    THIRD-PARTY PROPRIETARY DEPENDENCIES (optional):

    * odbAccess and abaqusConstants ( Abaqus )

Projects using BasicTools
=========================

[OpenPisco Home Page](https://gitlab.com/openpisco/openpisco), topology optimization using the level set method ([OpenPisco Documentation](https://openpisco.readthedocs.io)).
[GenericROM](https://gitlab.com/drti/genericrom), Reduced Order Modeling library  ([GenericROM Documentation](https://genericrom.readthedocs.io/en/latest/)).
