******************
What is BasicTools
******************

BasicTools was primary designed as a basic set of tools to work on meshes in the context of finite element computation.
The main functionalities of the library are:

* IO support: A set of classes to read and write meshes (and solutions fields) from/to a large variety of file formats. BasicTools does not have a proper file format mainly because existent formats provide most, if not all, of the functionalities needed.
* Mesh manipulation: Routines to filter, define, extract and manipulate meshes in many ways.
* Fields manipulation: Finite element fields can be defined using different kinds of interpolation (P0/P1/P2), in the full mesh or only in restricted zones, and also at integration points. This classes have overloaded operators to make computation of quantities of interest an easy task.
* Integration: Routines for the integration of weak formulations (tangent matrices, right hand terms, integral over only a part of a mesh).
* Field transfer: Basic routine to transfer field from one mesh to another.
* Finite element solver: Using all the previous tools, some basic finite element solvers are available to solve generic partial differential equations on unstructured meshes.

**************
Important URLs
**************

- Documentation: https://basictools.readthedocs.io/en/latest/
- Conda-forge Package: https://anaconda.org/conda-forge/basictools
- Sources: https://gitlab.com/drti/basic-tools
- Conda-forge feedstock: https://github.com/conda-forge/basictools-feedstock

*************************
Projects using BasicTools
*************************

OpenPisco, topology optimization using the level set method: https://gitlab.com/openpisco/openpisco, https://openpisco.readthedocs.io
GenericROM, Reduced Order Modeling library: https://gitlab.com/drti/genericrom, https://genericrom.readthedocs.io/en/latest/


************
Dependencies
************

    python minimal version: 3.8

    PYTHON OPEN-SOURCE DEPENDENCIES

    * numpy >= 1.20
    * scipy
    * sympy
    * pyparsing
    * cython
    * scikit-learn
    * scikit-sparse
    * vtk
    * eigency>=1.78
    * mkl
    * mkl-include
    * psutil
    * networkx


    Optionals Python packages (some functionalities may not be available without this packages):

    * matplotlib
    * pyamg
    * h5py
    * meshio
    * sphinx
    * sphinx-rtd-theme
    * setuptools-scm
    * pyvista
    * sksparse
    * CGNS
    * paraview
    * pywin32 [Only for windows]


    C++ OPEN-SOURCE DEPENDENCIES:

    * Eigen (http://eigen.tuxfamily.org)
      (the pypi eigency package has the Eigen library already inside the package)
      ( a conda-forge package is available for eigen)

    THIRD-PARTY PROPRIETARY DEPENDENCIES (optional):

    * odbAccess and abaqusConstants ( Abaqus )


    FOR WINDOWS:
      install Microsoft Visual C++ Build Tools to use eigen,

*********************
Installing BasicTools
*********************

Conda
-----

If you use conda, you can install BasicTools from the conda-forge channels [#anacondaurl]_:

    # Best practice, use an environment rather than install in the base env
    ``> conda create -n my-env``
    ``> conda activate my-env``
    # The actual install command
    ``> conda install -c conda-forge numpy``

PIP
---

If you use pip, you can install BasicTools with:

    ``> set BASICTOOLS_USE_EIGENCYEIGEN=True``
    ``> pip install  https://gitlab.com/drti/basic-tools/-/archive/1.9.1/basic-tools-1.9.1.tar.bz2``

Also when using pip, it is good practice to use a virtual environment.

For more complex installation (from sources) for developers please read the documentation.

***************
Asking for help
***************

All questions can be addressed using the Issues system of Gitlab https://gitlab.com/drti/basic-tools/-/issues.
