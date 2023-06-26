******************
What is BasicTools
******************

BasicTools was primary designed as a basic set of tools to work on meshes in the context of finite element computations.
The main functionalities of the library are:

* IO support: A set of classes to read and write meshes (and solutions fields) from/to a large variety of file formats. BasicTools does not have a proper file format mainly because existent formats provide most, if not all, of the functionalities needed.
* Mesh manipulation: Routines to filter, define, extract and manipulate meshes in many ways.
* Fields manipulation: Finite element fields can be defined using different kinds of interpolation (P0/P1/P2), in the full mesh or only in restricted zones, and also at integration points. This classes have overloaded operators to make computation of quantities of interest an easy task.
* Integration: Routines for the integration of weak formulations (tangent matrices, right hand terms, integral over only a part of a mesh).
* Field transfer: Basic routine to transfer field from one mesh to another.
* Finite element solver: Using all the previous tools, some basic finite element solvers are available to solve generic partial differential equations on unstructured meshes.

*********************
Installing BasicTools
*********************

Conda
-----

If you use conda, you can install BasicTools from the conda-forge channels [#anacondaurl]_.

A good practice is to use a virtual environment rather than modifying the base environment:

.. code-block::

    conda create -n my-env
    conda activate my-env

The actual install command is:

.. code-block::

    conda install -c conda-forge basictools

PIP
---

The pip installation requires a local compilation, so you need to have a C++ (C++17 compatible) compiler installed locally on your system.
To compile and install BasicTools (version 1.9.4 in this case) with pip:

.. code-block::

    set BASICTOOLS_USE_EIGENCYEIGEN=True                                # or "export" depending on your shell
    pip install eigency mkl numpy sympy mkl-include cython wheel
    pip install BasicTools@https://gitlab.com/drti/basic-tools/-/archive/1.9.4/basic-tools-1.9.4.tar.bz2

or for the latest master version:

.. code-block::

    set BASICTOOLS_USE_EIGENCYEIGEN=True                                # or "export" depending on your shell
    pip install eigency mkl numpy sympy mkl-include cython wheel
    pip install BasicTools@git+https://gitlab.com/drti/basic-tools.git

The user can set the environment variable `PREFIX` to point to external libraries (like mkl and eigen header). For advanced configuration please refer to the setup.py file on the git repository.

It is also a good practice to use a virtual environment when using pip.

.. note::
    We can not guarantee that every combination of operating system, Python version and packaging system works.

Manual installation (from sources) for developers
-------------------------------------------------

In the case you want to make changes to BasicTools or contribute new features, an installation from sources is mandatory.
The sources can be downloaded from gitlab.com [#gitlaburlpublic]_:

.. code-block::

    git clone https://gitlab.com/drti/basic-tools.git

Then inside the repository folder, the user must compile the C++ extensions to have the optimized algorithms available:

.. code-block::

    python setup.py build_clib
    python setup.py build_ext --inplace

Then the user must add the ``BASICTOOLS_REPOSITORY/src/`` folder to the ``PYTHONPATH`` environment variables (more information on [#pythonpathdoc]_).
Or using pip for development:

.. code-block::

    pip install -e .

The user can also install permanently using (Not recommended):

.. code-block::

    pip install .

The documentation for BasicTools can be compiled using sphinx:

.. code-block::

    python setup.py build_sphinx

***************
Asking for Help
***************

Questions can be submitted using the Issues system of Gitlab [#gitlaburlpublicissues]_.

Bugs should ideally be reported with a *minimal non working example* to make debugging easier for the developers.

**************************
Contributing to BasicTools
**************************

If you want to contribute some code you must:

*  clone the master branch of BasicTools from [#gitlaburlpublic]_
*  create a development branch
*  modify/create changes, commit changes
*  compile BasicTools
*  test your branch (see section :ref:`fordevs`)
*  accept the Contribution Agreement (see section :ref:`License`)
*  push your branch to Gitlab
*  create a merge request

************
Requirements
************

Python Dependencies
-------------------

Python minimal version: 3.8.
Some features may be unavailable when optional packages are not installed.

+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|Module Name     |Version|Compile|Run|Debug|Doc|Optional|Notes                                      |
+================+=======+=======+===+=====+===+========+===========================================+
|python          | >=3.8 |\*     |\* |\*   |\* |        |Supported distributions are: conda         |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|numpy           | >=1.20|\*     |\* |     |   |        |                                           |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|scipy           | >=1.9 |       |\* |     |   |        |sparse (coo_matrix),                       |
|                |       |       |   |     |   |        |spatial ( KDTree, delaunay, ConvexHull)    |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|sympy           |       |\*     |\* |     |   |        |matrices, Symbols, lambdify, Derivative,   |
|                |       |       |   |     |   |        |symplify                                   |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|cython          |       |\*     |   |     |   |        |Compilation of c++ extensions              |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|vtk             |       |       |\* |\*   |   |        |stlReader, UnstructuredMeshFieldOperations,|
|                |       |       |   |     |   |        |ImplicitGeometryObjects, vtkBridge         |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|eigency         | >=2   |\*     |\* |     |   |        |Compilation and run of c++ extensions      |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|mkl             |       |\*     |\* |     |   |        |Can be deactivated at compilation using    |
|                |       |       |   |     |   |        |the env variable : BASICTOOLS_DISABLE_MKL  |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|mkl-include     |       |\*     |   |     |   |        |Can be deactivated at compilation using    |
|                |       |       |   |     |   |        |the env variable : BASICTOOLS_DISABLE_MKL  |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|psutil          |       |       |   |\*   |\* |        |memory usagen and cpu_count()              |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|scikit-sparse   |       |       |\* |     |   |\*      |Linear solver: Cholesky "cholesky"         |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|matplotlib      |       |       |   |\*   |   |\*      |plot shape function for debugin            |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|pyamg           |       |       |\* |     |   |\*      |linear solver: Algebraic Multigrid "AMG"   |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|h5py            |       |       |\* |     |   |\*      |xdmf Reader/Writer                         |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|meshio          |       |       |\* |     |   |\*      |main usage in MeshIOBridge.py (derivated   |
|                |       |       |   |     |   |        |usage in Mesh File Converter)              |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|sphinx          |       |       |   |     |\* |\*      |Documentation Generation                   |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|sphinx-rtd-theme|       |       |   |     |\* |\*      |Documentation Generation                   |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|breathe         |       |       |   |     |\* |        |cmake documentation integration            |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|setuptools-scm  |       |\*     |   |     |\* |\*      |not sure we use it                         |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|pyvista         |       |       |\* |     |   |\*      |pyvista bridge                             |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|networkx        | >=3   |       |\* |     |   |        |only use in UnstructuredMeshGraphTools.py  |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+
|mpi4py          |       |       |   |     |   |\*      |only use in MPIInterface.py                |
+----------------+-------+-------+---+-----+---+--------+-------------------------------------------+

C++ Dependencies
----------------

+---------+-------+-------+---+-----+---+--------+----------------------------------------------------+
|Name     |Version|Compile|Run|Debug|Doc|Optional|Notes                                               |
+=========+=======+=======+===+=====+===+========+====================================================+
|eigen    | >=3.4 |\*     |   |     |   |        | For compilation of the C++ extensions              |
+---------+-------+-------+---+-----+---+--------+----------------------------------------------------+
|boost-cpp|       |\*     |   |     |   |        | For the compilation of the extension field transfer|
+---------+-------+-------+---+-----+---+--------+----------------------------------------------------+

External Dependencies
---------------------

+------+-------+-------+---+-----+---+--------+-------------------------------------------------+
|Name  |Version|Compile|Run|Debug|Doc|Optional|Notes                                            |
+======+=======+=======+===+=====+===+========+=================================================+
|cmake | >=3.8 |(\*)   |   |     |\* |        | for the cpp documentation generation            |
|      |       |       |   |     |   |        | (*) experimental cmake extensions compilation   |
+------+-------+-------+---+-----+---+--------+-------------------------------------------------+
|abaqus|       |       |   |     |   |\*      | odb reader. This feature is deprecated          |
|      |       |       |   |     |   |        | (only available on python 2.7, BasicTools 1.7.2)|
+------+-------+-------+---+-----+---+--------+-------------------------------------------------+

.. rubric:: Footnotes
.. [#gitlaburlpublic] https://gitlab.com/drti/basic-tools
.. [#gitlaburlpublicissues] https://gitlab.com/drti/basic-tools/-/issues
.. [#anacondaurl] https://anaconda.org/
.. [#scikitwindows] https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows
.. [#eigenurl] http://eigen.tuxfamily.org
.. [#pythonpathdoc] `https://docs.python.org/3/using/cmdline.html\#envvar-PYTHONPATH <https://docs.python.org/3/using/cmdline.html\#envvar-PYTHONPATH>`_
