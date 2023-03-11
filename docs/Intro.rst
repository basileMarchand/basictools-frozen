******************
What is BasicTools
******************

BasicTools was primary designed as a basic set of tools to work on meshes in the context of finite element computation.
The main functionalities of the library are:

* IO routine : A set of classes to read and write meshes (and solutions fields) from/to a large variety of file formats. BasicTools does not have a proper file format mainly because existent formats provide most, if not all, of the functionalities needed.
* Mesh manipulation: Routines to filter, define, extract and manipulate meshes in many ways.
* Fields manipulation: finite element fields can be defined using different kinds of interpolation (P0/P1/P2), in the full mesh or only in restricted zones, and also at integration points. This classes have overloaded operators to make computation of quantities of interest an easy task.
* Integration: Routines for the integration of weak formulations (tangent matrices, right hand terms, integral over only a part of a mesh).
* Field transfer: Basic routine to transfer field from one mesh to another.
* Finite element solver: using all the previous tools, some basic finite element solver are  available to solve generic partial differential equations on non structured meshes.

*********************
Installing BasicTools
*********************

To install BasicTools, we recommend using a scientific Python distribution [#anacondaurl]_.

If you already have Python, you can install BasicTools with:

    ``> conda install -c drti basictools``

If you don't have Python yet, you might want to consider using Anaconda.
It's the easiest way to get started.

We are currently working to include BasicTools package on conda-forge to increase the supported platform.

Manual installation (from sources) for developers
=================================================

In the case you want to make changes to BasicTools (and potentially contribute), an installation from sources is mandatory.
The sources can be downloaded from Gitlab.com [#gitlaburlpublic]_.
Then inside the repository folder the user must compile the c++ extensions to take profit of optimized algorithms.

    ``> python setup.py build_clib``

    ``> python setup.py build_ext --inplace``

Then the user is responsible to add the ``BASICTOOLS_REPOSITORY/src/`` folder to the ``PYTHONPATH`` environment variables (more information on [#pythonpathdoc]_).
Or using pip for development:

    ``> pip install -e .``

also you can do a permanent installation using (this is not recommended):

    ``> pip install .``


For Windows
^^^^^^^^^^^

Scikit-sparse package not available in anaconda for windows, should be able to compile it for windows following [#scikitwindows]_.
I good stating point for the installation is https://github.com/EmJay276/scikit-sparse .

Extra pre-requirement:
    - Microsoft Build Tools for C++.
    - suitesparse package


***************
Asking for Help
***************

All question can be addressed using the issues system of gitlab https://gitlab.com/drti/basic-tools/-/issues.

***************************
To contribute to BasicTools
***************************

To contribute to BasicTools you can fill a bug report on [#gitlaburlpublic]_ with a *minimal non working example*.
This makes the debugging easier in ours side.

If you want to contribute with code you must.

*  clone the master branch of BasicTools from [#gitlaburlpublic]_
*  make a development branch
*  modify/created changes, commit changes
*  compile BasicTools
*  test you branch (see section :ref:`fordevs` )
*  accept the Contribution Agreement (see section :ref:`License` )
*  push your branch to the server
*  create a merge request (on the web)

************
Requirements
************

Python minimal version: 3.8

Python packages:

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

Eigency is the only dependency not already on conda-forge.
A package is available in ours channel `drti` (use `conda install -c drti eigency`).
We are currently working to have only conda-forge dependencies

C++ packages:

* Eigen [#eigenurl]_

Optionals Python packages (some functionalities may not be available without these packages):

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

Optionals Proprietary packages (some functionalities may not be available without these packages) only for old version 1.7:

* odbAccess
* abaqusConstants


.. rubric:: Footnotes
.. [#gitlaburlpublic]  https://gitlab.com/drti/basic-tools
.. [#anacondaurl] https://anaconda.org/
.. [#scikitwindows] https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows
.. [#eigenurl] http://eigen.tuxfamily.org
.. [#pythonpathdoc] https://docs.python.org/3/using/cmdline.html\\#envvar-PYTHONPATH
