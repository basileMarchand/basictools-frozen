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

Conda
-----

If you use conda, you can install BasicTools from the conda-forge channels [#anacondaurl]_:

    # Best practice, use an environment rather than install in the base env
    ``conda create -n my-env``
    ``conda activate my-env``
    # The actual install command
    ``conda install -c conda-forge numpy``

PIP
---

If you use pip, you can install NumPy with:

    ``> set BASICTOOLS_USE_EIGENCYEIGEN=True``
    ``> pip install  https://gitlab.com/drti/basic-tools/-/archive/1.9.1/basic-tools-1.9.1.tar.bz2``

Also when using pip, it’s good practice to use a virtual environment

Manual installation (from sources) for developers
=================================================

In the case you want to make changes to BasicTools (and potentially contribute), an installation from sources is mandatory.
The sources can be downloaded from Gitlab.com [#gitlaburlpublic]_.

    ``> git clone https://gitlab.com/drti/basic-tools.git``

Then inside the repository folder, the user must compile the c++ extensions to take profit of optimized algorithms.

    ``> python setup.py build_clib``

    ``> python setup.py build_ext --inplace``

Then the user is responsible for adding the ``BASICTOOLS_REPOSITORY/src/`` folder to the ``PYTHONPATH`` environment variables (more information on [#pythonpathdoc]_).
Or using pip for development:

    ``> pip install -e .``

The user can also install permanently using (this is not recommended):

    ``> pip install .``

The documentation for BasicTools can be compiled using sphinx

    ``> python setup.py build_sphinx``

***************
Asking for Help
***************

All questions can be addressed using the issues system of gitlab https://gitlab.com/drti/basic-tools/-/issues.

***************************
To contribute to BasicTools
***************************

To contribute to BasicTools you can fill a bug report on [#gitlaburlpublic]_ with a *minimal non working example*.
This makes the debugging easier on our side.

If you want to contribute with your code you must:

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
* eigency
* mkl
* mkl-include
* psutil
* networkx

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

C++ OPEN-SOURCE DEPENDENCIES:

* Eigen (http://eigen.tuxfamily.org)
    (the pypi eigency package has the Eigen library already inside the package, need to set the env variable BASICTOOLS_USE_EIGENCYEIGEN=True)
    ( a conda-forge package is available for eigen)

Optionals Proprietary packages (some functionalities may not be available without these packages) only for old version 1.7:

* odbAccess
* abaqusConstants


.. rubric:: Footnotes
.. [#gitlaburlpublic]  https://gitlab.com/drti/basic-tools
.. [#anacondaurl] https://anaconda.org/
.. [#scikitwindows] https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows
.. [#eigenurl] http://eigen.tuxfamily.org
.. [#pythonpathdoc] https://docs.python.org/3/using/cmdline.html\\#envvar-PYTHONPATH
