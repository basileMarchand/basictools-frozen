---
title: 'BasicTools: a numerical simulation toolbox'
tags:
  - Python
  - C++
  - mesh
  - fields
  - finite elements
  - post-treatment
authors:
  - name: Felipe Bordeu^[corresponding author]
    orcid: 0000-0002-0640-5485
    affiliation: 1
  - name: Fabien Casenave
    orcid: 0000-0002-8810-9128
    affiliation: 1
  - name: Julien Cortial
    orcid: 0000-0002-7181-9561
    affiliation: 1
affiliations:
  - name: Safran Tech, Digital Sciences & Technologies Department, Rue des Jeunes Bois, Châteaufort, 78114 Magny-Les-Hameaux, France
    index: 1
date: 05 february 2023
bibliography: paper.bib

---

# Summary

Numerical simulations of physical phenomena can be computed by many (commercial/free) software packages, but despite the apparent variety, they rely on a relatively small set of operations during the preparation, exploitation and post-process of these simulations, e.g. handling and modifying meshes and fields. BasicTools is a Python library designed to address these supporting tasks. It features an efficient data model for meshes and field objects and input/output routines compatible with various formats. A finite element engine allows to assemble abstract variational formulations and integrate fields on volumes and surfaces.

BasicTools is actively used in artificial intelligence and model order reduction [@ROM-net; @mca26010017; @UQindustrialDesign; @datatargetVAE], topology optimization [@nardoni], and material sciences [@pymicro] projects.

# Statement of need

Industrial design tasks often rely on numerical simulation workflows involving different softwares, each providing its own specific post-processing tools. Common tasks like transferring computed fields from one simulation software to another must be routinely implemented, with subtle variations. This limits interoperability and increases complexity.

Our solution to these concerns is BasicTools, which introduces a data model for meshes and related physical fields. This model can be populated using different readers and exported using various writers: no new mesh or solution format is forced upon the user. The data-oriented design of BasicTools allows high performance operations using a high-level language (Python with NumPy). BasicTools contains tools to convert meshes to other "in-memory" formats (VTK [@VTK4], PyVista [@sullivan2019pyvista], MeshIO [@meshio], CGNS [@cgns], Gmsh [@gmsh] ). This allows mixing (and reusing) treatments available in other frameworks.
Additionally, BasicTools provides many algorithms including mesh manipulation, a finite element engine and field projection operators.


# State of the field

In the computational fluid dynamics community, the CFD General Notation System (CGNS) [@cgns] format is standard. In the solid mechanics community, to the authors' knowledge, no such standard exists.
One may consider respectively VTK and MeshIO for mesh manipulation and file format conversion, but the post-process of integration point data, which is a paramount important in solid mechanics, would not be possible. Most of the existent tool do the simple, but potentially dangerous approach, of extrapolating the integration point values to the nodes of the mesh or averaging in every cell. This can lead to misinterpretation of the solution and wrong engineering decisions. Finally, only a few finite element engines allow assembling abstract variational formulations on arbitrary geometries (e.g. FreeFem++ [@freefempp] and FEniCS [@fenics]).

# Overview

The main features of BasicTools are:

- Meshes (in module `Containers`):
  `ConstantRectilinearMesh` and `UnstructuredMesh` encapsulate respectively the data model for constant rectilinear and unstructured mesh types. Unstructured meshes are efficient: elements are stored using only one array for each element type. Both mesh types can feature nodes and element tags. Many functions are available for creating, cleaning and modifying meshes (e.g. field projection and mesh morphing).
- Filters (in module `Containers`):
  Various types of `ElementFilter`s and `NodeFilter`s allow to handle subparts of meshes by selecting element- and node-sets using threshold functions, tags, element types, element dimensionality and masks. Filters can be combined using Boolean operations (union, complementary...).
- A finite element engine (in module `FE`):
  A general weak formulation engine able to integrate fields over parts of the meshes is available. The `FETools` submodule contains specific functions for Lagrange P1 finite elements, including the computation of stiffness and mass matrices. The domain of integration is defined using `ElementFilter`s making the integration domain flexible. P0 and P2 Lagrange finite element spaces are implemented and tested. The framework is non isoparametric: the user can write weak formulations of mixes of P0, P1 and P2 fields on P1 or P2 meshes.
- Input/Output functions (in module `IO`):
  Various readers (respectively, writers) for importing (exporting) meshes and solution fields from (to) BasicTools' internal data model are available. See [BasicTools documentation](https://basictools.readthedocs.io/en/latest/_source/BasicTools.IO.html#submodules) for the complete list of supported files. Available formats include geo/geof (Z-set [@zset]), VTK, XDMF, SAMCEF, ABAQUS, and a bridge with MeshIO is provided. ABAQUS and SAMCEF solution (fields) readers depends on local installation due to proprietary format files.
- Implicit geometry engine (in module `ImplicitGeometry`):
  The classes are used to define arbitrary subdomains using implicit geometries (level-set function). Basic shapes (spheres, half-spaces, cylinders, cubes), transformations (symmetry, translation, rotation) and binary operators (union, difference and intersection) can be used to construct complex shapes. These shapes can be used to select elements (using `ElementFilter`), or be evaluated on point clouds to explicitly construct iso-zero surfaces.
- Linear algebra functions (in module `Linalg`):
  Some common operations on linear systems for finite elements: penalization, elimination, Lagrange multipliers and the Ainsworth [@AINSWORTH20016323] method to impose essential boundary conditions or linear multi-point constraints. The submodule `LinearSolver` offers an abstraction layer for sparse linear solvers, including: Cholesky of the `sksparse` package; factorized, CG, lsqr, gmres, lgmres of the `scipy.sparse.linalg` module; CG, LU, BiCGSTAB, SPQR of the C++ Eigen library; and the AMG solver of `pyamg` package.

The large majority of functions are illustrated in the same file where they are defined, in `CheckIntegrity` functions.

# Examples

We present two examples, see [BasicTools documentation](https://basictools.readthedocs.io/en/latest/Examples.html) for complete details.

## Pre/post deep learning

Convolution-based deep learning algorithms generally rely on structured data. BasicTools is used to transfer a field computed on an unstructured mesh using finite elements to a structured grid and vice versa. To validate the operation, the error on the final field is evaluated with respect to the original field.

![Deep learning workflow coupled to finite element simulator a) Initial field on unstructured mesh, b) transferred field into regular grid (projection step) c) inverse projection into original unstructured mesh d) projection error on unstructured mesh.\label{fig:DeepLearningPrepost}](DeepLearningPrepost.png)

## Mechanical analysis: Thick plate with two inclusions

Consider a thick plate with two inclusions, one softer and the other stiffer than the base material. The plate is clamped on the left side with a negative pressure applied on the right side. We compute the strain energy on only one inclusion. The linear elasticity problem is solved using P1 Lagrange finite elements on an unstructured mesh.

![Analysis of a mechanical thick plate with two inclusions.\label{fig:TwoInclusions}](TwoInclusions_img1.png)

# References
