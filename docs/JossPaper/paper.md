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

Numerical simulations of physical phenomena can be computed by many (commercial/free) software packages, but despite the apparent variety, they rely on a relatively small set of operations during the preparation, exploitation and post-process of these simulations. For instance, handling and modifying meshes and fields are common needs. BasicTools is a Python library designed to address these supporting tasks. It features an efficient data model for meshes and field objects and input/output routines compatible with various formats. A finite element engine allows to assemble abstract variational formulations, take into account differential operators and integrate fields on volumes and surfaces.

BasicTools is actively used in artificial intelligence and model order reduction for physical problems [@ROM-net; @mca26010017; @UQindustrialDesign; @datatargetVAE], topology optimization [@nardoni], and material sciences [@pymicro] projects.

# Statement of need

Industrial design tasks often rely on numerical simulation workflows involving different softwares, each providing its own set of specific post-processing tool. Common tasks like transferring computed fields from one simulation tool to another must be routinely implemented, often with subtle variations. This limits interoperability and increases the amount of work and the risk of introducing bugs.

Our solution to these concerns is BasicTools, which introduces a data model for meshes and related physical fields. This model can be populated using different readers and exported using various writers: no new mesh or solution format is forced upon the user. The data-oriented design of BasicTools, with no low-level objects, allows high performance operations using a high level language (Python with NumPy). BasicTools contains tools to convert meshes to other "in-memory" formats (VTK [@VTK4], PyVista [@sullivan2019pyvista], MeshIO [@meshio], CGNS [@cgns], Gmsh [@gmsh] ). This make possible to mix (and reuse) treatments already available in other frameworks.
Additionally, BasicTools provides many algorithms including mesh manipulation, a finite element engine and field projection operators.


# State of the field

In the computational fluid dynamics community, the CFD General Notation System (CGNS) [@cgns] standard has been proposed for data analysis. In the solid mechanics community, to the authors' knowledge, no such de facto standard exists. Concerning meshes, one may consider MeshIO for converting between various file formats, or VTK for the manipulation of meshes, although it lacks some capabilities of tremendous value for the solid mechanics community, for example the post-process of integration point data. Most, if not all, the existent tool do the simple, but potentially dangerous approach, of extrapolating the integration point values to the nodes of the mesh or averaging in every cell. This can lead to the misinterpretation of the solution and wrong engineering decisions. Finally, only a few finite element engines allow assembling abstract variational formulations on arbitrary geometries (e.g. FreeFem++ [@freefempp] and FEniCS [@fenics]).

# Overview

BasicTools was designed with the intention of creating a set of canonical tools to work on meshes and fields.

The main features of the library are:

- Meshes (in the module `Containers`):
  `ConstantRectilinearMesh` and `UnstructuredMesh` encapsulate respectively the data model for constant rectilinear and unstructured mesh types. Unstructured meshes are efficient in the sense that elements are stored using only one array for each element type. Both mesh types can feature nodes and element tags. Many functions are available for creating, cleaning and modifying meshes. In particular, field projection operations enable to project fields defined on a mesh onto a set of points, using various methods and options, with respect to the location of the destination points being inside or outside the origin mesh (finite element interpolation, extrapolation, clamped evaluations, nearest neighbors, zero fill). Mesh morphing capabilities are also included.
- Filters (in the module `Containers`):
  Various types of `ElementFilter`s and `NodeFilter`s allow to handle subparts of the meshes by selecting element- and node-sets using threshold functions, tags, element types, element dimensionality and masks. Arbitrary filters can be combined using Boolean operations (union, complementary...) to construct advanced filters on points and elements.
- A finite element engine (in the module `FE`):
  A general weak formulation engine able to integrate fields over any part of the considered mesh is available. The `FETools` submodule contains specific functions for Lagrange P1 finite elements, including the computation of stiffness and mass matrices. The domain of integration is defined using `ElementFilter`s to make the integration domain flexible. Depending on the parameter of the integration, the result can be a matrix (e.g. tangent operator), a vector (e.g. right hand side term), or a scalar (e.g. volume, energy). Also P0 and P2 Lagrange finite element spaces are fully implemented and tested. The framework of BasicTools is non isoparametric, meaning the user can write a weak formulation of a mix of P0, P1 and P2 fields on P1 or P2 meshes.
- Input/Output functions (in the module `IO`):
  Various readers (respectively, writers) for importing (exporting) meshes and solution fields from (to) the internal data model of BasicTools are available. See [BasicTools documentation](https://basictools.readthedocs.io/en/latest/_source/BasicTools.IO.html#submodules) for the complete list of supported files. Available formats include geo/geof (Z-set [@zset]), VTK, XDMF, SAMCEF, ABAQUS. The bridge with MeshIO brings extra import/export capabilities by wrapping the MeshIO readers/writers with the BasicTools API. The ABAQUS, and the SAMCEF solution (fields) readers depends on local installation of the respective software due to the proprietary format solution files.
- Implicit geometry engine (in the module `ImplicitGeometry`):
  The classes are used to define arbitrary subdomains using only implicit geometries (level-set function). Basic shapes (spheres, half-spaces, cylinders, cubes), transformations (symmetry, translation, rotation) as well as binary operators (union, difference and intersection) can be used to construct complex shapes. Then these shapes can be used to select elements (using an `ElementFilter`), or be evaluated on a point cloud (e.g. the points of a mesh) to explicitly construct the iso-zero surface.
- Linear algebra functions (in the module `Linalg`):
  Some common operation on linear systems in the domain of finite element. In particular, the user can choose between penalization, elimination, Lagrange multipliers and the Ainsworth [@AINSWORTH20016323] method to impose essential boundary conditions or linear multi-point constraints. The submodule `LinearSolver` offers an abstraction layer for sparse linear solvers, including: Cholesky of the `sksparse` package; factorized, CG, lsqr, gmres, lgmres of the `scipy.sparse.linalg` module; CG, LU, BiCGSTAB, SPQR of the C++ Eigen library; and the AMG solver of `pyamg` package.

The large majority of functions are illustrated in the same file where they are defined, in `CheckIntegrity` functions.

# Examples

We detail two examples illustrating some of the features mentioned above.
The complete examples can be found in the [BasicTools documentation](https://basictools.readthedocs.io/en/latest/Examples.html).

## Pre/post deep learning

Convolution-based deep learning algorithms generally rely on structured data. In this example we demonstrate the use of BasicTools to transfer a field classically computed on an unstructured mesh using finite elements to a structured grid and vice versa. To validate the operation, the error on the final field is evaluated with respect to the original field.

![Example of a deep learning workflow coupled to a finite element simulator a) Initial field on a unstructured mesh, b) transferred field into a regular grid (projection step) c) inverse projection into the original unstructured mesh d) projection error on the unstructured mesh .\label{fig:DeepLearningPrepost}](DeepLearningPrepost.png)

## Mechanical analysis: Thick plate with two inclusions

We consider a thick plate with two inclusions, one softer and the other stiffer than the base material. The plate is clamped on the left side with a negative pressure applied on the right side. Then we compute the strain energy on only one inclusion. The linear elasticity problem is solved using P1 Lagrange finite elements on an unstructured mesh.

![Analysis of a mechanical thick plate with two inclusions.\label{fig:TwoInclusions}](TwoInclusions_img1.png)

# References
