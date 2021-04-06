************
Introduction
************

BasicTools was primary designed as a basic set of tools to work on meshes in the context of finite element computation.
The main functionalities of the library are:

* IO routine : A set of classes to read and write meshes (and solutions fields) from/to a large variety of file formats. BasicTools does not have a proper file format mainly because existent formats provide most, if not all, of the functionalities needed.
* Mesh manipulation: Routines to filter, define, extract and manipulate meshes in many ways. 
* Fields manipulation: finite element fields can be defined using different kinds of interpolation (P0/P1/P2), in the full mesh or only in restricted zones, and also at integration points. This classes have overloaded operators to make computation of quantities of interest an easy task.
* Integration: Routines for the integration of weak formulations (tangent matrices, right hand terms, integral over only a part of a mesh). 
* Field transfer: Basic routine to transfer field from one mesh to another.
* Finite element solver: using all the previous tools, some basic finite element solver are  available to solve generic partial differential equations on non structured meshes.
