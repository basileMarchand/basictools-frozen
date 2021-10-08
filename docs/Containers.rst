**********
Containers
**********

A finite element mesh is a discretisation of a domain :math:`\Omega` using different type of elements.
The next figure is an example of a mesh composed of different type of elements (quadrangles, triangles, bars and point elements).

.. figure:: images/Mesh.svg
    :width: 600px
    :align: center
    :height: 300px
    :alt: alternate text
    :figclass: align-center
    
    A finite element mesh.
    
Many real life finite element problem use mixed elements to correctly discreticide the domain.
Also to correctly represent all the quantities needed for the definition of a problem (like material properties and external force) it is necessary to define subgroups of element and points.
So a real mesh structure can not be reduced to just the points and the elements.
Also to keep track of mesh transformations (e.i. extraction, mirroring) a suitable structured is necessary. 

In general the operations on meshes (e.i. extraction, integration, export to file) can be grouped by type of elements so it is not necessary to have a very low level classes like "Point class" or "Element Class".
In BasicTools all homogeneous entities are grouped in high level classes.
This enable us to work on big chunks of data and leverage the full power of a high level language.

The :py:class:`BasicTools.Containers.UnstructuredMesh.UnstructuredMesh` class of BasicTools keep all this information in a compact compact but relatible simple structure.

    
Properties
##########


The mesh class contains a dictionary named :py:attr:`BasicTools.Containers.MeshBase.MeshBase.props` to store user data related to the mesh itself. No function must relies in the presence of a specific key in the ``props``. The values of the dictionary are copied as derivative meshed are generated. 

Tags
####

Tags are used to track subsets of points or elements in the mesh.
Each tag is defined by a string (the name) and the indices of the entities to be capture by the tag.
Internally the class :py:class:`BasicTools.Containers.MeshBase.Tag`  uses a memory buffer and a counter to reduce the number of memory allocation in the case of incrementally adding ids to the tag.

The class :py:class:`BasicTools.Containers.MeshBase.Tags` offer dictionary like interface to easy interact with a group of tags.


Nodes
#####

A :py:attr:`BasicTools.Containers.UnstructuredMesh.UnstructuredMesh.nodes` in a mesh is defined by its 3 "coordinates" (float) (2 if the mesh is a 2D mesh) and a :py:attr:`BasicTools.Containers.UnstructuredMesh.UnstructuredMesh.originalIDNodes` (int) to keep track its origin. The ``nodes`` attribute store the positions of the points (one point per row). The size of ``originalIdNodes`` must be equal to the number of points in the mesh.

The role of the originalIdNodes is double. 
In the case the origin of the mesh is a file from disk, it is used to store the "node number" (a concept present in many file formats). 
In the case the mesh is a transformation of a initial mesh (extraction for example) the originalIdNodes is used to track the origin of each entity. 


Elements
########


The storage of elements (:py:class:`BasicTools.Containers.UnstructuredMesh.ElementsContainer`) is more complex only because the information is heterogeneous (e.g. the number of point per element is different for every type of element).
The mesh in BasicTools does not store any information about the physics or any finite element implementation detail (e.i. reduce integration, degree of interpolation of the solution field).
This make the definition of the mesh relatively simple.

Each type of element is only defined by it geometrical definition (e.i. linear triangle with 3 point).

A variable named :py:attr:`BasicTools.Containers.UnstructuredMesh.UnstructuredMesh.elements`` (of type :py:class:`BasicTools.Containers.UnstructuredMesh.AllElements`), analog to a dictionary in python, store for each element type all the information.

.. figure:: images/Elements.svg
    :width: 1600px
    :align: center
    :alt: Element names and numbering 
    :figclass: align-center
    
    Elements name and numbering (the numbering is the same as vtk for compatibillity).