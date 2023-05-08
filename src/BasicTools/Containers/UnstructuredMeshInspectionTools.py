# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Dict, List, Optional, Tuple
import numpy as np

from BasicTools.Containers.Filters import ElementFilter
from BasicTools.FE.Fields.FEField import FEField


from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType, ArrayLike
import BasicTools.Containers.ElementNames as ElementNames
import BasicTools.Containers.Filters as Filters
from BasicTools.Containers.UnstructuredMesh import ElementsContainer, UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes
from BasicTools.NumpyDefs import PBasicIndexType


def GetDataOverALine(startPoint: ArrayLike, stopPoint: ArrayLike, nbPoints: int, mesh: UnstructuredMesh, fields: List, method: Optional[str]=None)-> UnstructuredMesh:
    """Compute the values of a mesh/field over a line
    startPoint and stopPoint are used to construct a line with nbPoints points.
    the data in mesh.nodeFields, mesh.elemFields and fields are evaluated at every point of the line

    Parameters
    ----------
    startPoint : ArrayLike
        Start point. this is passed to the CreateUniformMeshOfBars function
    stopPoint : ArrayLike
        Stop point. this is passed to the CreateUniformMeshOfBars function
    nbPoints : int
        this argument is passed to the CreateUniformMeshOfBars function
    mesh : UnstructuredMesh
        the mesh to extract mesh.nodeFields, mesh.elemFields
    fields : List
        a list containing only FEFields
    method : Optional[str], optional
        this argument is passed to the GetFieldTransferOp function, by default None

    Returns
    -------
    UnstructuredMesh
        a unstructured mesh with nodeFields populated with the data
        nodeFields with contain a dictionary for every field name and the numpy vector with
        the values at the points. If two or more fields have the same name only one is recovered,
        the priority order staring with the hightest priority are : fields, elemFields and nodeFields

    Raises
    ------
    Exception
        In the case an IPField is present in the 'fields' argument
    """
    res: Dict[str,np.ndarray] = {}
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars

    lineMesh = CreateUniformMeshOfBars(startPoint,stopPoint,nbPoints)

    from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
    from BasicTools.FE.Fields.FieldTools import NodeFieldToFEField, ElemFieldsToFEField

    #transfer of nodal fields
    nodeFieldsAsFEFields = NodeFieldToFEField(mesh)
    if len(nodeFieldsAsFEFields) != 0:
        firstField = next(iter(nodeFieldsAsFEFields.values()))
        op,status = GetFieldTransferOp(firstField, lineMesh.nodes, method=method)
        for name, data in mesh.nodeFields.items():
            res[name] = op.dot(data)


    #transfer of cell data
    cellFieldsAsFEFields = ElemFieldsToFEField(mesh)
    if len(cellFieldsAsFEFields) != 0:
        firstField = next(iter(cellFieldsAsFEFields.values()))
        op,status = GetFieldTransferOp(firstField, lineMesh.nodes, method=method)
        for name, data in mesh.elemFields.items():
            res[name] = op.dot(data)

    #transfer of FEFields
    for efField in fields:
        if isinstance(efField,FEField):
            op, status = GetFieldTransferOp(efField, lineMesh.nodes, method=method)
            res[efField.name] = op.dot(efField.data)
        else:# pragma: no cover
            raise Exception(f"Don't know how to treat field of type ({type(efField)})")

    lineMesh.nodeFields = res
    return lineMesh

def GetElementsFractionInside(field:ArrayLike, points:ArrayLike, name:str, elements:ElementsContainer, ids:ArrayLike)-> np.ndarray:
    """Compute the volume fraction for each element in ids, based in the level-set field field.
    Function valid for simplices (bars, triangles, tetrahedrons)

    Parameters
    ----------
    field : ArrayLike
        the levelset field
    points : ArrayLike
        the points of the mesh
    name : str
        the element name
    elements : ElementsContainer
        the element data
    ids : ArrayLike
        the ids to treat

    Returns
    -------
    np.ndarray
        the volume fraction for each element in ids
    """

    from scipy.spatial import Delaunay

    def triangle_lob_fraction(index, phis):
        div = phis-phis[index]
        div[index] = phis[index]
        val = phis[index]/(div)
        return np.prod(val)

    def extract_signs(phis):
        signs = np.sign(phis)
        dominant_sign = np.sign(np.sum(signs))
        if dominant_sign == -1:
            signs[signs != 1] = -1
        else:
            signs[signs != -1] = 1
        return signs

    def split_signs(signs):
        minuses = np.nonzero(signs == -1)[0]
        pluses = np.nonzero(signs == 1)[0]
        return minuses, pluses


    # helper function for tets

    def tetrahedron_volumic_fraction(phis):
        assert phis.size == 4
        assert np.count_nonzero(phis) > 0
        minuses, pluses = sort_by_sign(phis)
        if minuses.size == 0:
            return 0.0
        elif minuses.size == 1:
            return tetrahedron_volumic_fraction_dominated(minuses, pluses)
        elif minuses.size == 2:
            return tetrahedron_volumic_fraction_non_dominated(minuses, pluses)
        elif minuses.size == 3:
            return 1.0 - tetrahedron_volumic_fraction_dominated(pluses, minuses)
        else:
            return 1.0

    def tetrahedron_volumic_fraction_dominated(phi_in, phis_out):
        return np.prod(phi_in / (phi_in - phis_out))

    def tetrahedron_volumic_fraction_non_dominated(phis_in, phis_out):
        phis_in_r = np.reshape(phis_in, (2, 1))
        phis_out_r = np.reshape(phis_out, (1, 2))
        ratios = phis_in_r / (phis_in_r - phis_out_r)
        result = ratios[0, 0] * ratios[0, 1] * (1.0 - ratios[1, 0]) + \
                ratios[1, 0] * ratios[1, 1] * (1.0 - ratios[0, 1]) + \
                ratios[1, 0] * ratios[0, 1]
        return result

    def sort_by_sign(phis):
        signs = np.sign(phis)
        dominant_sign = np.sign(np.sum(signs))
        if dominant_sign == -1:
            minuses = phis[signs != 1]
            pluses = phis[signs == 1]
        else:
            minuses = phis[signs == -1]
            pluses = phis[signs != -1]
        return minuses, pluses

    #-----------------------------------------------------------------------------
    res = np.zeros(len(ids)) -1

    for cpt,id in enumerate(ids):
        coon = elements.connectivity[id,:]
        vals = field[coon]
        if np.all(vals<0)  :
            res[cpt] = 1
        elif np.all(vals>0) :
            res[cpt] = 0
        else:
            elemPoints = points[coon,:]
            simplex = Delaunay(points).simplices
            if ElementNames.dimension[name] == 1 and len(vals) == 2:
                #bar2
                if vals[0] < 0 :
                    res[cpt] = abs(vals[0])/np.sup(abs(vals))
                else:
                    res[cpt] = 1-abs(vals[0])/np.sup(abs(vals))
            elif ElementNames.dimension[name] == 2 and len(vals) == 3:
               # triangles
                signs = extract_signs(vals)
                minuses, pluses = split_signs(signs)

                if minuses.size == 2:
                    res[cpt] =  1.0 - triangle_lob_fraction(pluses[0], vals)

                elif minuses.size == 1:
                    res[cpt] =  triangle_lob_fraction(minuses[0], vals)
            elif ElementNames.dimension[name] == 3 and len(vals) == 4:
                # tets
                res[cpt] = tetrahedron_volumic_fraction(vals)

    return res

def GetVolumePerElement(inmesh: UnstructuredMesh, elementFilter: Optional[ElementFilter]=None ) -> np.ndarray:
    """Compute the volume (surface for 2D element and length for 1D elements) for each element selected by the elementFilter

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the mesh to extract elements
    elementFilter : Optional[ElementFilter], optional
        filter to select some elements, if None the volume of all the element are computed

    Returns
    -------
    np.ndarray
        a numpy array of size number of "element selected by the elementFilter" with the volume
    """
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.FE.Integration import IntegrateGeneral

    numbering = ComputeDofNumbering(inmesh,LagrangeSpaceP0,elementFilter=elementFilter)

    wform = GetField("F",1).T*GetTestField("T",1)

    F = FEField("F",inmesh,LagrangeSpaceP0,numbering)
    F.Allocate(1.)
    unkownFields = [ FEField("T",mesh=inmesh,space=LagrangeSpaceP0,numbering=numbering) ]
    _,f  = IntegrateGeneral( mesh=inmesh, wform=wform, constants={}, fields=[F], unkownFields=unkownFields,elementFilter=elementFilter)
    return f

def GetMeasure(inmesh: UnstructuredMesh, elementFilter:Optional[ElementFilter] = None):
    """Compute the measure of the mesh (volume in 3D, surface in 2D, ...)
    Only elements of the higher dimensionality are taken into account if no element filter
    is supplied

    Attention: if the element filter contain dimensionality mixed elements the result is the
    sum of the measure (volume + surface for example)

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the mesh to use for the computation

    elementFilter : ElementFilter
        element to take into account in the computation of the measure, if none
        all the element of the higher dimensionality are used.

    Returns
    -------
    PBasicFloatType
        the volume if the mesh contains 3D elements
        the surface if the mesh contains 2D elements and no 3D elements
        the length if the mesh contains 1D elements and no 3D elements nor 2D elements
    """
    if elementFilter is None:
        elementFilter = ElementFilter(inmesh,dimensionality=inmesh.GetElementsDimensionality())

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo, ConstantSpaceGlobal
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.FE.Integration import IntegrateGeneral

    numbering = ComputeDofNumbering(inmesh,LagrangeSpaceGeo,fromConnectivity=True)

    wform = GetField("F",1).T*GetTestField("T",1)

    F = FEField("F",inmesh,LagrangeSpaceGeo,numbering)
    F.Allocate(1.)
    gNumbering = ComputeDofNumbering(inmesh,ConstantSpaceGlobal)
    unkownFields = [ FEField("T",mesh=inmesh,space=ConstantSpaceGlobal,numbering=gNumbering) ]
    _,f  = IntegrateGeneral( mesh=inmesh, wform=wform, constants={}, fields=[F], unkownFields=unkownFields, elementFilter=elementFilter)
    return f[0]

def GetVolume(inmesh: UnstructuredMesh) -> PBasicFloatType:
    """Compute the volume of the mesh
    Only element of the bigger dimensionality are taken into account
    if no elements on mesh we return zero.

    for a more general function see: GetMeasure

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the mesh to use for the computation

    Returns
    -------
    PBasicFloatType
        the volume if the mesh contains 3D elements
        the surface if the mesh contains 2D elements and no 3D elements
        the length if the mesh contains 1D elements and no 3D elements nor 2D elements


    """
    return GetMeasure(inmesh)

def ComputeNodeToElementConnectivity(inmesh: UnstructuredMesh, maxNumConnections:int=200)-> Tuple[np.ndarray,np.ndarray]:
    """Compute the node to connectivity connection for every point

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the mesh
    maxNumConnections : int, optional
        the number of potential connection to a single node, normally 200 is more than enough, by default 200

    Returns
    -------
    Tuple
        numpy.ndarray, (nb point, max nb connections) item, the connectivity for each node -> element
        numpy.array,  the number of connections per point
    """
    # generation of the dual graph
    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=PBasicIndexType )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType )

    cpt =0

    for name, elements in inmesh.elements.items():
        for i in range(elements.GetNumberOfElements()):
            coon = elements.connectivity[i,:]
            for j in coon:
                dualGraph[j,usedPoints[j]] =  cpt
                usedPoints[j] += 1
            cpt += 1

    #we crop the output data
    maxsize = np.max(np.sum(dualGraph>=0,axis=1))
    dualGraph = dualGraph[:,0:maxsize]
    return dualGraph,usedPoints

def ComputeNodeToNodeConnectivity(inmesh: UnstructuredMesh, maxNumConnections=200)-> Tuple[np.ndarray,np.ndarray]:
    """Compute the node to node connection for every node in the mesh

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the mesh
    maxNumConnections : int, optional
        the number of potential connection to a single node, normally 200 is more than enough, by default 200

    Returns
    -------
    Tuple
        numpy.ndarray, (nb point, max nb connections) item, the connectivity for each node -> node
        numpy.array,  the number of connections per point
    """
    # generation of the dual graph
    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=PBasicIndexType )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType )

    for name, elements in inmesh.elements.items():
        size = elements.GetNumberOfNodesPerElement()
        for i in range(elements.GetNumberOfElements()):
            coon = elements.connectivity[i,:]
            for j in range(size):
                myIndex = coon[j]
                for k in range(size):
                    if k == j:
                        continue
                    dualGraph[myIndex,usedPoints[myIndex]] =  coon[k]
                    usedPoints[myIndex] += 1
                    # we reached the maximum number of connection
                    # normally we have some data duplicated
                    # we try to shrink the vector
                    if usedPoints[myIndex] == maxNumConnections :# pragma: no cover
                        # normally we shouldn't pas here
                        c = np.unique(dualGraph[myIndex,:])
                        dualGraph[myIndex,0:len(c)] = c
                        usedPoints[myIndex] = len(c)

    maxsize = 0
    # we finish now we try to compact the structure
    for i in range(inmesh.GetNumberOfNodes()):
        c = np.unique(dualGraph[i,0:usedPoints[i]])
        dualGraph[i,0:len(c)] = c
        dualGraph[i,len(c):] = -1
        usedPoints[i] = len(c)
        maxsize = max(len(c),maxsize)

    #we crop the output data
    dualGraph = dualGraph[:,0:maxsize]
    return dualGraph,usedPoints

def ComputeElementToElementConnectivity(inmesh:UnstructuredMesh, dimensionality:int = None, connectivityDimension:int = None, maxNumConnections:int = 200)-> Tuple[np.ndarray,np.ndarray]:
    """Generates the elements to element connectivity of an UnstructuredMesh in the following sense:
    an element is linked to another in the graph if they share:
        - a face   if connectivityDimension = 2
        - an edge  if connectivityDimension = 1
        - a vertex if connectivityDimension = 0
        (if connectivityDimension is initialized to None, it will be set to dimensionality - 1)

    Also provides the array of number of connection per cell.

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the input mesh
    dimensionality : int, optional
        dimensionality filter, by default inmesh.GetDimensionality()
    connectivityDimension : int, optional
        type of connectivity, by default dimensionality - 1
    maxNumConnections : int, optional
        the number of potential connection to a single node, normally 200 is more than enough, by default 200

    Returns
    -------
    Tuple
        numpy.ndarray, (nb elements, max nb connections) item, the connections for each element -> element
        numpy.array,  the number of connections per elements
    """


    if dimensionality == None:
        dimensionality = inmesh.GetDimensionality()

    if connectivityDimension == None:
        connectivityDimension = dimensionality - 1

    assert dimensionality > connectivityDimension

    elFilter = Filters.ElementFilter(inmesh, dimensionality = dimensionality)

    meshGraph, dump = ComputeNodeToElementConnectivity(inmesh) # build node connectivity table first
    elemGraph = np.zeros((inmesh.GetNumberOfElements(), maxNumConnections), dtype = PBasicIndexType) - 1
    usedCells = np.zeros(inmesh.GetNumberOfElements(), dtype = PBasicIndexType)

    for groupName, elemGroup, ids in elFilter:
        neighboringElems = np.zeros((elemGroup.connectivity.shape[1], meshGraph.shape[1]), dtype = PBasicIndexType)
        nbNodesPerFace = []

        if connectivityDimension == 2 or (connectivityDimension == 1 and dimensionality == 2):
            connectingObjects = ElementNames.faces[groupName]
        elif connectivityDimension == 1 and dimensionality == 3:
            connectingObjects = ElementNames.faces2[groupName]
        else: # hence connectivityDimension == 0
            connectingObjects = [(ElementNames.Point_1,[0])]

        for connectingObject in connectingObjects: # evaluate connectivity
            nbNodesPerFace.append(len(connectingObject[1]))

        for element in ids: # for each element in the filtered groups
            elementNodes = elemGroup.connectivity[element,:]
            neighboringElems.fill(0)
            for j, node in enumerate(elementNodes):
                neighboringElems[j] = meshGraph[node]
            unique, counts = np.unique(neighboringElems, return_counts=True) # count how many time each other element is connected to a node

            yindex = 0
            for j in range(len(unique)):
                condition1 = (counts[j] in nbNodesPerFace and connectivityDimension > 0)# if nb of connections correspond to face connectivity, add it to graph
                condition2 = (element!=unique[j] and j>0 and connectivityDimension == 0)
                if condition1 or condition2:
                    elemGraph[element,yindex] = unique[j]
                    usedCells[element] += 1
                    yindex += 1

    maxsize = np.max(np.sum(elemGraph>=0,axis=1)) # crop output data
    elemGraph = elemGraph[:,0:maxsize]
    return elemGraph, usedCells

def ExtractElementsByElementFilter(inmesh: UnstructuredMesh, elementFilter: ElementFilter, copy:bool=True) -> UnstructuredMesh:
    """Create a new mesh with the selected element by elementFilter
    For the moment this function make a copy of nodes

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the input mesh
    elementFilter : ElementFilter
        the ElementFilter to select the element to extract
    copy : bool
        if true this function will try to make reuse as mush as possible the memory of the input mesh.
        Warning!!! Making modification of one of the two meshes will potentially invalidate the other mesh

    Returns
    -------
    UnstructuredMesh
        the mesh with the extracted elements
    """

    inmesh.ComputeGlobalOffset()
    outMesh = type(inmesh)()
    if copy:
        outMesh.CopyProperties(inmesh)
        outMesh.nodes = np.copy(inmesh.nodes)
        outMesh.nodesTags = inmesh.nodesTags.Copy()
    else:
        outMesh.props = inmesh.props
        outMesh.nodes = inmesh.nodes
        outMesh.nodesTags = inmesh.nodesTags

    outMesh.originalIDNodes = np.arange(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType)


    elementFilter.mesh = inmesh
    for name,data,ids in elementFilter:
        if len(ids) == data.GetNumberOfElements() and copy == False:
            outElements = type(data)(name)
            outElements.connectivity = data.connectivity
            outElements.originalIds = np.arange(data.GetNumberOfElements(), dtype=PBasicIndexType)
            outElements.cpt = data.GetNumberOfElements()
            outElements.tags = data.tags
            outMesh.elements[name] = outElements
        else:
            outMesh.elements[name] = ExtractElementsByMask(data,ids)

    outMesh.PrepareForOutput()
    return outMesh

def ExtractElementsByMask(inElementContainer: ElementsContainer, mask: ArrayLike) -> ElementsContainer:
    """Create a new ElementContainer with the element selected by the mask

    Parameters
    ----------
    inelems : ElementsContainer
        _description_
    mask : ArrayLike
        a vector of bool (a boolean mask) or a vector with the indices to extract

    Returns
    -------
    ElementsContainer
        a new container with the extracted elements (the tags are updated)
    """

    outElements = ElementsContainer(inElementContainer.elementType)

    newIndex = np.empty(inElementContainer.GetNumberOfElements(),dtype=PBasicIndexType)

    mask = np.asarray(mask)
    if mask.dtype == bool:
        nbElements =0
        for i in range(inElementContainer.GetNumberOfElements()):
           newIndex[i] = nbElements
           nbElements += 1 if mask[i] else 0
        iMask = mask
    else:
        nbElements = len(mask)
        iMask = np.zeros(inElementContainer.GetNumberOfElements(),dtype=bool)
        cpt =0
        for index in mask:
           newIndex[index ] = cpt
           iMask[index] = True
           cpt += 1

    outElements.Allocate(nbElements)
    outElements.connectivity = inElementContainer.connectivity[iMask,:]
    outElements.originalIds = np.where(iMask)[0]

    for tag in inElementContainer.tags  :
       temp = np.extract(iMask[tag.GetIds()],tag.GetIds())
       newid = newIndex[temp]
       outElements.tags.CreateTag(tag.name).SetIds(newid)

    return outElements

def ExtractElementByTags(inmesh: UnstructuredMesh, tagsToKeep:List[str], allNodes:bool=False, dimensionalityFilter:int= None, cleanLonelyNodes:bool=True)-> UnstructuredMesh:
    """Extract element by tags. create a new UnstructuredMesh with element present in the tags

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the input mesh
    tagsToKeep : List[str]
        List of Tag name to extract element can contain tag name for nodes tags
    allNodes : bool, optional
        if tagsToKeep point to a node tags then this option controls if all the nodes of the element must be in the tag to select the element , by default False
    dimensionalityFilter : int, optional
        keep all the element of dimensionality dimensionalityFilter, by default None
    cleanLonelyNodes : bool, optional
        true to apply a after the extraction, by default True

    Returns
    -------
    UnstructuredMesh
        a new UnstructuredMesh with the extracted elements (the tags are updated)
    """


    outMesh = UnstructuredMesh()
    outMesh.CopyProperties(inmesh)

    outMesh.nodes = np.copy(inmesh.nodes)

    import copy
    for tag in inmesh.nodesTags:
        outMesh.nodesTags.AddTag(copy.deepcopy(tag) )


    nodalMask = np.zeros(inmesh.GetNumberOfNodes(),dtype = bool)
    for name, elements in inmesh.elements.items():

       #if dimensionalityFilter is not None:
       #    if dimensionalityFilter !=  ElementNames.dimension[name]:
       #        continue

       if (np.any([x in elements.tags for x in tagsToKeep] ) == False) and (dimensionalityFilter is None) :
           if np.any([x in inmesh.nodesTags for x in tagsToKeep]) == False:
               continue# pragma: no cover


       toKeep = np.zeros(elements.GetNumberOfElements(), dtype=bool)
       # check elements tags
       for tagToKeep in tagsToKeep:
           if tagToKeep in elements.tags:
               toKeep[elements.tags[tagToKeep].GetIds()] = True

       # check for nodes tags
       for tagToKeep in tagsToKeep:
           if tagToKeep in inmesh.nodesTags:
             nodalMask.fill(False)
             tag = inmesh.GetNodalTag(tagToKeep)
             nodalMask[tag.GetIds()] = True
             elemMask = np.sum(nodalMask[elements.connectivity],axis=1)
             if allNodes :
                 toKeep[elemMask == elements.GetNumberOfNodesPerElement()] = True
             else:
                 toKeep[elemMask > 0] = True

       # if dimensionality is ok and no tag to keep, we keep all elements
       if dimensionalityFilter is not None:
           if dimensionalityFilter ==  ElementNames.dimension[name]:
               toKeep[:] = True
           #if len(tagsToKeep)  == 0:
           #    toKeep[:] = True

       newIndex = np.empty(elements.GetNumberOfElements(), dtype=PBasicIndexType )
       cpt =0
       for i in range(elements.GetNumberOfElements()):
           newIndex[i] = cpt
           cpt += 1 if toKeep[i] else 0





       outElements = outMesh.GetElementsOfType(name)
       nbToKeep = np.sum(toKeep)
       outElements.Allocate(nbToKeep)
       outElements.connectivity = elements.connectivity[toKeep,:]
       outElements.originalIds = np.where(toKeep)[0]

       for tag in elements.tags  :
           temp = np.extract(toKeep[tag.GetIds()],tag.GetIds())
           newId = newIndex[temp]
           outElements.tags.CreateTag(tag.name,errorIfAlreadyCreated=False).SetIds(newId)

    if cleanLonelyNodes:
        outMesh.originalIDNodes = np.arange(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType)
        CleanLonelyNodes(outMesh)
    else:
        outMesh.originalIDNodes = np.copy(inmesh.originalIDNodes)

    outMesh.PrepareForOutput()
    return outMesh

def ComputeMeshDensityAtNodes(mesh: UnstructuredMesh)-> np.ndarray:
    """Function to compute the mesh size at each point
    This is done using the length of the bars of each element and averaging the data at nodes.
    All the element are treated, even the lower dimensionality element (not the 0D elements of course).
    This means the triangles in the surface can affect the values of the nodes on the surface
    No weighting is done.

    Parameters
    ----------
    mesh : UnstructuredMesh
        The input mesh

    Returns
    -------
    np.ndarray
        a numpy array (node field) with average of the mesh size at each point
    """
    nbNodes = mesh.GetNumberOfNodes()
    result = np.zeros(nbNodes,dtype=PBasicFloatType)
    cpt = np.zeros(nbNodes,dtype=PBasicIndexType)

    posOfNodes = mesh.GetPosOfNodes()

    def AddBarSizesToOutput(connectivityOfBars):

        dx = posOfNodes[connectivityOfBars[:,0],0] - posOfNodes[connectivityOfBars[:,1],0]
        dy = posOfNodes[connectivityOfBars[:,0],1] - posOfNodes[connectivityOfBars[:,1],1]
        if posOfNodes.shape[1] == 3:
            dz = posOfNodes[connectivityOfBars[:,0],2] - posOfNodes[connectivityOfBars[:,1],2]
        else:
            dz = 0
        lengths  = np.sqrt(dx**2+dy**2+dz**2)
        for i in range(2):
            result[connectivityOfBars[:,i]] += lengths
            cpt[connectivityOfBars[:,i]] += 1


    for name, data, ids in ElementFilter(mesh=mesh,dimensionality=3):
        faces = ElementNames.faces2[name]
        for fname, facesConnectivity in faces:
            AddBarSizesToOutput(data.connectivity[:,facesConnectivity][:,0:2])

    for name, data, ids in ElementFilter(mesh=mesh,dimensionality=2):
        faces = ElementNames.faces[name]
        for fname, facesConnectivity in faces:
            AddBarSizesToOutput(data.connectivity[:,facesConnectivity][:,0:2])

    for name, data, ids in ElementFilter(mesh=mesh,dimensionality=1):
        AddBarSizesToOutput(data.connectivity[:,0:2])

    cpt[cpt==0] = 1
    result /= cpt

    return result

def EnsureUniquenessElements(mesh: UnstructuredMesh)->None:
    """Ensure that every element if present only once on the mesh

    Parameters
    ----------
    mesh : UnstructuredMesh
        input mesh

    Raises
    ------
    Exception
        if 2 element (event if the connectivity is permuted) a exception is raised
    """
    cpt = 0
    for name,data in mesh.elements.items():
        dd = dict()
        for el in range(data.GetNumberOfElements()):
            n = tuple(np.sort(data.connectivity[el,:]))
            if len(n) != len(np.unique(n)):
                raise Exception("("+str(name)+") element " + str(el) +" (global["+str(el+cpt)+"])" + " use a point more than once ("+str(data.connectivity[el,:])+")"  )
            if n in dd.keys():
                raise Exception("("+str(name)+") element " + str(el) +" (global["+str(el+cpt)+"])" + " is a duplication of element " + str(dd[n]) +" (global["+str(dd[n]+cpt)+"])" )
            dd[n] = el
        cpt = data.GetNumberOfElements()

def MeshQualityAspectRatioBeta(mesh: UnstructuredMesh):
    """experimental mesh quality only available for tets

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh

    Raises
    ------
    Exception
        raise if the quality is larger than 1000
    """

    #https://cubit.sandia.gov/public/15.2/help_manual/WebHelp/mesh_generation/mesh_quality_assessment/tetrahedral_metrics.htm
    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] == 0:
            pass
        elif ElementNames.dimension[name] == 1:
            pass
        elif ElementNames.dimension[name] == 2:
            pass
        elif ElementNames.dimension[name] == 3:
            if name == ElementNames.Tetrahedron_4:
               mMax = 0
               for el in range(data.GetNumberOfElements()):
                   n = data.connectivity[el,:]
                   nodes = mesh.nodes[n,:]
                   p0 = nodes[0,:]
                   p1 = nodes[1,:]
                   p2 = nodes[2,:]
                   p3 = nodes[3,:]

                   a = np.linalg.norm(p1-p0)
                   b = np.linalg.norm(p2-p0)
                   normal = np.cross(p1-p0,p2-p0)
                   #https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
                   base_area = 0.5*a*b*np.sqrt(1-(np.dot(p1-p0,p2-p0)/(a*b))**2)

                   normal /= np.linalg.norm(normal)

                   # distance to the opposite point
                   d = np.dot(normal,p3-p0)

                   volume = base_area*d
                   inscribed_sphere_radius = d/3
                   #https://math.stackexchange.com/questions/2820212/circumradius-of-a-tetrahedron

                   c = np.linalg.norm(p3-p0)

                   A = np.linalg.norm(p3-p2)
                   B = np.linalg.norm(p1-p3)
                   C = np.linalg.norm(p2-p1)

                   circumRadius_sphere_radius = np.sqrt((a*A+b*B+c*C)*
                                                         (a*A+b*B-c*C)*
                                                         (a*A-b*B+c*C)*
                                                         (-a*A+b*B+c*C))/(24*volume)

                   AspectRatioBeta = circumRadius_sphere_radius/(3.0 * inscribed_sphere_radius)

                   mMax = max([mMax,AspectRatioBeta])
                   if AspectRatioBeta > 1000:
                       raise Exception("Element " +str(el) + " has quality of " +str(AspectRatioBeta))
        else:
            raise

def ComputeMeshMinMaxLengthScale(mesh) -> Tuple[PBasicFloatType,PBasicFloatType]:
    """Compute a estimation of the minimal and maximal length scale of the elements
    for this we compute the centroid of the element and compute the min and max distance
    between the centroid and each node of the mesh.
    Then return a tuple with (2*min,2*max) of all the distances on the mesh


    Parameters
    ----------
    mesh : UnstructuredMesh
        The input mesh

    Returns
    -------
    Tuple[PBasicFloatType,PBasicFloatType]
        tuple with (2*min_dist,2*max_dist) for all the elements in the mesh
    """
    (resMin,resMax) = (None,None)
    for name,data in mesh.elements.items():
        if data.GetNumberOfNodesPerElement() < 2: continue
        if data.GetNumberOfElements() == 0: continue
        posX = mesh.nodes[data.connectivity,0]
        posY = mesh.nodes[data.connectivity,1]
        meanX = np.sum(posX,axis=1)/data.GetNumberOfNodesPerElement()
        meanY = np.sum(posY,axis=1)/data.GetNumberOfNodesPerElement()
        meanX.shape = (len(meanX),1)
        meanY.shape = (len(meanX),1)

        distToBarycenter2 = (posX-meanX)**2 + (posY-meanY)**2

        if mesh.nodes.shape[1] == 3:
            posZ = mesh.nodes[data.connectivity,2]
            meanZ = np.sum(posZ,axis=1)/data.GetNumberOfNodesPerElement()
            meanZ.shape = (len(meanX),1)
            distToBarycenter2 += (posZ-meanZ)**2

        distToBarycenter = np.sqrt(distToBarycenter2)
        mMin = np.min(distToBarycenter)
        mMax = np.max(distToBarycenter)
        if resMin is None:
            resMin = mMin
        else:
            resMin = min(resMin,mMin)

        if resMax is None:
            resMax = mMax
        else:
            resMax = min(resMax,mMax)
    return (2*resMin,2*resMax)

def PrintMeshInformation(mesh: UnstructuredMesh):
    """Print mesh information to the screen

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh

    """

    from BasicTools.Helpers.TextFormatHelper import TFormat as TF

    def L25(text):
        return TF.Left(text ,fill=" ",width=25)
    def L15(text):
        return TF.Left(text ,fill=" ",width=15)

    mesh.VerifyIntegrity()

    print(TF.Center("Nodes information"))
    print(L25("nodes.shape:") , str(mesh.nodes.shape) )
    print(L25("nodes.dtype:") , str(mesh.nodes.dtype) )
    mesh.ComputeBoundingBox()
    print(L25("boundingMin:") + str(mesh.boundingMin) )
    print(L25("boundingMax:") + str(mesh.boundingMax) )
    print(L25("originalIDNodes.shape:") , str(mesh.originalIDNodes.shape) )
    print(L25("min(originalIDNodes):") , str(min(mesh.originalIDNodes) )  )
    print(L25("max(originalIDNodes):") , str(max(mesh.originalIDNodes) )  )
    print("---- nodes tags ---- " )
    print(mesh.nodesTags)
    for tag in mesh.nodesTags:
        res = L15(tag.name) + L15(" size:"+ str(len(tag)))
        if len(tag):
            res +=L15(" min:"+str(min(tag.GetIds()) ) )+ "  max:"+str(max(tag.GetIds()) )
        print(res)


    print(TF.Center("Elements information"))

    for name,data in mesh.elements.items():
        res = L25(str(type(data)).split("'")[1].split(".")[-1]+ " " + name+" ")
        res += L15("size:" + str(data.GetNumberOfElements()))
        res += L25("min(connectivity):" + str(min(data.connectivity.ravel())))
        res += L25("max(connectivity):" +  str(max(data.connectivity.ravel())))
        print(res)
        for tag in data.tags:
            res = L15(tag.name) + L15(" size:"+ str(len(tag)))
            if len(tag):
                res +=L15(" min:"+str(min(tag.GetIds()) ) )+ "  max:"+str(max(tag.GetIds()) )
            print("  Tag: " +res)

    print(TF.Center(""))
#------------------------- CheckIntegrity ------------------------
def CheckIntegrity_ExtractElementByTags(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.AddElementToTagUsingOriginalId(0,"first")
    res.AddElementToTagUsingOriginalId(1,"second")
    res.AddElementToTagUsingOriginalId(0,"all")
    res.AddElementToTagUsingOriginalId(1,"all")
    res.GetNodalTag("Point3").AddToTag(3)

    if ExtractElementByTags(res,["first"] ).GetNumberOfElements() != 1:
        raise# pragma: no cover
    if ExtractElementByTags(res,["Point3"] ).GetNumberOfElements() != 1:
        raise# pragma: no cover
    if ExtractElementByTags(res,["Point3"],allNodes=True ).GetNumberOfElements() != 0:
        raise# pragma: no cover

    if ExtractElementByTags(res,[],dimensionalityFilter=2 ).GetNumberOfElements() != 2:
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_GetVolume(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf,CreateMeshFromConstantRectilinearMesh

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    vol = GetVolume(mesh)
    if vol != (1./6.):
        raise Exception('Error en the calculation of the volume')# pragma: no cover

    vol1elem = GetVolumePerElement(mesh,ElementFilter(mesh,elementType=ElementNames.Tetrahedron_4))
    print(vol1elem)
    if sum(vol1elem) != vol:
        raise Exception("incompatible solution of GetVolume vs GetVolumePerElement")

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([3,3,3])
    myMesh.SetSpacing([0.5, 0.5,0.5])
    print(myMesh)
    vol  = GetVolume(CreateMeshFromConstantRectilinearMesh(myMesh))
    if abs(vol-1.) > 1e-8 :
        raise Exception('Error en the calculation of the volume')# pragma: no cover

    return "ok"

def CheckIntegrity_EnsureUniquenessElements(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    # this function must work
    EnsureUniquenessElements(mesh)

    tets = [[0,1,2,3],[3,1,0,2]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    # This function must not work
    try:
       EnsureUniquenessElements(mesh)
       return "No ok"
    except :
       pass
    return "ok"

def CheckIntegrity_GetElementsFractionInside(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles,CreateMeshOf

    meshTriangles = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    field = np.array([-1, -1, -1, 1])
    for name,elements in meshTriangles.elements.items():
        res = GetElementsFractionInside(field,meshTriangles.nodes,name,elements,range(elements.GetNumberOfElements()))
        print(res)

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh3D = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    for name,elements in mesh3D.elements.items():
        res = GetElementsFractionInside(field,mesh3D.nodes,name,elements,range(elements.GetNumberOfElements()))
        print(res)

    return "ok"

def CheckIntegrity_ComputeNodeToNodeConnectivity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    dg, unUsed = ComputeNodeToNodeConnectivity(res)

    return "ok"

def CheckIntegrity_ComputeMeshMinMaxLengthScale(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],[3,1,0,2]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    print(ComputeMeshMinMaxLengthScale(mesh))
    PrintMeshInformation(mesh)
    return "ok"

def CheckIntegrity_ExtractElementsByMask(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("tri1").AddToTag(0)

    #tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([False,True]))
    #print(tri.connectivity)
    print(res.GetElementsOfType(ElementNames.Triangle_3).connectivity)

    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0],dtype=PBasicIndexType))
    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0,1],dtype=bool))

    print(tri.connectivity)
    print(tri.originalIds)
    return "ok"

def CheckIntegrity_MeshQualityAspectRatioBeta(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh  = CreateCube(dimensions = [20,20,20], origin = [0.1,0.1,0.,], spacing=[0.9/19]*3)
    MeshQualityAspectRatioBeta(mesh)
    return "ok"

def CheckIntegrity_GetDataOverALine(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh  = CreateCube(dimensions = [20,20,20], origin = [0.1,0.1,0.,], spacing=[0.9/19]*3)
    mesh.nodeFields["xpos"] = mesh.nodes[:,0]

    from BasicTools.FE.Fields.FieldTools import CreateFieldFromDescription
    field = CreateFieldFromDescription(mesh,[(ElementFilter(mesh,zone = lambda x :  x[:,2]>0.5),1)])
    field.name = "biMaterial"
    mesh.elemFields["biMaterialAtElem"] = field.GetCellRepresentation()

    print("---")
    print(mesh.elemFields["biMaterialAtElem"])
    print("---")

    res = GetDataOverALine([0,0,0], [1,1,1], 100, mesh, [field])
    print(res)
    if res.GetNumberOfNodes()  != 100:
        raise Exception("Wrong number of point in the line")

    if len(res.nodeFields) != 3:
        raise Exception("Wrong number of fields")

    if GUI:
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh,filename="CheckIntegrity_GetDataOverALine_bulk.xdmf",run =False)
        OpenInParaView(res,filename="CheckIntegrity_GetDataOverALine_line.xdmf")

    return "ok"

def CheckIntegrity(GUI=False):
    toTest= [
    CheckIntegrity_GetDataOverALine,
    CheckIntegrity_MeshQualityAspectRatioBeta,
    CheckIntegrity_EnsureUniquenessElements,
    CheckIntegrity_ExtractElementsByMask,
    CheckIntegrity_GetVolume,
    CheckIntegrity_ComputeNodeToNodeConnectivity,
    CheckIntegrity_ComputeMeshMinMaxLengthScale,
    CheckIntegrity_GetElementsFractionInside,
    CheckIntegrity_ExtractElementByTags,
    ]

    from BasicTools.Helpers.Tests import RunListOfCheckIntegrities
    return RunListOfCheckIntegrities(toTest, GUI)

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover