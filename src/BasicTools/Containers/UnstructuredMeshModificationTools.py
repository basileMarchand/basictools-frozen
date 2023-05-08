# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import copy
from typing import Optional, Tuple
import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.MeshBase import Tags
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType, ArrayLike

def CleanDoubleNodes(mesh: UnstructuredMesh, tol: Optional[PBasicFloatType]=None, nodesToTestMask: ArrayLike=None):
    """Remove double nodes for the input mesh

    Parameters
    ----------
    mesh : UnstructuredMesh
        the in put mesh
    tol : Optional[PBasicFloatType], optional
        the tolerance, by default value is = np.linalg.norm(mesh.boundingMax - mesh.boundingMin)*1e-7
        if tol is zero a faster algorithm is used
    nodesToTestMask : ArrayLike, optional
        a mask of len number of nodes with the value True for nodes to be tested (this can increase the speed of the algorithm), by default None

    """

    mesh.ComputeBoundingBox()
    if tol is None:
        tol = np.linalg.norm(mesh.boundingMax - mesh.boundingMin)*1e-7

    nbNodes = mesh.GetNumberOfNodes()
    toKeep = np.zeros(nbNodes, dtype=bool )
    newIndex = np.zeros(nbNodes, dtype=PBasicIndexType)
    tol2 = tol**2
    def dist2(array,value):
        d = array-value
        return np.inner(d,d)

    if nodesToTestMask is None and tol == 0:
        # optimized version for tol = 0.0
        database = {}
        cpt  = 0
        for i in range(nbNodes):
            point = tuple(mesh.nodes[i,:])
            ni = database.get(point)
            if ni is None:
                database[point] = cpt
                newIndex[i] = cpt
                cpt +=1
                toKeep[i] = True
            else:
                newIndex[i] = ni

    elif nodesToTestMask is None :

        from BasicTools.Containers.Octree import Octree
        cpt  = 0

        Ma = mesh.boundingMax
        Mi = mesh.boundingMin
        diag2 = np.linalg.norm(mesh.boundingMax-mesh.boundingMin)**2
        tree = Octree(Ma[0],Ma[1],Ma[2], Mi[0], Mi[1], Mi[2])

        for i in range(nbNodes):
            point = tuple(mesh.nodes[i,:])
            entries = tree.find_within_range_cube(point, tol)

            oDist = diag2
            if len(entries):
                for entry in entries:
                    dist = (point[0] - entry[0][0])**2   +(point[1] - entry[0][1])**2+(point[2] - entry[0][2])**2
                    if dist < oDist:
                        oDist = dist
                        index = entry[1]

            if oDist <= tol2:
                newIndex[i] = index
            else:
                tree.add_item(cpt, point )
                newIndex[i] = cpt
                cpt += 1
                toKeep[i] = True


    else:
        cpt =0
        for i in range(nbNodes):
            if not nodesToTestMask[i]:
                newIndex[i] = cpt
                cpt += 1
                toKeep[i] = True
                continue

            posI = mesh.nodes[i,:]

            for j in range(i):
                if not nodesToTestMask[j]:
                    continue
                if toKeep[j]:
                    dist = dist2(posI,mesh.nodes[j,:] )
                    if dist < tol2 :
                        newIndex[i] = newIndex[j]
                        break
            else:
                newIndex[i] = cpt
                cpt += 1
                toKeep[i] = True

    mesh.nodes = mesh.nodes[toKeep,:]
    mesh.originalIDNodes = np.where(toKeep)[0]

    for tag in mesh.nodesTags :
        tag.SetIds(np.unique(newIndex[tag.GetIds()]) )

    for elementName in mesh.elements:
        elements = mesh.elements[elementName]
        elements.connectivity = newIndex[elements.connectivity]

    mesh.PrepareForOutput()

def CleanLonelyNodes(mesh:UnstructuredMesh, inPlace:bool=True)-> Tuple[np.ndarray, UnstructuredMesh]:
    """Remove nodes not used by the elements

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh
    inPlace : bool, optional
        if true the nodes are removed in place, by default True

    Returns
    -------
    (usedNodes, cleanMesh)
        usedNodes: a mask of type ndarray (size the number of nodes in mesh)  with the nodes present on the cleanMesh
        cleanMesh: UnstructuredMesh the mesh without lonely nodes
    """
    usedNodes = np.zeros(mesh.GetNumberOfNodes(),dtype=bool )
    for data in mesh.elements.values():
        usedNodes[data.connectivity.ravel()] = True

    cpt = 0
    NewIndex =  np.zeros(mesh.GetNumberOfNodes(),dtype=PBasicIndexType )-1
    for n in range(mesh.GetNumberOfNodes()):
        if usedNodes[n]:
            NewIndex[n] = cpt
            cpt += 1
    originalIDNodes = np.where(usedNodes)[0]
    #filter the nodes
    if inPlace:
        mesh.nodes = mesh.nodes[usedNodes ,:]
        mesh.originalIDNodes = mesh.originalIDNodes[originalIDNodes]
        newTags = Tags()
        #node tags
        for tag in mesh.nodesTags :
            newTags.CreateTag(tag.name).SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])
        mesh.nodesTags = newTags

        #renumbering the connectivity matrix
        for elementName in mesh.elements:
            elements = mesh.elements[elementName]
            elements.connectivity = NewIndex[elements.connectivity]

        for name,data in mesh.nodeFields.items():
            mesh.nodeFields[name] = mesh.nodeFields[name][usedNodes]
        return usedNodes, mesh
    else:
        res = UnstructuredMesh()
        res.nodes = mesh.nodes[usedNodes ,:]
        res.originalIDNodes = originalIDNodes
        #node tags
        for tag in mesh.nodesTags :
            outTag = res.GetNodalTag(tag.name)
            outTag.SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])

        for elementName in mesh.elements:
            elements = mesh.elements[elementName]
            outElements = res.GetElementsOfType(elementName)
            outElements.connectivity = NewIndex[elements.connectivity]

            outElements.tags = elements.tags.Copy()

        for name,data in mesh.nodeFields.items():
            mesh.nodeFields[name] = mesh.nodeFields[name][usedNodes]
        return usedNodes, res

def CleanDoubleElements(mesh: UnstructuredMesh):
    """Remove double elements on the mesh (inplace), even if a permutation exist

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh
    """
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    elemFieldMask = np.zeros((mesh.GetNumberOfElements(),), dtype=PBasicIndexType)
    maskCpt = 0 # mask counter
    elementCpt = 0 # element counter
    newElementContainer = type(mesh.elements)()
    for elemName, data in mesh.elements.items():
        if data.mutable == False:
            newElementContainer[elemName] = data
            elemFieldMask[maskCpt:maskCpt+data.GetNumberOfElements()] = np.arange(data.GetNumberOfElements())+elementCpt
            maskCpt  += data.GetNumberOfElements()
            elementCpt += data.GetNumberOfElements()
            continue

        nbe = data.GetNumberOfElements()
        if  nbe == 0:
            continue
        data.tighten()
        _, index, = np.unique(np.sort(data.connectivity,axis=1),return_index=True,axis=0)
        mask = np.zeros(nbe, dtype=bool)
        mask[index] = True
        newElements = ExtractElementsByMask(data, mask)
        newElementContainer[elemName] = newElements
        elemFieldMask[maskCpt:maskCpt+newElements.GetNumberOfElements()] = newElements.originalIds+elementCpt
        newElementContainer[elemName].originalIds = data.originalIds[newElementContainer[elemName].originalIds]

        maskCpt  += newElements.GetNumberOfElements()
        elementCpt += data.GetNumberOfElements()


    mesh.elements = newElementContainer
    elemFieldMask = elemFieldMask[0:maskCpt]
    for name, data  in  mesh.elemFields.items():
        if len(mesh.elemFields[name].shape) == 1:
            mesh.elemFields[name] = mesh.elemFields[name][elemFieldMask]
        else:
            mesh.elemFields[name] = mesh.elemFields[name][elemFieldMask,:]

def CopyElementTags(sourceMesh: UnstructuredMesh, targetMesh: UnstructuredMesh, extendTags: bool=False):
    """Copy tags from sourceMesh to the targetMesh
    We use the connectivity to identify the elements

    Parameters
    ----------
    sourceMesh : UnstructuredMesh
        The source
    targetMesh : UnstructuredMesh
        _description_
    extendTags : bool, optional
        if False will overwrite the tags on targetMesh
        if  True will extend the tags on targetMesh
        by default False
    """

    for elemName, data1 in sourceMesh.elements.items():
        data1.tighten()

        nbe1 = data1.GetNumberOfElements()
        if nbe1 == 0 :
            continue
        if elemName not in targetMesh.elements :
            continue
        data2 = targetMesh.elements[elemName]
        data2.tighten()
        nbe2 = data2.GetNumberOfElements()
        if nbe2 == 0 :
            continue

        connectivity = np.vstack((data2.connectivity,data1.connectivity))
        _,index  = np.unique(np.sort(connectivity,axis=1),return_inverse=True,axis=0)
        ids = index[nbe2:]
        for tag1 in data1.tags:
            id2 = ids[tag1.GetIds()]
            mask = np.zeros(_.shape[0])
            mask[id2] = True
            id2 = np.where(mask[index[:nbe2]])[0]

            tag2 = data2.tags.CreateTag(tag1.name,False)
            if extendTags:
                tag2.AddToTag(id2)
            else:
                tag2.SetIds(id2)

            tag2.RemoveDoubles()

def DeleteElements(mesh: UnstructuredMesh, mask: ArrayLike):
    """Delete elements on the input mesh (inplace)

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh
    mask : ArrayLike
        the mash of elements to delete
    """
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    OriginalNumberOfElements = mesh.GetNumberOfElements()
    dataToChange = {}

    offsets = mesh.ComputeGlobalOffset()

    for name, data in mesh.elements.items():
        offset = offsets[name]
        localMask = mask[offset:offset+data.GetNumberOfElements()]
        dataToChange[name] = ExtractElementsByMask(data, np.logical_not(localMask))

    for name, data in dataToChange.items():
        mesh.elements[name] = data

    if OriginalNumberOfElements != mesh.GetNumberOfElements():
        print("Number Of Elements Changed: Dropping elemFields")
        mesh.elemFields = {}

    mesh.PrepareForOutput()

def DeleteInternalFaces(mesh:UnstructuredMesh):
    """Delete faces not present on the skin (internal faces) (inplace)

    Parameters
    ----------
    mesh : UnstructuredMesh
        the mesh

    """
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    OriginalNumberOfElements = mesh.GetNumberOfElements()
    skin = ComputeSkin(mesh)
    dataToChange = {}

    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] != 2:
            continue
        ne = data.GetNumberOfElements()
        mask = np.zeros(ne, dtype=bool)

        data2 = skin.elements[name]
        ne2 =data2.GetNumberOfElements()
        surf2 = {}
        key = np.array([ne**x for x in range(ElementNames.numberOfNodes[name]) ])

        for i in range(ne2):
            cc = data2.connectivity[i,:]
            lc = np.sort(cc)
            elementHash = np.sum(lc*key)
            surf2[elementHash] = [1,cc]

        for i in range(ne):
            cc = data.connectivity[i,:]
            lc = np.sort(cc)

            elementHash = np.sum(lc*key)
            if elementHash in surf2:
                mask[i] = True

        if np.any(mask):
            dataToChange[name] = ExtractElementsByMask(data,mask)

    for name, data in dataToChange.items():
        mesh.elements[name] = data

    if OriginalNumberOfElements != mesh.GetNumberOfElements():
        print("Number Of Elements Changed: Dropping elemFields")
        mesh.elemFields = {}

    mesh.PrepareForOutput()

def ComputeSkin(mesh:UnstructuredMesh, md:Optional[int]=None, inPlace:bool=False)->UnstructuredMesh:
    """Compute the skin of a mesh (mesh), if md (mesh dimensionality) is None
    the mesh.GetDimensionality() is used to filter the element to compute
    the skin.

    Warning if some elements are duplicated the behavior of computeSkin with
    the option inPlace=True is not defined (Use CleanDoubleElements before
    ComputeSkin)

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh
    md : Optional[int], optional
        _description_, by default None
    inPlace : bool, optional
        if True the skin is added to the original mesh, a tag "Skin" is created.
        if False a new mesh is returned with the elements of the skin (even if
            some element of the skin are already present on the original mesh).
            If the user merged the two meshes, you must call CleanDoubleElements to clean
            the mesh
        , by default False

    Returns
    -------
    UnstructuredMesh
        The mesh containing the skin, (the mesh if inPlace == True, or only the skin if inPlace ==False)
    """

    if md is None:
        md = mesh.GetDimensionality()

    # first we count the total number of potential individual skin elemnts
    totalCpt = {}
    DFilter = ElementFilter(mesh,dimensionality = md)
    for elemName, data, ids in DFilter:
        faces = ElementNames.faces[elemName]
        nbElements = data.GetNumberOfElements()
        for faceType,localFaceConnectivity in faces:
            cpt = totalCpt.setdefault(faceType,0) + nbElements
            totalCpt[faceType] =  cpt

    if inPlace :
        # we add the element already present on the mesh
        D1filter = ElementFilter(mesh,dimensionality = md-1)
        for elemName, data ,ids  in D1filter:
            nbElements = data.GetNumberOfElements()
            cpt = totalCpt.setdefault(elemName,0) + nbElements
            totalCpt[elemName] =  cpt

    #  Allocation of matrices to store all skin elements
    storage = {}
    for k,v in totalCpt.items():
        nn = ElementNames.numberOfNodes[k]
        storage[k] = np.empty((v,nn),dtype=int)
        totalCpt[k] = 0

    if inPlace :
        # fill storage with skin element already on the mesh
        for elemName, data, ids in D1filter:
            cpt = totalCpt[elemName]
            nbElements = data.GetNumberOfElements()
            storage[elemName][cpt:cpt+nbElements,:] = data.connectivity
            totalCpt[elemName] += nbElements

    # fill storage with skin of elements of dimensionality D
    for elemName,data,ids in DFilter:
        faces = ElementNames.faces[elemName]
        nbElements = data.GetNumberOfElements()
        for faceType,localFaceConnectivity in faces:
            globalFaceConnectivity = data.connectivity[:,localFaceConnectivity]
            cpt = totalCpt[faceType]
            storage[faceType][cpt:cpt+nbElements,:] = globalFaceConnectivity
            totalCpt[faceType] += nbElements

    if inPlace:
        res = mesh
    else:
        res = UnstructuredMesh()
        res.nodesTags =  mesh.nodesTags
        res.originalIDNodes =  mesh.originalIDNodes
        res.nodes = mesh.nodes
        res.elements = type(mesh.elements)()

    # recover only the unique elements
    for elemName,cpt in totalCpt.items():
        if cpt == 0:
            continue
        store = storage[elemName]
        _,index,counts = np.unique(np.sort(store,axis=1),return_index=True,return_counts=True,axis=0)

        # the elements
        if inPlace:
            uniqueElements = index[counts==2]
            if elemName in mesh.elements:
                meshSkinElements = mesh.elements.GetElementsOfType(elemName)
                nbMElements = meshSkinElements.GetNumberOfElements()
                ids =  uniqueElements[uniqueElements < nbMElements]
                # tag element already present on the original mesh
                meshSkinElements.tags.CreateTag("Skin",False).SetIds(ids)
            else:
                nbMElements = 0
        else:
            nbMElements = 0

        # all the elements present only 1 time not in the original mesh
        uniqueElements = index[counts==1]
        ids =  uniqueElements[uniqueElements >= nbMElements]
        if len(ids) == 0 :
            continue

        if inPlace:
            elementContainer = mesh.elements.GetElementsOfType(elemName)
        else:
            elementContainer = res.GetElementsOfType(elemName)


        newElements = store[ids,:]
        # add and tag new elements into the output
        elementContainer.tags.CreateTag("Skin",False).AddToTag(np.arange(newElements.shape[0])+elementContainer.GetNumberOfElements())
        elementContainer.AddNewElements(newElements)

    return res

def ComputeFeatures(inputMesh: UnstructuredMesh, featureAngle:PBasicFloatType=90., skin:Optional[UnstructuredMesh]=None) -> Tuple[UnstructuredMesh,UnstructuredMesh]:
    """Compute features, element of dimensionality dim-2 (edges) for a given angle

    Parameters
    ----------
    inputMesh : UnstructuredMesh
        the input mesh
    featureAngle : PBasicFloatType, optional
        the detection angle: edges with faces with an angle grater than featureAngle will be consider as ridges   , by default 90.
    skin : Optional[UnstructuredMesh], optional
        if the user can provide the skin (to save time), by default None

    Returns
    -------
    Tuple[UnstructuredMesh,UnstructuredMesh]
       edgeMesh: a mesh containing the edges
       skinMesh: a mesh containing the skin
    """


    if skin is None:
        skinMesh = ComputeSkin(inputMesh)
        skinMeshSave = copy.deepcopy(skinMesh)
        # we have to merge all the 2D elements form the original mesh to the skin mesh
        for name,data,ids in ElementFilter(inputMesh, dimensionality = 2):
            skinMesh.elements[name].Merge(data)
        CleanDoubleElements(skinMesh)
    else:
        skinMesh = skin
        skinMeshSave = skinMesh

    md = skinMesh.GetPointsDimensionality()

    nex = skinMesh.GetNumberOfElements()

    #we use the original id to count the number of time the faces is used
    surf = {}
    for name,data in skinMesh.elements.items():
        if ElementNames.dimension[name] != md-1:
            continue
        faces = ElementNames.faces[name]
        numberOfNodes = ElementNames.numberOfNodes[name]

        ne = data.GetNumberOfElements()

        for faceType,localFaceConnectivity in faces:
            globalFaceConnectivity = data.connectivity[:,localFaceConnectivity]
            if not faceType in surf:
                surf[faceType] = {}
            surf2 = surf[faceType]
            key = np.array([nex**x for x in range(ElementNames.numberOfNodes[faceType]) ])
            for i in range(ne):
                bariCentre = np.sum(skinMesh.nodes[data.connectivity[i,:] ,:],axis=0)/numberOfNodes
                cc = globalFaceConnectivity[i,:]
                lc = np.sort(cc)

                elementHash = np.sum(lc*key)

                edgeVector = skinMesh.nodes[lc[0],:] - skinMesh.nodes[lc[1],:]
                planeVector = bariCentre - skinMesh.nodes[lc[1],:]
                normal = np.cross(edgeVector, planeVector)
                normal /= np.linalg.norm(normal)

                if elementHash in surf2:
                    surf2[elementHash][0] +=1
                    normal1 = surf2[elementHash][1]
                    cross = np.cross(normal, normal1)
                    norm = np.linalg.norm(cross)
                    norm =  norm if norm <=1 else 1.
                    angle = np.arcsin(norm)
                    surf2[elementHash][2] = 180*angle/np.pi
                else:
                    #[number of of used, normal of the first insertion,angle,   connectivity
                    surf2[elementHash] = [1,normal,None,cc]



    edgeMesh = UnstructuredMesh()
    edgeMesh.nodes = inputMesh.nodes
    edgeMesh.originalIDNodes = inputMesh.originalIDNodes
    numberOfNodes = inputMesh.GetNumberOfNodes()

    for name, surf2 in surf.items():
        if len(surf2) == 0:
            continue
        data = edgeMesh.GetElementsOfType(name)
        for data2 in surf2.values():
            if data2[0] == 1 or data2[0] > 2:
                data.AddNewElement(data2[3],-1)
            elif data2[0] == 2 and data2[2]  >= featureAngle:
                data.AddNewElement(data2[3],-1)

    for elType in [ElementNames.Bar_2, ElementNames.Bar_3]:
        if elType not in edgeMesh.elements:
            continue
        bars = edgeMesh.GetElementsOfType(elType)
        bars.tags.CreateTag("Ridges").SetIds(np.arange(bars.GetNumberOfElements()))

    skinMesh.PrepareForOutput()
    edgeMesh.PrepareForOutput()

    #Now we compute the corners

    #The corner are the points:
    #    1) where 3 Ridges meet
    #    2) or touched by only one Ridge
    #    3) or the angle of 2 ridges is bigger than the featureAngle


    extractedBars = edgeMesh.GetElementsOfType(ElementNames.Bar_2)
    extractedBars.tighten()

    originalBars = inputMesh.GetElementsOfType(ElementNames.Bar_2)
    originalBars.tighten()

    #corners
    mask = np.zeros((numberOfNodes), dtype=bool )

    almanac = {}
    for bars in [originalBars,extractedBars]:
        intMask = np.zeros((numberOfNodes), dtype=np.int8 )

        np.add.at(intMask, bars.connectivity.ravel(),1)

        mask[intMask > 2 ] = True
        mask[intMask == 1 ] = True

        for bar in range(bars.GetNumberOfElements()):
            id0,id1 =  bars.connectivity[bar,:]
            p0 = inputMesh.nodes[id0,:]
            p1 = inputMesh.nodes[id1,:]
            vec1 = (p1-p0)
            vec1 /= np.linalg.norm(vec1)
            for p,sign in zip([id0,id1],[1,-1]):
                if mask[p]:
                    #already in the mask no need to treat this point
                    continue

                if p in almanac:
                    vec2 = almanac[p][1]
                    norm = np.dot(vec1, vec2*sign)
                    norm = norm if norm <= 1 else 1
                    angle = (180/np.pi)*np.arccos( norm )
                    almanac[p][2] = angle
                    almanac[p][0] +=1
                else:
                    almanac[p] = [1,vec1*sign,id0]

    for p,data in almanac.items():
        if (180-data[2]) >= featureAngle:
            mask[p] = True


    edgeMesh.nodesTags.CreateTag("Corners").AddToTag(np.where(mask)[0])

    edgeMesh.PrepareForOutput()
    skinMesh.PrepareForOutput()
    return (edgeMesh,skinMeshSave)

def NodesPermutation(mesh: UnstructuredMesh, per:ArrayLike):
    """Function to do a permutation of the nodes in a mesh (in place)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh to be modified
    per : ArrayLike
        the permutation vector ([1,0,2,3] to permute first an second node)
    """

    newNodes = mesh.nodes[per,:]

    perII =np.argsort(per)

    mesh.nodes = newNodes
    for tag in mesh.nodesTags:
        ids = tag.GetIds()
        newIds = perII[ids]
        tag.SetIds(newIds)

    for data in mesh.elements.values():
        data.connectivity = perII[data.connectivity]


    for name,data in mesh.nodeFields.items():
        mesh.nodeFields[name] = mesh.nodeFields[name][per]

def AddTagPerBody(inmesh:UnstructuredMesh)-> np.ndarray:
    """Generate nodal tag (in the form of "Body_"+int) per body
     a body is defined by all the nodes connected by the elements (connectivity filter in vtk )

    Parameters
    ----------
    inmesh : UnstructuredMesh
        the input mesh

    Returns
    -------
    np.ndarray
        pointsPerBody : a vector with the number of point in every body.
        len(pointsPerBody) is the number of unconnected bodies in the mesh

    Raises
    ------
    Exception
        in the case of an internal error
    """
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ComputeNodeToNodeConnectivity
    dualGraph,usedPoints = ComputeNodeToNodeConnectivity(inmesh)

    # Connectivity walk
    nbOfNodes = inmesh.GetNumberOfNodes()
    treated = np.zeros(inmesh.GetNumberOfNodes(), dtype=bool)
    nextPoint = np.zeros(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType)

    # we start from the first point and body number 0
    nextPointCpt = 0
    bodyCpt = 0

    cpt = 0
    pointsPerBody = []
    while True:
        # we finish all the points
        if cpt == nbOfNodes :
            break

        # we already explored all the point in this body
        if cpt == nextPointCpt:
            initialPoint = np.argmax(treated == False)
            #(edge case) in the case we have only Trues in the treated
            if treated[initialPoint]:# pragma: no cover
                break
            treated[initialPoint] = True
            nextPoint[nextPointCpt] = initialPoint
            nextPointCpt +=1

            tagName = "Body_"+str(bodyCpt)
            tag = inmesh.GetNodalTag(tagName)
            tag.AddToTag(initialPoint)
            bodyCpt += 1
            pointsInThisBody = 1


        else:
            raise Exception ("Error in this function")# pragma: no cover


        while cpt < nextPointCpt:
            workingIndex = nextPoint[cpt]
            indexes = dualGraph[workingIndex,0:usedPoints[workingIndex]]

            for index in indexes:
                if not treated[index] :
                    treated[index] = True
                    nextPoint[nextPointCpt] = index
                    nextPointCpt +=1
                    tag.AddToTag(index)
                    pointsInThisBody +=1

            cpt += 1
        pointsPerBody.append(pointsInThisBody)
    return pointsPerBody

def Morphing(mesh: UnstructuredMesh, targetDisplacement, targetDisplacementMask, radius:Optional[PBasicFloatType]=None, forceOneStep=False)-> np.ndarray:
    """method for computing the deform mesh knowing displacement of some nodes.
    the user can push the morphed point back to the mesh by doing : mesh.node = morphedPoints

    note: https://www.researchgate.net/publication/288624175_Mesh_deformation_based_on_radial_basis_function_interpolation_computers_and_structure

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh, the use only the point position information
    targetDisplacement : _type_
        is the known displacement in a numpy array (shape [number_of_of_known_nodes,3])
    targetDisplacementMask : _type_
        array containing the ids of known nodes (list of ids or boolean array) in the same order as targetDisplacement
    rayon : Optional[PBasicFloatType], optional
        you can choose a radius by setting rayon to a value, by default np.linalg.norm(mesh.boundingMax-mesh.boundingMin)/2
    forceOneStep : bool, optional
        if True do the morphing in one step, by default False

    Returns
    -------
    np.ndarray
        morphedPoints: the morphed points
    """

    def dPhi(x):
        table=x>=1
        y=-20*x*(1-x)**3
        if len(np.argwhere(table))>1:
            y[table]=np.zeros(len(np.argwhere(table)))
        return y

    def Phi(x):
        table=x>=1
        y=(1-x)**4*(4*x+1)
        if len(np.argwhere(table))>1:
            y[table]=np.zeros(len(np.argwhere(table)))
        return y

    grad_max=0.8
    max_step=10
    nb_step=1
    if forceOneStep:
        nb_step=1
    step=0
    ##################################RBF###############################################

    nb_nodes = mesh.GetNumberOfNodes()
    new_nodes =  np.copy(mesh.GetPosOfNodes())

    if radius is None:
        mesh.ComputeBoundingBox()
        r=np.linalg.norm(mesh.boundingMax-mesh.boundingMin)/2
    else:
        r=radius
    if r==0:
        r=1
    while step<nb_step:
        border_nodes=new_nodes[targetDisplacementMask,:]

        rhs=targetDisplacement/nb_step
        M=np.eye(np.shape(border_nodes)[0])
        for j in range(np.shape(M)[0]):
            d=np.linalg.norm(border_nodes-border_nodes[j],axis=1)
            M[:,j]=Phi(d/r)
        op=M
        del M
        try:
            ab=np.linalg.solve(op,rhs)
        except np.linalg.LinAlgError:
            #print('Bad conditioning of RBF operator, using least square solution')
            ab=np.linalg.lstsq(op,rhs,cond=10**(9))
        del op
        alpha=ab
        ds=np.zeros((nb_nodes,mesh.GetDimensionality()))
        s=np.zeros((nb_nodes,mesh.GetDimensionality()))
        for j in range(np.shape(border_nodes)[0]):
            d=np.array([np.linalg.norm(new_nodes-border_nodes[j],axis=1)])
            s=s+(Phi(d/r)).T*alpha[j]
            ds=ds+(dPhi(d/r)).T*alpha[j]/r
        if step==0 and forceOneStep:
            nb_step=min(max_step,int(np.floor(np.max(np.linalg.norm(ds,axis=1)/grad_max)))+1)
            s=s/nb_step

        new_nodes += s
        step+=1

    return new_nodes

def LowerNodesDimension(mesh:UnstructuredMesh) -> UnstructuredMesh:
    """Remove the last columns of the mesh nodes coordinate.
    This is done in place

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh

    Returns
    -------
    UnstructuredMesh
        the mesh
    """
    newDim = mesh.GetPointsDimensionality() - 1
    mesh.nodes = np.ascontiguousarray(mesh.nodes[:,0:newDim])
    return mesh

def RigidBodyTransformation(mesh: UnstructuredMesh, rotationMatrix: np.ndarray, translationVector: np.ndarray):
    """In place rigid  body transformation.
    new pos = Q.pos + translation

    the rotation matrix Q should verify:
    QtQ = I = QQt et det Q = 1

    Parameters
    ----------
    mesh : UnstructuredMesh
        The input mesh
    rotationMatrix : np.ndArray
        a 3x3 rotation matrix
    translationVector : np.ndArray
        a vector of size  3 with the translation
    """

    assert np.linalg.norm(np.dot(rotationMatrix.T, rotationMatrix) - np.eye(3)) < 1.e-12
    assert np.linalg.norm(np.dot(rotationMatrix, rotationMatrix.T) - np.eye(3)) < 1.e-12
    assert (np.linalg.det(rotationMatrix) - 1) < 1.e-12

    mesh.nodes = np.dot(rotationMatrix, mesh.nodes.T).T
    for i in range(3):
        mesh.nodes[:,i] += translationVector[i]

def ComputeRigidBodyTransformationBetweenTwoSetOfPoints(setPoints1: np.ndarray, setPoints2: np.ndarray)-> Tuple[np.ndarray, np.ndarray]:
    """ Compute the rotation and the translation operator from two sets of points
    setPoints1 and setPoints2 have the same dimension (nbeOfPoints,dimension)

    Parameters
    ----------
    setPoints1 : np.ndarray
        First set of points
    setPoints2 : np.ndarray
        Second set of points

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A, b such that setPoints1.T approx A*setPoints2.T + b
    """

    assert setPoints1.shape == setPoints2.shape

    dim = setPoints1.shape[1]
    setPoints1 = np.hstack((setPoints1, np.ones((setPoints1.shape[0],1))))

    res = np.linalg.lstsq(setPoints1, setPoints2, rcond=None)[0]

    A = res[:dim,:dim].T
    b = res[dim,:].T

    return A, b

def ConvertNTagsToETags(mesh: UnstructuredMesh, prefix="NTag", targetDim:Optional[int]=None):
    """Create eTags from nTags on the mesh. Skin (surface element) and edges are created to hold the eTag
    information. The name of the created element tags is constructed using :

    prefix+str(d)+"D_"+nTagName

    Where d is the dimensionality of the target elements (1 for edges, 2 for surfaces)
    I.e.: "NTag1D_ForceEdgeNodes"

    Parameters
    ----------
    mesh : UnstructuredMesh
        The input mesh
    prefix : str, optional
        The prefix used for the creation of eTags, by default "NTag"
    targetDim : int | None, optional
        if None, surface and edge element with be created, by default None
        if 1 only edges with be created
        if 2 only surface elements are created
    """
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    elementFilter = ElementFilter(mesh,dimensionality=targetDim)

    nTags = mesh.nodesTags.keys()

    edges, skin = ComputeFeatures(mesh, featureAngle=80)

    toTreat = [mesh]
    if targetDim == 1 or targetDim == None:
        toTreat.append(edges)
    if targetDim == 2 or targetDim == None:
        toTreat.append(skin)

    for nTag in nTags:
        #mask = mesh.nodesTags[nTag].GetIdsAsMask(mesh.GetNumberOfNodes())

        for m in toTreat:
            elementFilter.mesh = m
            elementFilter.nTags = [nTag]
            for name, data, ids in elementFilter:
                d = ElementNames.dimension[name]
                name = prefix+str(d)+"D_"+nTag
                data.tags.CreateTag(name).SetIds(ids)

    toMerge =[]
    if targetDim == 1 or targetDim == None:
        toMerge.append(edges)
    if targetDim == 2 or targetDim == None:
        toMerge.append(skin)

    for m in toMerge:
        for name,data in m.elements.items():
            nbe  = data.GetNumberOfElements()
            if nbe ==0:
                continue
            mask = np.zeros(nbe, dtype=bool)
            for etag in data.tags.keys():
                np.logical_or(mask, data.tags[etag].GetIdsAsMask(nbe), out=mask)

            if not np.any(mask):
                continue

            subData = ExtractElementsByMask(data,mask)
            mesh.GetElementsOfType(name).Merge(subData)
    CleanDoubleElements(mesh)
#------------------------- CheckIntegrity ------------------------

def CheckIntegrity_AddTagPerBody(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
    from BasicTools.Containers.UnstructuredMeshCreationTools import MirrorMesh

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    resII = MirrorMesh(res,x=0,y=0,z=0)
    print( resII.nodes)
    print( resII.GetElementsOfType(ElementNames.Triangle_3))
    print(AddTagPerBody(resII))

    if len(resII.nodesTags) != 8:
        raise # pragma: no cover
    print(resII.nodesTags)

    return "ok"

def CheckIntegrity_CleanDoubleNodes(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh)

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points, tets, ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("OnePoint").SetIds([0,1])

    CleanDoubleNodes(mesh,tol =0)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points, tets, ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh, nodesToTestMask= np.array([True,False,True,False,True], dtype=bool) )
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover

    return "ok"

def CheckIntegrity_CleanLonelyNodes(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("OnePoint").SetIds([0,1])
    CleanLonelyNodes(mesh)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("OnePoint").SetIds([0,1])
    usedNodes, out  = CleanLonelyNodes(mesh,inPlace=False)

    if out.GetNumberOfNodes() != 4:
        print(out)
        raise# pragma: no cover
    return "ok"


def CheckIntegrity_ComputeFeatures(GUI =False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([2,3,4])
    myMesh.SetOrigin([-1.0,-1.0,-1.0])
    myMesh.SetSpacing([2., 2.,2]/myMesh.GetDimensions())
    print("ConstantRectilinearMesh")
    print(myMesh)
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex=True)
    print("UnstructuredMesh")
    print(res2)

    edges, skin = ComputeFeatures(res2, featureAngle=80)
    print("edges mesh")
    print(edges)
    print("SkinMesh")
    print(skin)
    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView

        OpenInParaView(mesh=edges,filename="edges.xmf")
        OpenInParaView(mesh=skin,filename="skin.xmf")

        for name,data in edges.elements.items():
           res2.GetElementsOfType(name).Merge(data)

        OpenInParaView(res2,filename="all+edges.xmf")
        print(res2)

    print(edges)
    return "ok"

def CheckIntegrity_ComputeSkin(GUI=False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([2,3,4])
    myMesh.SetOrigin([-1.0,-1.0,-1.0])
    myMesh.SetSpacing([2., 2.,2]/myMesh.GetDimensions())
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex=True)

    skin = ComputeSkin(res2,inPlace=False)
    print(skin)
    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(skin,filename="CheckIntegrity_ComputeSkin_InPlace_False.xmf")


    print(res2)
    ComputeSkin(res2,inPlace=True)
    print(res2)
    if GUI :
        OpenInParaView(res2,filename="CheckIntegrity_ComputeSkin_InPlace_True.xmf")


    return "ok"

def CheckIntegrity_Morphing(GUI = False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)
    targetDisplacement = np.empty((mesh.GetNumberOfNodes(),3),dtype=float)
    targetDisplacementMask = np.empty(mesh.GetNumberOfNodes(),dtype=int)
    cpt = 0
    print(mesh)
    for name,data in mesh.elements.items():

        if ElementNames.dimension[name] != 2:
            continue

        ids = data.GetNodesIdFor(data.GetTag("X0").GetIds())
        print(ids)
        targetDisplacementMask[cpt:cpt+len(ids)] = ids
        targetDisplacement[cpt:cpt+len(ids),:] = 0
        cpt += len(ids)

        ids = data.GetNodesIdFor(data.GetTag("X1").GetIds())
        print(ids)
        targetDisplacementMask[cpt:cpt+len(ids)] = ids
        targetDisplacement[cpt:cpt+len(ids),:] = [[0,0,10]]
        cpt += len(ids)

    targetDisplacement = targetDisplacement[0:cpt,:]
    targetDisplacementMask = targetDisplacementMask[0:cpt]

    new_p1 = Morphing(mesh, targetDisplacement,targetDisplacementMask)
    new_p2 = Morphing(mesh, targetDisplacement,targetDisplacementMask, radius= 20.)

    new_p0 = np.copy(mesh.nodes)
    new_p0[targetDisplacementMask,:] += targetDisplacement
    mesh.nodeFields["morph0"] = new_p0
    mesh.nodeFields["morph1"] = new_p1
    mesh.nodeFields["morph2"] = new_p2


    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

    return "ok"

def CheckIntegrity_NodesPermutation(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
    points = [[1,2,3],[4,5,6],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],[3,1,0,2]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    NodesPermutation(mesh, [1,0,2,3])

    print( mesh.nodes[0,:] )
    if np.any( mesh.nodes[0,:]-[4,5,6] ):
        raise(Exception("error in the permutation"))

    return "ok"

def CheckIntegrity_DeleteInternalFaces(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    ## add 2 triangles
    triangles = mesh.GetElementsOfType(ElementNames.Triangle_3)
    triangles.AddNewElement([0,1,2],0)
    triangles.AddNewElement([3,0,2],1)
    triangles.AddNewElement([3,0,2],2)
    #print(mesh)

    DeleteInternalFaces(mesh)
    print(mesh)
    return "ok"

def CheckIntegrity_RigidBodyTransformation(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    ## add 2 triangles
    triangles = mesh.GetElementsOfType(ElementNames.Triangle_3)
    triangles.AddNewElement([0,1,2],0)
    triangles.AddNewElement([3,0,2],1)
    triangles.AddNewElement([3,0,2],2)
    #print(mesh)

    theta = np.pi/2.
    A = np.array([[1, 0, 0],
                  [0, np.cos(theta), -np.sin(theta)],
                  [0, np.sin(theta), np.cos(theta)]])
    b = np.array([1,0,0])
    RigidBodyTransformation(mesh, A, b)
    return "ok"




def CheckIntegrity_ComputeRigidBodyTransformationBetweenTwoSetOfPoints(GUI=False):

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points1 = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4]]
    mesh1 = CreateMeshOf(points1,tets,ElementNames.Tetrahedron_4)
    mesh2 = CreateMeshOf(points1,tets,ElementNames.Tetrahedron_4)

    theta = np.pi/3.
    A = np.array([[1, 0, 0],
                  [0, np.cos(theta), -np.sin(theta)],
                  [0, np.sin(theta), np.cos(theta)]])
    b = np.array([1,0,0])

    RigidBodyTransformation(mesh2, A, b)

    A, b = ComputeRigidBodyTransformationBetweenTwoSetOfPoints(mesh1.nodes, mesh1.nodes)
    return "ok"

def CheckIntegrity_DeleteElements(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
    points1 = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2],[1,2,3],[2,3,4]]
    mesh1 = CreateMeshOf(points1,tets,ElementNames.Triangle_3)
    mask = np.zeros(3)
    mask[1,2]
    DeleteElements(mesh1, mask )

def CheckIntegrity_CopyElementTags(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
    points1 = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2],[1,2,3],[2,3,4]]
    mesh1 = CreateMeshOf(points1,tets,ElementNames.Triangle_3)
    mesh1.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("fromMesh1").SetIds([0,2])
    mesh1.GetElementsOfType(ElementNames.Bar_2).AddNewElement([0,1],0)
    mesh1.GetElementsOfType(ElementNames.Bar_3)

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
    points1 = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,3],[2,3,4],[1,2,3],]
    mesh2 = CreateMeshOf(points1,tets,ElementNames.Triangle_3)
    mesh2.GetElementsOfType(ElementNames.Bar_2)

    CopyElementTags(UnstructuredMesh(),mesh2)
    CopyElementTags(mesh1,UnstructuredMesh())
    CopyElementTags(mesh1,mesh2)
    CopyElementTags(mesh1,mesh2,extendTags=True)

    ids2 = mesh2.GetElementsOfType(ElementNames.Triangle_3).tags["fromMesh1"].GetIds()
    print("ids2  -> ", ids2)
    if len(ids2) != 1 or ids2[0] != 1:
        print(mesh1)
        print(mesh2)
        raise Exception("KO")
    return "ok"

def CheckIntegrity_CleanDoubleElements(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points1 = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4],[3,0,2,4], [1,0,2,4] ]
    mesh1 = CreateMeshOf(points1,tets,ElementNames.Tetrahedron_4)
    mesh1.GetElementsOfType(ElementNames.Tetrahedron_4).originalIds = np.array([5,6,7,8])
    mesh1.GetElementsOfType(ElementNames.Triangle_3)
    mesh1.elemFields["testField"] = np.arange(4)
    mesh1.elemFields["testField4x3"] = np.zeros((4,3))

    CleanDoubleElements(mesh1)
    if mesh1.GetNumberOfElements() != 3:
        raise Exception("ko")
    print(mesh1.GetElementsOfType(ElementNames.Tetrahedron_4).originalIds)

    return "ok"

def CheckIntegrity(GUI=False):
    functionsToTest= [
        CheckIntegrity_NodesPermutation,
        CheckIntegrity_Morphing,
        CheckIntegrity_ComputeFeatures,
        CheckIntegrity_ComputeSkin,
        CheckIntegrity_DeleteInternalFaces,
        CheckIntegrity_CleanDoubleNodes,
        CheckIntegrity_CleanLonelyNodes,
        CheckIntegrity_AddTagPerBody,
        CheckIntegrity_RigidBodyTransformation,
        CheckIntegrity_ComputeRigidBodyTransformationBetweenTwoSetOfPoints,
        CheckIntegrity_CleanDoubleElements,
        CheckIntegrity_CopyElementTags
    ]
    for f in functionsToTest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
