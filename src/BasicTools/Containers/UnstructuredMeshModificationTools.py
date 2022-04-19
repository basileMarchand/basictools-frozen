# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import copy

import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.MeshBase import Tags
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.NumpyDefs import PBasicIndexType

def CleanDoubleNodes(res, tol = None, nodesToTestMask= None):

    res.ComputeBoundingBox()
    if tol is None:
        tol = np.linalg.norm(res.boundingMax - res.boundingMin)*1e-7

    nbnodes = res.GetNumberOfNodes()
    toKeep = np.zeros(nbnodes, dtype=bool )
    newindex = np.zeros(nbnodes, dtype=PBasicIndexType)
    tol2 = tol**2
    def dist2(array,value):
        d = array-value
        return np.inner(d,d)

    if nodesToTestMask is None and tol == 0:
        # optimized version for tol = 0.0
        database = {}
        cpt  = 0
        for i in range(nbnodes):
            point = tuple(res.nodes[i,:])
            ni = database.get(point)
            if ni is None:
                database[point] = cpt
                newindex[i] = cpt
                cpt +=1
                toKeep[i] = True
            else:
                newindex[i] = ni

    elif nodesToTestMask is None :

        from BasicTools.Containers.Octree import Octree
        cpt  = 0

        Ma = res.boundingMax
        Mi = res.boundingMin
        diag2 = np.linalg.norm(res.boundingMax-res.boundingMin)**2
        tree = Octree(Ma[0],Ma[1],Ma[2], Mi[0], Mi[1], Mi[2])

        for i in range(nbnodes):
            point = tuple(res.nodes[i,:])
            entries = tree.find_within_range_cube(point, tol)

            odist = diag2
            if len(entries):
                for entry in entries:
                    dist = (point[0] - entry[0][0])**2   +(point[1] - entry[0][1])**2+(point[2] - entry[0][2])**2
                    if dist < odist:
                        odist = dist
                        index = entry[1]

            if odist <= tol2:
                newindex[i] = index
            else:
                tree.add_item(cpt, point )
                newindex[i] = cpt
                cpt += 1
                toKeep[i] = True


    else:
        cpt =0
        for i in range(nbnodes):
            if not nodesToTestMask[i]:
                newindex[i] = cpt
                cpt += 1
                toKeep[i] = True
                continue

            posi = res.nodes[i,:]

            for j in range(i):
                if not nodesToTestMask[j]:
                    continue
                if toKeep[j]:
                    dist = dist2(posi,res.nodes[j,:] )
                    if dist < tol2 :
                        newindex[i] = newindex[j]
                        break
            else:
                newindex[i] = cpt
                cpt += 1
                toKeep[i] = True

    res.nodes = res.nodes[toKeep,:]
    res.originalIDNodes = np.where(toKeep)[0]

    for tag in res.nodesTags :
        tag.SetIds(np.unique(newindex[tag.GetIds()]) )

    for elementName in res.elements:
        elements = res.elements[elementName]
        elements.connectivity = newindex[elements.connectivity]

def CleanLonelyNodes(res,out=None):

    usedNodes = np.zeros(res.GetNumberOfNodes(),dtype=bool )
    for data in res.elements.values():
        usedNodes[data.connectivity.ravel()] = True

    cpt = 0
    NewIndex =  np.zeros(res.GetNumberOfNodes(),dtype=PBasicIndexType )-1
    for n in range(res.GetNumberOfNodes()):
        if usedNodes[n]:
            NewIndex[n] = cpt
            cpt += 1
    originalIDNodes = np.where(usedNodes)[0]
    #filter the nodes
    #inplace
    if out is None:
        res.nodes = res.nodes[usedNodes ,:]
        res.originalIDNodes = originalIDNodes
        newTags = Tags()
        #node tags
        for tag in res.nodesTags :
            newTags.CreateTag(tag.name).SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])
        res.nodesTags = newTags

        #renumbering the connectivity matrix
        for elementName in res.elements:
            elements = res.elements[elementName]
            elements.connectivity = NewIndex[elements.connectivity]

    else:
        out.nodes = res.nodes[usedNodes ,:]
        out.originalIDNodes = originalIDNodes
        #node tags
        for tag in res.nodesTags :
            outtag = out.nodesTags.CreateTag(tag.name)
            outtag.SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])

        for elementName in res.elements:
            elements = res.elements[elementName]
            outelements = out.GetElementsOfType(elementName)
            outelements.connectivity = NewIndex[elements.connectivity]

            outelements.tags = copy.deepcopy(elements.tags)

    for name,data in res.nodeFields.items():
        res.nodeFields[name] = res.nodeFields[name][usedNodes]

    return usedNodes

def CleanEmptyTags(inmesh):
    for data in inmesh.elements.values():
        data.tags.RemoveEmptyTags()
    inmesh.nodesTags.RemoveEmptyTags()

def CleanDoubleElements(mesh, preserveOriginalIds=False):
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    elemFieldMask = np.zeros((mesh.GetNumberOfElements(),),dtype=int)
    mcpt = 0 # mask counter
    ecpt = 0 # element counter
    newElements = type(mesh.elements)()
    for elemName, data in mesh.elements.items():

        nbe = data.GetNumberOfElements()
        if  nbe == 0:
            continue
        data.tighten()
        _, index, = np.unique(np.sort(data.connectivity,axis=1),return_index=True,axis=0)
        mask = np.zeros(nbe,dtype=bool)
        mask[index] = True
        newElems = ExtractElementsByMask(data, mask)
        newElements[elemName] = newElems
        elemFieldMask[mcpt:mcpt+newElems.GetNumberOfElements()] = newElems.originalIds+ecpt
        mcpt  += newElems.GetNumberOfElements()

        if preserveOriginalIds:
            newElements[elemName].originalIds = data.originalIds[newElements[elemName].originalIds]

    mesh.elements = newElements
    elemFieldMask = elemFieldMask[0:mcpt]
    for name, data  in  mesh.elemFields.items():
        if len(mesh.elemFields[name].shape) == 1:
            mesh.elemFields[name] = mesh.elemFields[name][elemFieldMask]
        else:
            mesh.elemFields[name] = mesh.elemFields[name][elemFieldMask,:]

def CopyElementTags(mesh1,mesh2,extendTags=False):
    """Copy tags from one mesh1 to mesh2. We use the connectivity to identify the elements
    extendTags == False will overwrite the thags on mesh2
    extendTags == True will extend the thags on mesh2
    """
    for elemName, data1 in mesh1.elements.items():
        data1.tighten()

        nbe1 = data1.GetNumberOfElements()
        if nbe1 == 0 :
            continue
        if elemName not in mesh2.elements :
            continue
        data2 = mesh2.elements[elemName]
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

def DeleteElements(mesh,mask):
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    OriginalNumberOfElements = mesh.GetNumberOfElements()
    datatochange = {}

    mesh.ComputeGlobalOffset()

    for name,data in mesh.elements.items():
        offset = data.globaloffset
        localmask = mask[offset:offset+data.GetNumberOfElements()]
        datatochange[name] = ExtractElementsByMask(data, np.logical_not(localmask))

    for name, data in datatochange.items():
        mesh.elements[name] = data

    if OriginalNumberOfElements != mesh.GetNumberOfElements():
        print("Number Of Elements Changed: Droping elemFields")
        mesh.elemFields = {}

    mesh.PrepareForOutput()

def DeleteInternalFaces(mesh):
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    OriginalNumberOfElements = mesh.GetNumberOfElements()
    skin = ComputeSkin(mesh)
    datatochange = {}

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
            ehash = np.sum(lc*key)
            surf2[ehash] = [1,cc]

        for i in range(ne):
            cc = data.connectivity[i,:]
            lc = np.sort(cc)

            ehash = np.sum(lc*key)
            if ehash in surf2:
                mask[i] = True

        if np.any(mask):
            datatochange[name] = ExtractElementsByMask(data,mask)

    for name, data in datatochange.items():
        mesh.elements[name] = data

    if OriginalNumberOfElements != mesh.GetNumberOfElements():
        print("Number Of Elements Changed: Droping elemFields")
        mesh.elemFields = {}

    mesh.PrepareForOutput()
    return mesh

def ComputeSkin(mesh, md=None ,inplace=False):
    """ Compute the skin of a mesh (mesh), if md (mesh dimensionality) (md) is None
    the mesh.getdimensionality() is used to filter the element to compute
    the skin.

    if inplace=True the skin is added to the original mesh, a tag "skin" is created.

    if inplace==False a new mesh is returned with the element of the skin (even if
    some element of the skin are already present on the original mesh).
    If the user merged the two meshes, you must call CleanDoubleElements to clean
    the mesh

    Warning if some elements are duplicated the behaviour of computeSkin with
    the option inplace==True is not defined (Use CleanDoubleElements before
    ComputeSkin)
    """

    if md is None:
        md = mesh.GetDimensionality()

    # first we count the total number of potential individual skin elemnts
    totalcpt = {}
    Dfilter = ElementFilter(mesh,dimensionality = md)
    for elemName, data, ids in Dfilter:
        faces = ElementNames.faces[elemName]
        nbelements = data.GetNumberOfElements()
        for faceType,localFaceConnectivity in faces:
            cpt = totalcpt.setdefault(faceType,0) + nbelements
            totalcpt[faceType] =  cpt

    if inplace :
        # we add the element already present on the mesh
        D1filter = ElementFilter(mesh,dimensionality = md-1)
        for elemName, data ,ids  in D1filter:
            nbelements = data.GetNumberOfElements()
            cpt = totalcpt.setdefault(elemName,0) + nbelements
            totalcpt[elemName] =  cpt

    #  Allocation of matrices to store all skin elements
    #sortedstorage = {}
    storage = {}
    for k,v in totalcpt.items():
        nn = ElementNames.numberOfNodes[k]
        storage[k] = np.empty((v,nn),dtype=int)
        totalcpt[k] = 0

    if inplace :
        # fill storage with skin element aready on the mesh
        for elemName, data, ids in D1filter:
            cpt = totalcpt[elemName]
            nbelements = data.GetNumberOfElements()
            storage[elemName][cpt:cpt+nbelements,:] = data.connectivity
            totalcpt[elemName] += nbelements

    # fill storage with skin of elements of dimensionality D
    for elemName,data,ids in Dfilter:
        faces = ElementNames.faces[elemName]
        nbelements = data.GetNumberOfElements()
        for faceType,localFaceConnectivity in faces:
            globalFaceConnectivity = data.connectivity[:,localFaceConnectivity]
            cpt = totalcpt[faceType]
            storage[faceType][cpt:cpt+nbelements,:] = globalFaceConnectivity
            totalcpt[faceType] += nbelements

    if inplace:
        res = mesh
    else:
        res = UnstructuredMesh()
        res.nodesTags =  mesh.nodesTags
        res.originalIDNodes =  mesh.originalIDNodes
        res.nodes = mesh.nodes
        res.elements = type(mesh.elements)()


    # recover only the unique elements
    for elemName,cpt in totalcpt.items():
        if cpt == 0:
            continue
        store = storage[elemName]
        _,index,counts = np.unique(np.sort(store,axis=1),return_index=True,return_counts=True,axis=0)

        # the elements
        if inplace:

            uniqueelems = index[counts==2]

            if elemName in mesh.elements:
                melements = mesh.elements.GetElementsOfType(elemName)
                nbmelems = melements.GetNumberOfElements()
                ids =  uniqueelems[uniqueelems < nbmelems]
                # tag element already present on the original mesh
                melements.tags.CreateTag("Skin",False).SetIds(ids)
            else:
                nbmelems = 0
        else:
            nbmelems = 0

        # all the elements present only 1 time not in the original mesh
        uniqueelems = index[counts==1]
        ids =  uniqueelems[uniqueelems >= nbmelems]
        if len(ids) == 0 :
            continue

        if inplace:
            elems = mesh.elements.GetElementsOfType(elemName)
        else:
            elems = res.GetElementsOfType(elemName)


        newelements = store[ids,:]
        # add and tag new elements into the output
        elems.tags.CreateTag("Skin",False).AddToTag(np.arange(newelements.shape[0])+elems.GetNumberOfElements())
        elems.AddNewElements(newelements)

    return res

def ComputeFeatures(inputmesh, featureAngle=90, skin=None):
    if skin is None:
        skinmesh = ComputeSkin(inputmesh)
        skinmeshSave = copy.deepcopy(skinmesh)
        # we have to merge all the 2D elements form the original mesh to the skinmesh
        for name,data,ids in ElementFilter(inputmesh, dimensionality = 2):
            skinmesh.elements[name].Merge(data)
        CleanDoubleElements(skinmesh)
    else:
        skinmesh = skin
        skinmeshSave = skinmesh

    md = skinmesh.GetDimensionality()

    nex = skinmesh.GetNumberOfElements()

    #we use the original id to count the number of time the faces is used
    surf = {}
    for name,data in skinmesh.elements.items():
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
                baricentre = np.sum(skinmesh.nodes[data.connectivity[i,:] ,:],axis=0)/numberOfNodes
                cc = globalFaceConnectivity[i,:]
                lc = np.sort(cc)

                ehash = np.sum(lc*key)

                edgeVector = skinmesh.nodes[lc[0],:] - skinmesh.nodes[lc[1],:]
                planeVector = baricentre - skinmesh.nodes[lc[1],:]
                normal = np.cross(edgeVector, planeVector)
                normal /= np.linalg.norm(normal)

                if ehash in surf2:
                    surf2[ehash][0] +=1
                    normal1 = surf2[ehash][1]
                    cross = np.cross(normal, normal1)
                    angle = np.arcsin(np.linalg.norm(cross))
                    surf2[ehash][2] = 180*angle/np.pi
                else:
                    #[number of of used, normal of the first insertion,angle,   connectivity
                    surf2[ehash] = [1,normal,None,cc]



    edgemesh = UnstructuredMesh()
    edgemesh.nodes = inputmesh.nodes
    edgemesh.originalIDNodes = inputmesh.originalIDNodes
    numberOfNodes = inputmesh.GetNumberOfNodes()

    for name, surf2 in surf.items():
        if len(surf2) == 0:
            continue
        data = edgemesh.GetElementsOfType(name)
        for data2 in surf2.values():
            if data2[0] == 1 or data2[0] > 2:
                data.AddNewElement(data2[3],-1)
            elif data2[0] == 2 and data2[2]  >= featureAngle:
                data.AddNewElement(data2[3],-1)

    for eltype in [ElementNames.Bar_2, ElementNames.Bar_3]:
        if eltype not in edgemesh.elements:
            continue
        bars = edgemesh.GetElementsOfType(eltype)
        bars.tags.CreateTag("Ridges").SetIds(np.arange(bars.GetNumberOfElements()))
    skinmesh.PrepareForOutput()
    edgemesh.PrepareForOutput()

    """
    Now we compute the corners

    The corner are the points:
        1) where 3 Ridges meet
        2) or touched by only one Ridge
        3) or the angle of 2 ridges is bigger than the featureAngle
    """

    extractedBars = edgemesh.GetElementsOfType(ElementNames.Bar_2)
    extractedBars.tighten()

    originalBars = inputmesh.GetElementsOfType(ElementNames.Bar_2)
    originalBars.tighten()

    #corners
    mask = np.zeros((numberOfNodes), dtype=bool )

    almanac = {}
    for bars in [originalBars,extractedBars]:
        intmask = np.zeros((numberOfNodes), dtype=int )

        np.add.at(intmask, bars.connectivity.ravel(),1)

        mask[intmask > 2 ] = True
        mask[intmask == 1 ] = True

        for bar in range(bars.GetNumberOfElements()):
            id0,id1 =  bars.connectivity[bar,:]
            p0 = inputmesh.nodes[id0,:]
            p1 = inputmesh.nodes[id1,:]
            vec1 = (p1-p0)
            vec1 /= np.linalg.norm(vec1)
            for p,sign in zip([id0,id1],[1,-1]):
                if mask[p]:
                    #already in the mask no need to treat this point
                    continue

                if p in almanac:
                    vec2 = almanac[p][1]
                    angle = (180/np.pi)*np.arccos(np.dot(vec1, vec2*sign) )
                    almanac[p][2] = angle
                    almanac[p][0] +=1
                else:
                    almanac[p] = [1,vec1*sign,id0]

    for p,data in almanac.items():
        if (180-data[2]) >= featureAngle:
            mask[p] = True


    edgemesh.nodesTags.CreateTag("Corners").AddToTag(np.where(mask)[0])

    edgemesh.PrepareForOutput()
    skinmesh.PrepareForOutput()
    return (edgemesh,skinmeshSave)

def NodesPermutation(mesh,per):
    """
    Function to do a permutation of the nodes in a mesh (inplace)

    mesh : (UnstructuredMesh) mesh to be modified
    per : the permutation vector ([1,0,2,3] to permute first an secod node)
    """
    nnodes = mesh.nodes[per,:]

    perII =np.argsort(per)

    mesh.nodes = nnodes
    for tag in mesh.nodesTags:
        ids = tag.GetIds()
        nids = perII[ids]
        tag.SetIds(nids)

    for data in mesh.elements.values():
        data.connectivity = perII[data.connectivity]

# to generate one tag per body
# a body is defines by all the nodes connected by the elements
def AddTagPerBody(inmesh):
    from BasicTools.Containers.UnstructuredMeshInspectionTools import GetDualGraph
    dualGraph,usedPoints = GetDualGraph(inmesh)

    # Connectivity walk
    nbOfNodes = inmesh.GetNumberOfNodes()
    treated = np.zeros(inmesh.GetNumberOfNodes(),dtype=bool)
    nextpoint = np.zeros(inmesh.GetNumberOfNodes(),dtype=np.int_)

    # we start from the first point and body number 0
    nextpointcpt = 0
    bodyCpt = 0

    cpt = 0
    pointsPerBody = []
    while True:
        # we finish all the points
        if cpt == nbOfNodes :
            break

        # we already explored all the point in this body
        if cpt == nextpointcpt:
            initialPoint = np.argmax(treated == False)
            #(edge case) in the case we have only Trues in the treated
            if treated[initialPoint]:# pragma: no cover
                break
            treated[initialPoint] = True
            nextpoint[nextpointcpt] = initialPoint
            nextpointcpt +=1

            tagName = "Body_"+str(bodyCpt)
            tag = inmesh.GetNodalTag(tagName)
            tag.AddToTag(initialPoint)
            bodyCpt += 1
            pointsInThisBody = 1


        else:
            raise Exception ("Error in this function")# pragma: no cover


        while cpt < nextpointcpt:
            workingIndex = nextpoint[cpt]
            indexes = dualGraph[workingIndex,0:usedPoints[workingIndex]]

            for index in indexes:
                if not treated[index] :
                    treated[index] = True
                    nextpoint[nextpointcpt] = index
                    nextpointcpt +=1
                    tag.AddToTag(index)
                    pointsInThisBody +=1

            cpt += 1
        pointsPerBody.append(pointsInThisBody)
    return pointsPerBody

def Morphing(mesh,BC,bord_tot,rayon=None,tridim=False,IDW=False):
    """
    method for computing the deform mesh knowing displacement of some nodes
    BC is the known displacement in a numpy array (shape [number_of_of_known_nodes,3])
    bord_tot contains the ids of known nodes (list of ids or boolean array)  in same order as BC
    you can choose a radius by setting rayon to a value
    put tridim to True if you want to do the morphing in one step
    put IDW to True if you want to use another morphing method without system inversion

    https://www.researchgate.net/publication/288624175_Mesh_deformation_based_on_radial_basis_function_interpolation_computers_and_structure
    """
    def dphi(x):
        table=x>=1
        y=-20*x*(1-x)**3
        if len(np.argwhere(table))>1:
            y[table]=np.zeros(len(np.argwhere(table)))
        return y

    def phi(x):
        table=x>=1
        y=(1-x)**4*(4*x+1)
        if len(np.argwhere(table))>1:
            y[table]=np.zeros(len(np.argwhere(table)))
        return y

    grad_max=0.8
    max_step=10
    nb_step=1
    if tridim:
        nb_step=1
    step=0
    ##################################RBF###############################################

    nb_nodes = mesh.GetNumberOfNodes()
    new_nodes =  np.copy(mesh.GetPosOfNodes())

    if rayon is None:
        mesh.ComputeBoundingBox()
        r=np.linalg.norm(mesh.boundingMax-mesh.boundingMin)/2
    else:
        r=rayon
    if r==0:
        r=1
    while step<nb_step:
        border_nodes=new_nodes[bord_tot,:]

        rhs=BC/nb_step
        M=np.eye(np.shape(border_nodes)[0])
        #print('Building RBF operator {}/{}(shape {})'.format(step+1,nb_step,np.shape(M)))
        for j in range(np.shape(M)[0]):
            d=np.linalg.norm(border_nodes-border_nodes[j],axis=1)
            M[:,j]=phi(d/r)
        op=M
        del M
        #print('Solving RBF problem (shape {})'.format(np.shape(op)))
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
            s=s+(phi(d/r)).T*alpha[j]
            ds=ds+(dphi(d/r)).T*alpha[j]/r
        if step==0 and tridim:
            nb_step=min(max_step,int(np.floor(np.max(np.linalg.norm(ds,axis=1)/grad_max)))+1)
            s=s/nb_step

        new_nodes += s
        step+=1

    return new_nodes

def LowerNodesDimension(mesh):
    newDim = mesh.GetDimensionality() - 1
    newNodes = np.empty((mesh.GetNumberOfNodes(),newDim))
    for i in range(mesh.GetNumberOfNodes()):
        newNodes[i,:] = mesh.nodes[i,0:newDim]
    mesh.nodes = newNodes
    return mesh


def RigidBodyTransformation(mesh, rotationMatrix, translationVector):
    """
    the rotation matrix Q should verify:
    QtQ = I = QQt et det Q = 1
    """
    assert np.linalg.norm(np.dot(rotationMatrix.T, rotationMatrix) - np.eye(3)) < 1.e-12
    assert np.linalg.norm(np.dot(rotationMatrix, rotationMatrix.T) - np.eye(3)) < 1.e-12
    assert (np.linalg.det(rotationMatrix) - 1) < 1.e-12

    mesh.nodes = np.dot(rotationMatrix, mesh.nodes.T).T
    for i in range(3):
        mesh.nodes[:,i] += translationVector[i]


def ComputeRigidBodyTransformationBetweenTwoSetOfPoints(setPoints1, setPoints2):
    """
    setPoints1 and setPoints2 have the same dimension (nbeOfPoints,dimension)

    returns A, b such that setPoints1.T approx A*setPoints2.T + b
    """
    assert setPoints1.shape == setPoints2.shape

    dim = setPoints1.shape[1]
    setPoints1 = np.hstack((setPoints1, np.ones((setPoints1.shape[0],1))))

    res = np.linalg.lstsq(setPoints1, setPoints2, rcond=None)[0]

    A = res[:dim,:dim].T
    b = res[dim,:].T

    return A, b

def ConvertNTagsToETags(mesh,prefix="NTag",targetDim=None):
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByMask

    filt = ElementFilter(mesh,dimensionality=targetDim)

    ntags = mesh.nodesTags.keys()

    edges, skin = ComputeFeatures(mesh, featureAngle=80)

    for ntag in ntags:
        mask = mesh.nodesTags[ntag].GetIdsAsMask(mesh.GetNumberOfNodes())

        for m in [edges, skin, mesh]:
            filt.mesh = m
            for name,data,ids in filt:
                a = np.sum(mask[data.connectivity],axis=1)
                elids = np.where(a==data.GetNumberOfNodesPerElement())[0]
                d = ElementNames.dimension[name]
                name = prefix+str(d)+"D_"+ntag
                data.tags.CreateTag(name).SetIds(elids)

    tomerge =[]
    if targetDim == 1 or targetDim == None:
        tomerge.append(edges)
    if targetDim == 2 or targetDim == None:
        tomerge.append(skin)

    for m in tomerge:
        for name,data in m.elements.items():
            ne  = data.GetNumberOfElements()
            mask = np.zeros(ne,dtype=bool)
            for etag in data.tags.keys():
                np.logical_or(mask, data.tags[etag].GetIdsAsMask(ne), out=mask)

            if not np.any(mask):
                continue

            if data.GetNumberOfElements()==0:
                continue

            subdata = ExtractElementsByMask(data,mask)
            mesh.GetElementsOfType(name).Merge(subdata)

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

def CheckIntegrity_CleanEmptyTags(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetNodalTag("Point0").AddToTag(0)
    res.GetNodalTag("EmptyTagN")

    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("EmptyTagE")
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("Element0").AddToTag(0)

    CleanEmptyTags(res)
    if len(res.nodesTags) != 1 :
        raise# pragma: no cover

    if "EmptyTagE" in  res.GetNamesOfElemTags() :
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_CleanDoubleNodes(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh)

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("OnePoint").SetIds([0,1])

    CleanDoubleNodes(mesh,tol =0)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh, nodesToTestMask= np.array([True,False,True,False,True],dtype=bool) )
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

    out = UnstructuredMesh()
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("OnePoint").SetIds([0,1])
    CleanLonelyNodes(mesh,out=out)

    if out.GetNumberOfNodes() != 4:
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
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=True)
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
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=True)

    skin = ComputeSkin(res2,inplace=False)
    print(skin)
    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(skin,filename="CheckIntegrity_ComputeSkin_Inplace_False.xmf")


    print(res2)
    ComputeSkin(res2,inplace=True)
    print(res2)
    if GUI :
        OpenInParaView(res2,filename="CheckIntegrity_ComputeSkin_Inplace_True.xmf")


    return "ok"

def CheckIntegrity_Morphing(GUI = False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)
    BC = np.empty((mesh.GetNumberOfNodes(),3),dtype=float)
    bord_tot = np.empty(mesh.GetNumberOfNodes(),dtype=int)
    cpt = 0
    print(mesh)
    for name,data in mesh.elements.items():

        if ElementNames.dimension[name] != 2:
            continue

        ids = data.GetNodesIdFor(data.GetTag("X0").GetIds())
        print(ids)
        bord_tot[cpt:cpt+len(ids)] = ids
        BC[cpt:cpt+len(ids),:] = 0
        cpt += len(ids)

        ids = data.GetNodesIdFor(data.GetTag("X1").GetIds())
        print(ids)
        bord_tot[cpt:cpt+len(ids)] = ids
        BC[cpt:cpt+len(ids),:] = [[0,0,10]]
        cpt += len(ids)

    BC = BC[0:cpt,:]
    bord_tot = bord_tot[0:cpt]

    new_p1 = Morphing(mesh, BC,bord_tot)
    new_p2 = Morphing(mesh, BC,bord_tot,rayon= 20. )

    new_p0 = np.copy(mesh.nodes)
    new_p0[bord_tot,:] += BC
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
        raise(Exception("error in the pertmutation"))

    return "ok"

def CheckIntegrity_DeleteInternalFaces(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    ## add 2 tris
    tris = mesh.GetElementsOfType(ElementNames.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([3,0,2],1)
    tris.AddNewElement([3,0,2],2)
    #print(mesh)

    DeleteInternalFaces(mesh)
    print(mesh)
    return "ok"

def CheckIntegrity_RigidBodyTransformation(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    ## add 2 tris
    tris = mesh.GetElementsOfType(ElementNames.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([3,0,2],1)
    tris.AddNewElement([3,0,2],2)
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
    mesh1.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("frommesh1").SetIds([0,2])
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

    ids2 = mesh2.GetElementsOfType(ElementNames.Triangle_3).tags["frommesh1"].GetIds()
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

    CleanDoubleElements(mesh1, preserveOriginalIds=False)
    if mesh1.GetNumberOfElements() != 3:
        raise Exception("ko")
    print(mesh1.GetElementsOfType(ElementNames.Tetrahedron_4).originalIds)


    points1 = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4],[3,0,2,4], [1,0,2,4] ]
    mesh1 = CreateMeshOf(points1,tets,ElementNames.Tetrahedron_4)
    mesh1.GetElementsOfType(ElementNames.Tetrahedron_4).originalIds = np.array([5,6,7,8])

    CleanDoubleElements(mesh1, preserveOriginalIds=True)
    if mesh1.GetNumberOfElements() != 3:
        raise Exception("ko")
    print(mesh1.GetElementsOfType(ElementNames.Tetrahedron_4).originalIds)

    return "ok"

def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_NodesPermutation,
    CheckIntegrity_Morphing,
    CheckIntegrity_ComputeFeatures,
    CheckIntegrity_ComputeSkin,
    CheckIntegrity_DeleteInternalFaces,
    CheckIntegrity_CleanDoubleNodes,
    CheckIntegrity_CleanLonelyNodes,
    CheckIntegrity_CleanEmptyTags,
    CheckIntegrity_AddTagPerBody,
    CheckIntegrity_RigidBodyTransformation,
    CheckIntegrity_ComputeRigidBodyTransformationBetweenTwoSetOfPoints,
    CheckIntegrity_CleanDoubleElements,
    CheckIntegrity_CopyElementTags
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
