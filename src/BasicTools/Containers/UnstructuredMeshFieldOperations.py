# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Tuple, Optional, Union

import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType, ArrayLike
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshCreationTools import QuadToLin
from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByElementFilter
from BasicTools.Linalg.Transform import Transform
from BasicTools.FE.Fields.FEField import FEField

def ApplyRotationMatrixTensorField(fields,fieldstoTreat, baseNames=["v1","v2"],inplace=False,prefix="new_",inverse=False):
    nbentries = fields[fieldstoTreat[0][0]].shape[0]

    bs =  Transform()
    bs.keepNormalised = True
    v1 = fields[baseNames[0]]
    v2 = fields[baseNames[1]]
    tempdata = np.zeros((nbentries,3,3))
    for i in range(3):
        for j in range(3):
            tempdata[:,i,j] = fields[fieldstoTreat[i][j]][:]

    for i in range(nbentries):
        bs.SetFirst(v1[i,:])
        bs.SetSecond(v2[i,:])
        if inverse:
            tempdata[i,:,:] = bs.ApplyInvTransformTensor(tempdata[i,:,:])
        else:
            tempdata[i,:,:] = bs.ApplyTransformTensor(tempdata[i,:,:])

    output = {}
    for i in range(3):
        for j in range(3):
            output[fieldstoTreat[i][j]]  = tempdata[:,i,j]

    output = {(prefix+x):v for x,v in output.items()}

    if inplace:
        fields.update(output)

    return output

def CopyFieldsFromOriginalMeshToTargetMesh(inmesh,outmesh):
    """ Function to copy fields for the original mesh to the
        derivated mesh (mesh generated from one operations)"""
    def Work(indic,outdic,oid):
        for k,d in indic.items():
            outdic[k] = d[oid]
    Work(inmesh.nodeFields,outmesh.nodeFields, outmesh.originalIDNodes )
    Work(inmesh.elemFields,outmesh.elemFields, outmesh.GetElementsOriginalIDs() )


def GetFieldTransferOp(inputField: FEField, targetPoints:ArrayLike, method:Union[str,None]=None, verbose:bool=False, elementfilter:Optional[ElementFilter]=None)-> Tuple[np.ndarray,np.ndarray]:
    """Compute the transfer operator from the inputField to the target points so:
    valueAtTargetPoints = op.dot(FEField.data)


        possibleMethods =["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"]
    possibleMethodsDict = {"Nearest":0, "Interp":1, "Extrap":2, "Clamp":3, "ZeroFill":4  }

    Parameters
    ----------
    inputField : FEField
        the FEField to be transferred
    targetPoints : ArrayLike
        Numpy array of the target points. Position to extract the values
    method : Union[str,None], optional
        A couple for the algorithm used when the point is inside the mesh :
            "Interp"  -> to use the interpolation of the FEField to extract the values
            "Nearest" -> to use the closest point to extract the values
        Possible options are:
            "Interp/Nearest", "Nearest/Nearest", "Interp/Clamp", "Interp/Extrap", "Interp/ZeroFill"
        If None is provided then "Interp/Clamp" is used
    verbose : bool, optional
        Print a progress bar, by default False
    elementfilter : Optional[ElementFilter], optional
        ElementFilter to extract the information from only a subdomain, by default None

    Returns
    -------
    Tuple [np.ndarray,np.ndarray]
        a tuple with 2 object containing:
            op, sparse matrix with the operator to make the transfer
            status: vector of ints with the status transfer for each target point :
               0: "Nearest"
               1: "Interp"
               2: "Extrap"
               3: "Clamp"
               4: "ZeroFill"
        return op, status

    """

    possibleMethods =["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"]
    possibleMethodsDict = {"Nearest":0, "Interp":1, "Extrap":2, "Clamp":3, "ZeroFill":4  }

    from scipy.spatial import cKDTree
    from scipy.sparse import coo_matrix
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

    if method is None:
        method = possibleMethods[2]
    elif method not in possibleMethods:
        raise(Exception(f"Method for transfert operator not know '{method}' possible options are : {possibleMethods}" ))

    insideMethod = possibleMethodsDict[method.split("/")[0]]
    outsideMethod = possibleMethodsDict[method.split("/")[1]]

    imeshdim = 0
    for i in range(4):
        if inputField.mesh.GetNumberOfElements(i):
            imeshdim = i

    originalmesh = inputField.mesh

    if elementfilter == None:
        elementfilter = ElementFilter(mesh = originalmesh, dimensionality=imeshdim)

    from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes
    imesh = ExtractElementsByElementFilter(originalmesh, elementfilter )
    CleanLonelyNodes(imesh)

    inodes = imesh.nodes

    numbering = inputField.numbering
    space = inputField.space
    inodes = imesh.nodes
    nbtp = targetPoints.shape[0]


    kdt = cKDTree(inodes)

    if insideMethod == 0 and  outsideMethod == 0:
        dist, ids = kdt.query(targetPoints)
        ids = imesh.originalIDNodes[ids]

        if numbering is None or numbering["fromConnectivity"]:
            cols = ids
        else:
            cols = [ numbering.GetDofOfPoint(pid) for pid in ids ]
        row = np.arange(nbtp)
        data = np.ones(nbtp)
        return coo_matrix((data, (row, cols)), shape=(nbtp , inodes.shape[0])), np.zeros(nbtp)

    # be sure to create the spaces
    for elementType, data in inputField.mesh.elements.items():
        space[elementType].Create()
        LagrangeSpaceGeo[elementType].Create()

        facestoTreat = [ElementNames.faces[elementType],
                        ElementNames.faces2[elementType],
                        ElementNames.faces3[elementType]]
        for tt in facestoTreat:
            for elementType2, num in tt:
                LagrangeSpaceGeo[elementType2].Create()

    # we build de Dual Coonectivity
    from BasicTools.Containers.UnstructuredMeshInspectionTools import GetDualGraphNodeToElement
    dualGraph,nused = GetDualGraphNodeToElement(imesh)

    from BasicTools.Containers.MeshTools import  GetElementsCenters

    centers = GetElementsCenters(imesh)
    kdtcenters = cKDTree(centers)

    # 30 to be sure to hold exa27 coefficients
    cols = np.empty(nbtp*30, dtype=int)
    rows = np.empty(nbtp*30, dtype=int)
    datas = np.empty(nbtp*30 )
    fillcpt = 0
    cood = [cols,rows,datas, fillcpt]

    def AddToOutput(l,col,row,dat,cood):
        fillcpt = cood[3]
        cood[0][fillcpt:fillcpt+l] = col
        cood[1][fillcpt:fillcpt+l] = row
        cood[2][fillcpt:fillcpt+l] = dat
        cood[3] += l

    def GetElement(imesh,enb):
        for name,data in imesh.elements.items():
            if enb < data.GetNumberOfElements():
                return originalmesh.elements[name], data.originalIds[enb], data, enb
            else:
                enb -= data.GetNumberOfElements()

        raise(Exception("Element not found"))

    distTP, idsTP = kdt.query(targetPoints)
    distTPcenters, idsTPcenters = kdtcenters.query(targetPoints)

    if verbose:
        from BasicTools.Helpers.ProgressBar import printProgressBar
        printProgressBar(0, nbtp, prefix = 'Building Transfer '+method+':', suffix = 'Complete', length = 50)
        verbosecpt = 0

    status = np.zeros(nbtp)
    ones = np.ones(50)
    for p in range(nbtp):
        if verbose:
            nvc = int(p/nbtp*1000)
            if verbosecpt != nvc:
                printProgressBar(p+1, nbtp, prefix = 'Building Transfer '+method+':', suffix = 'Complete', length = 50)
                verbosecpt = nvc

        TP = targetPoints[p,:]  # target point posicion
        CP =  idsTP[p]          # closest point posicion

        ## we use the closest element (in the sens of cells center )
        origial_data, lenb, imesh_data, imesh_elnb = GetElement(imesh,idsTPcenters[p])
        ## construct the potentialElements list (all element touching the closest element)
        potentialElements = []
        #Element connected to the closest point
        potentialElements.extend(dualGraph[idsTP[p],0:nused[idsTP[p]]])
        #Elements connected to the closest element (bases on the element center)
        for elempoint in  imesh_data.connectivity[imesh_elnb,:]:
            potentialElements.extend(dualGraph[elempoint,0:nused[elempoint]])
        potentialElements = np.unique(potentialElements)
        # compute distance to elements
        # for the moment we use the distance to the center, this gives a good estimate
        # of the order to check the elements
        dist = np.empty(len(potentialElements))
        for cpt,e in enumerate(potentialElements):
            data, lenb, imesh_data, imesh_elnb = GetElement(imesh,e)
            diff = centers[e,:]-TP
            dist[cpt] = diff.dot(diff)

        # order the element to test, closest element first
        index = np.argsort(dist)
        dist = dist[index]
        potentialElements = potentialElements[index]
        distmem = 1e10

        for cpt,e in enumerate(potentialElements):
            origial_data, lenb, imesh_data, imesh_elnb = GetElement(imesh,e)
            localnumbering = numbering[origial_data.elementType]
            localspace = space[origial_data.elementType]

            posnumbering = origial_data.connectivity
            posspace = LagrangeSpaceGeo[origial_data.elementType]
            coordAtDofs = originalmesh.nodes[posnumbering[lenb,:],:]

            inside, distv, bary, baryClamped = ComputeInterpolationExtrapolationsBarycentricCoordinates(TP, origial_data.elementType ,coordAtDofs, LagrangeSpaceGeo)

            #update the distance**2 with a *exact* distance
            dist[cpt] = dmin2 = distv.dot(distv)

            #compute shape function of the incomming space using the xi eta phi
            shapeFunc = localspace.GetShapeFunc(bary)
            shapeFuncClamped = localspace.GetShapeFunc(baryClamped)

            if inside:
                sF =  shapeFunc
                status[p] = 1
                break

            # store the best element (closest)
            if dist[cpt] < distmem:
                distmem = dist[cpt]
                if inside is None:
                    # non is error in the inverse mapping, only clamp is available
                    memshapeFunc = None
                else:
                    memshapeFunc = shapeFunc
                memshapeFuncClamped =  shapeFuncClamped
                memlenb = lenb
                memlocalnumbering = localnumbering
        else:
            # we are outside
            if distmem == 1e10 or outsideMethod == 0 or (outsideMethod == 2 and memshapeFunc is None):
                # not a single element found
                col = [CP]
                row = [p]
                dat = [1.]
                AddToOutput(len(col),col,row,dat,cood)
                continue

            if outsideMethod == 2 and memshapeFunc is not None:
                sF = memshapeFunc
                status[p] = 2
            elif outsideMethod == 3:
                sF = memshapeFuncClamped
                status[p] = 3
            else:
                status[p] = 4
                continue

            lenb = memlenb
            localnumbering = memlocalnumbering

        col = localnumbering[lenb,:]
        l = len(col)
        row = p*ones[0:l]
        dat = sF
        AddToOutput(l,col,row,dat,cood)

    if verbose:
        printProgressBar(nbtp, nbtp, prefix = 'Building Transfer '+method+':', suffix = 'Complete', length = 50)

    return coo_matrix((cood[2][0:cood[3]], (cood[1][0:cood[3]], cood[0][0:cood[3]])), shape=(nbtp , inputField.numbering["size"])), status

def ComputeInterpolationExtrapolationsBarycentricCoordinates(TP,elementType,coordAtDofs,posspace):
    if ElementNames.dimension[elementType]==0:
        distv = coordAtDofs[0,:] - TP
        inside = np.all(distv == 0)
        bary = np.array([])
        baryClamped = np.array([])
        return inside, distv, bary, baryClamped

    localspace = posspace[elementType]

    # we compute the baricentric coordinate of the target point (TP) on the current element
    inside, bary, baryClamped = ComputeBarycentricCoordinateOnElement(coordAtDofs,localspace,TP,elementType)
    if inside :
        # the point is inside, we compute the distance vector
        distv = localspace.GetShapeFunc(bary).dot(coordAtDofs) - TP
        return inside, distv, bary, baryClamped
    else :
        if inside is None:
            bary = np.zeros_like(bary)
        #compute distance to the points
        dist2 = np.sum((coordAtDofs - TP)**2,axis=1)
        minidx = np.argmin(dist2)
        mask = (-dist2 + dist2[minidx]) >= 0
        facestoTreat = [ElementNames.faces[elementType],
                        ElementNames.faces2[elementType],
                        ElementNames.faces3[elementType]]

        for faces in facestoTreat:
            if len(faces) == 0:
                break
            if faces[0][0] == ElementNames.Point_1:
                return False, coordAtDofs[minidx] - TP, bary, localspace.posN[minidx]
            faceInside, faceDistv, fbaryClamped =  ComputeInterpolationCoefficients(mask, TP, faces, coordAtDofs, posspace, localspace)
            if faceInside:
                return False, faceDistv, bary, fbaryClamped
        raise
    return None, 0, 0,0

def ComputeInterpolationCoefficients(mask, TP, elementsdata, coordAtDofs, posspace, localspace):
    ft = True
    for cpt, (faceElementType, faceLocalConnectivity) in enumerate(elementsdata):

        # work only on element touching the mask
        if not np.any(mask[faceLocalConnectivity]):
            continue

        faceInside, fbary, fbaryClamped = ComputeBarycentricCoordinateOnElement(coordAtDofs[faceLocalConnectivity,:],posspace[faceElementType],TP,faceElementType)
        faceDistv = posspace[faceElementType].GetShapeFunc(fbary).dot(coordAtDofs[faceLocalConnectivity,:]) - TP
        faceDist = faceDistv.dot(faceDistv)
        if faceInside:
            if ft:
                ft =False
                faceDistvMem = faceDistv
                faceDistMem = faceDist
                fbaryMem = fbary
                faceLocalConnectivityMem = faceLocalConnectivity
                faceElementTypeMem = faceElementType
            else:
                if faceDist < faceDistMem:
                    faceDistMem = faceDist
                    faceDistvMem = faceDistv
                    fbaryMem = fbary
                    faceLocalConnectivityMem = faceLocalConnectivity
                    faceElementTypeMem = faceElementType
    if ft:
        return False, None,None
    else:
        flocalspace = posspace[faceElementTypeMem]
        fshapeFunc = flocalspace.GetShapeFunc(fbaryMem)
        baryClamped = fshapeFunc.dot(localspace.posN[faceLocalConnectivityMem,:])
        return True, faceDistvMem,  baryClamped

def ComputeShapeFunctionsOnElement(coordAtDofs,localspace,localnumbering,point,faceElementType):
    inside, xietaphi, xietaphiClamped = ComputeBarycentricCoordinateOnElement(coordAtDofs,localspace,point,faceElementType)
    N = localspace.GetShapeFunc(xietaphi)
    NClamped = localspace.GetShapeFunc(xietaphiClamped)
    return inside, N , NClamped

def ddf(f,xietaphi,dN,GetShapeFuncDerDer,coordAtDofs,linear):
    dNX = dN.dot(coordAtDofs)
    res = dNX.dot(dNX.T)
    # After some investigation the Newton is more stable using only the first part
    # of the hessian we comment this part only for historic reasons. of if someone
    # can implement a better (without bugs) version.
    return res

    if linear :
        return res



    ddN = GetShapeFuncDerDer(xietaphi)
    # we work for every coordinate
    for ccpt in range(coordAtDofs.shape[1]):
        #error of the component ccpt
        fx = f[ccpt]
        coordx = coordAtDofs[:,ccpt]
        # we build the derivative of the pos field with respect of xi chi
        # d2f_i/dxidchi * pos[i,ccpt]
        pp = [ ddNi.dot(xi) for ddNi,xi in zip(ddN,coordx ) ]
        p = np.add.reduce(pp)
        res -= fx*p

    return res

def df(f,dN,coordAtDofs):
    dNX = dN.dot(coordAtDofs)
    res = -f.dot(dNX.T)
    return res

def vdet(A):
    detA = A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) -\
           A[0, 1] * (A[2, 2] * A[1, 0] - A[2, 0] * A[1, 2]) +\
           A[0, 2] * (A[1, 0] * A[2, 1] - A[2, 0] * A[1, 1])
    return detA

def hdinv(A):
    invA = np.zeros_like(A)
    detA = vdet(A)

    invA[0, 0] = (-A[1, 2] * A[2, 1] +
                  A[1, 1] * A[2, 2]) / detA
    invA[1, 0] = (A[1, 2] * A[2, 0] -
                  A[1, 0] * A[2, 2]) / detA
    invA[2, 0] = (-A[1, 1] * A[2, 0] +
                  A[1, 0] * A[2, 1]) / detA
    invA[0, 1] = (A[0, 2] * A[2, 1] -
                  A[0, 1] * A[2, 2]) / detA
    invA[1, 1] = (-A[0, 2] * A[2, 0] +
                  A[0, 0] * A[2, 2]) / detA
    invA[2, 1] = (A[0, 1] * A[2, 0] -
                  A[0, 0] * A[2, 1]) / detA
    invA[0, 2] = (-A[0, 2] * A[1, 1] +
                  A[0, 1] * A[1, 2]) / detA
    invA[1, 2] = (A[0, 2] * A[1, 0] -
                  A[0, 0] * A[1, 2]) / detA
    invA[2, 2] = (-A[0, 1] * A[1, 0] +
                  A[0, 0] * A[1, 1]) / detA
    return invA

def inv22(A):
    a = A[0,0]
    b = A[0,1]
    c = A[1,0]
    d = A[1,1]
    invA = np.zeros_like(A)
    invA[0,0] = d
    invA[0,1] = -b
    invA[1,0] = -c
    invA[1,1] = a
    invA /= (a*d-b*c)
    return invA

def ComputeBarycentricCoordinateOnElement(coordAtDofs,localspace,targetPoint,elementType):

    linear = ElementNames.linear[elementType]
    spacedim = localspace.GetDimensionality()

    xietaphi = np.array([0.5]*spacedim)
    N = localspace.GetShapeFunc(xietaphi)
    currentPoint = N.dot(coordAtDofs)
    f = targetPoint - currentPoint
    for x in range(10):
        dN = localspace.GetShapeFuncDer(xietaphi)
        df_num = df(f,dN,coordAtDofs)
        H = ddf(f, xietaphi, dN, localspace.GetShapeFuncDerDer, coordAtDofs, linear)
        if spacedim == 2:
            dxietaphi = inv22(H).dot(df_num)
        elif spacedim == 3:
            dxietaphi = hdinv(H).dot(df_num)
        else:
            dxietaphi = df_num/H[0,0]
        xietaphi -= dxietaphi

        # if the cell is linear only one iteration is needed
        if linear :
            break

        N = localspace.GetShapeFunc(xietaphi)
        f = targetPoint - N.dot(coordAtDofs)

        if normsquared(dxietaphi) < 1e-3 and normsquared(f) < 1e-3 :
            break
    else:
        return None, xietaphi,localspace.ClampParamCoorninates(xietaphi)

    xichietaClamped = localspace.ClampParamCoorninates(xietaphi)
    # we treat very closes point as inside
    inside = normsquared(xichietaClamped-xietaphi) < 1e-5
    return inside, xietaphi, xichietaClamped

def normsquared(v):
    return sum([x*x for x in v])

def TransportPosToPoints(imesh,points,method="Interp/Clamp",verbose=None):
    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, offset, NGauss = PrepareFEComputation(imesh,numberOfComponents=1)
    field = FEField("",mesh=imesh,space=space,numbering=numberings[0])
    op, status = GetFieldTransferOp(field,points,method=method, verbose=verbose, elementfilter = ElementFilter(imesh) )
    return op.dot(imesh.nodes)

def TransportPos(imesh,tmesh,tspace,tnumbering,method="Interp/Clamp",verbose=None):
    """
    Function to transport the position from the input mesh (imesh)
    to a target FEField target mesh, target space, tnumbering
    """
    from BasicTools.FE.Fields.FEField import FEField
    pos = TransportPosToPoints(imesh, tmesh.nodes,method=method,verbose=verbose)
    names = ["x","y","z"]
    posFields = np.array([ FEField("ipos_"+names[x],tmesh, space=tspace, numbering=tnumbering,data=pos[:,x]) for x in [0,1,2] ])
    return posFields


def PointToCellData(mesh,pointfield,dim=None):

    nbelemtns = 0
    filt = ElementFilter(mesh,dimensionality=dim)
    for name,data,ids in filt:
        nbelemtns +=  len(ids)

    if len(pointfield.shape) == 2:
        ncols = pointfield.shape[1]
        res = np.zeros((nbelemtns,ncols),dtype=float)
    else:
        ncols  = 1
        res = np.zeros((nbelemtns),dtype=float)

    cpt = 0
    for name,data,ids in filt:
        if len(pointfield.shape) == 1:
            valAtCenter = (np.sum(pointfield[data.connectivity],axis=1)/data.GetNumberOfNodesPerElement()).flatten()
            res[cpt:cpt+data.GetNumberOfElements()] = valAtCenter
        else:
            for i in range(ncols):
                valAtCenter = (np.sum(pointfield[data.connectivity,i],axis=1)/data.GetNumberOfNodesPerElement()).flatten()
                res[cpt:cpt+data.GetNumberOfElements(),i] = valAtCenter
        cpt += len(ids)
    return res

def QuadFieldToLinField(quadMesh, quadField, linMesh = None):

    if linMesh == None:
        linMesh = QuadToLin(quadMesh)

    extractIndices = np.arange(quadMesh.GetNumberOfNodes())[linMesh.originalIDNodes]

    return(quadField[extractIndices])

def GetValueAtPosLinearSymplecticMesh(fields,mesh,constantRectilinearMesh, verbose=False):
    """
    Works only for linear symplectic meshes
    """
    import math
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
    from BasicTools.FE.DofNumbering import ComputeDofNumbering

    if verbose:
        from BasicTools.Helpers.ProgressBar import printProgressBar

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity =True)

    mesh.ComputeBoundingBox()

    origin = constantRectilinearMesh.GetOrigin()
    spacing = constantRectilinearMesh.GetSpacing()
    dimensions = constantRectilinearMesh.GetDimensions()

    kmin, kmax = 0, 1

    shapeRes = [fields.shape[0]]
    for d in dimensions:
        shapeRes.append(d)
    result = np.zeros(tuple(shapeRes))
    mask = np.zeros(tuple(dimensions))

    numbering = ComputeDofNumbering(mesh, LagrangeSpaceGeo,fromConnectivity = True)

    for name, data in mesh.elements.items():
        if (ElementNames.dimension[name] == mesh.GetDimensionality() and ElementNames.linear[name] == True):

            if verbose:
                printProgressBar(0, data.GetNumberOfElements(), prefix = 'Field transfert element '+name+':', suffix = 'Complete', length = 50)

            for el in range(data.GetNumberOfElements()):
                if verbose:
                    printProgressBar(el, data.GetNumberOfElements(), prefix = 'Field transfert element '+name+':', suffix = 'Complete', length = 50)

                localNumbering = numbering[name][el,:]

                localNodes = mesh.nodes[data.connectivity[el,:]]
                nodesCoords = localNodes - mesh.boundingMin
                localBoundingMin = np.amin(localNodes, axis=0)
                localBoundingMax = np.amax(localNodes, axis=0)

                imin, imax = max(int(math.floor((localBoundingMin[0]-origin[0])/spacing[0])),0),min(int(math.floor((localBoundingMax[0]-origin[0])/spacing[0])+1),dimensions[0])
                jmin, jmax = max(int(math.floor((localBoundingMin[1]-origin[1])/spacing[1])),0),min(int(math.floor((localBoundingMax[1]-origin[1])/spacing[1])+1),dimensions[1])

                if mesh.GetDimensionality()>2:
                    kmin, kmax = min(int(math.floor((localBoundingMin[2]-origin[2])/spacing[2])),0),max(int(math.floor((localBoundingMax[2]-origin[2])/spacing[2])+1),dimensions[2])

                """imin, imax = math.floor((localBoundingMin[0])/spacing[0]),math.floor((localBoundingMax[0])/spacing[0])+1
                jmin, jmax = math.floor((localBoundingMin[1])/spacing[1]),math.floor((localBoundingMax[1])/spacing[1])+1

                if mesh.GetDimensionality()>2:
                    kmin, kmax = math.floor((localBoundingMin[2])/spacing[2]),math.floor((localBoundingMax[2])/spacing[2])+1"""

                for i in range(imin,imax):
                    for j in range(jmin,jmax):
                        for k in range(kmin,kmax):
                            if mesh.GetDimensionality()==2:
                                point = np.asarray([i*spacing[0],j*spacing[1]]) + origin - mesh.boundingMin
                            else:
                                point = np.asarray([i*spacing[0],j*spacing[1],k*spacing[2]]) + origin - mesh.boundingMin

                            rhs = np.hstack((point,np.asarray([1.])))
                            M = np.vstack((nodesCoords.T,np.ones(ElementNames.numberOfNodes[name])))
                            qcoord = np.linalg.solve(M,rhs)        # coordonnees barycentriques pour evaluer les fct de forme
                            if (qcoord>=-1.e-12).all() == True:
                                if mesh.GetDimensionality()==2:
                                    mask[i,j] = 1.
                                    for l in range(fields.shape[0]):
                                        result[l,i,j] = np.dot(qcoord,fields[l][localNumbering])
                                else:
                                    mask[i,j,k] = 1.
                                    for l in range(fields.shape[0]):
                                        result[l,i,j,k] = np.dot(qcoord,fields[l][localNumbering])
    return result, mask

#------------------------- CheckIntegrity ------------------------
def CheckIntegrityApplyRotationMatrixTensorField(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars
    inputmesh = CreateUniformMeshOfBars(0,1,10)
    nn = inputmesh.GetNumberOfNodes()
    inputmesh.nodeFields["Sxx"] = np.ones(nn)
    inputmesh.nodeFields["Syy"] = np.ones(nn)*2
    inputmesh.nodeFields["Szz"] = np.ones(nn)*3
    inputmesh.nodeFields["Sxy"] = np.ones(nn)*0
    inputmesh.nodeFields["Sxz"] = np.ones(nn)*0
    inputmesh.nodeFields["Syz"] = np.ones(nn)*0

    inputmesh.nodeFields["v1"] = np.zeros((nn,3))
    inputmesh.nodeFields["v1"][:,0] = 1
    inputmesh.nodeFields["v2"] = np.zeros((nn,3))
    inputmesh.nodeFields["v2"][:,2] = 1


    res = ApplyRotationMatrixTensorField(inputmesh.nodeFields,[["Sxx","Sxy","Sxz"],["Sxy","Syy","Syz"],["Sxz","Syz","Szz"]], baseNames=["v1","v2"],inplace=False,prefix="new_",inverse=False)

    for k,v in res.items():
        print(k,v)

    return "ok"

def RunTransfer(inputFEField,data,outmesh):
    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    WriteMeshToXdmf(tempdir+"GetFieldTransferOp_Original"+inputFEField.name+".xdmf",
                    inputFEField.mesh,
                    PointFields = [data],
                    PointFieldsNames = ["OriginalData"])

    PointFields = []
    PointFieldsNames = []
    for method in ["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap"]:
        op,status = GetFieldTransferOp(inputFEField,outmesh.nodes,method = method,verbose=True)
        newdata = op.dot(data)
        PointFieldsNames.append(method)
        PointFields.append(newdata)
        PointFieldsNames.append(method+"Status")
        PointFields.append(status)

    WriteMeshToXdmf(tempdir+"GetFieldTransferOp_TargetMesh"+inputFEField.name+".xdmf",
                    outmesh,
                    PointFields = PointFields,
                    PointFieldsNames = PointFieldsNames)

def CheckIntegrity1D(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars
    inputmesh = CreateUniformMeshOfBars(0,1,10)

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, _offset, _NGauss = PrepareFEComputation(inputmesh)

    b = CreateUniformMeshOfBars(-0.1,1.1,15)

    from BasicTools.FE.Fields.FEField import FEField
    inputFEField = FEField(name="",mesh=inputmesh,space=space,numbering=numberings[0])

    #for el,data in inputmesh.elements.items():
    #    print(data.connectivity)
    if GUI:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(nrows=2, ncols=3, constrained_layout=True)
        axis = axs.flat
    else:
        axis = [None]*5
    data = (inputmesh.nodes[:,0] -0.5)**2


    for method,ax in zip(["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"],axis):

        op = GetFieldTransferOp(inputFEField,b.nodes,method = method,verbose=True)[0]
        result = op.dot(data)
        if GUI:
            ax.plot(inputmesh.nodes[:,0],data,"X-",label="Original Data")
            ax.plot(b.nodes[:,0],result,"o:",label=method)
            legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')

            ax.set(xlabel='x', ylabel='data', title=method)

    if GUI:
        fig.savefig("test.png")
        plt.show()

    return "ok"

def CheckIntegrity2D(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    inputmesh = CreateSquare(dimensions = [5,5],origin=[0,0],spacing=[1./4,1./4.])

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, _offset, _NGauss = PrepareFEComputation(inputmesh)
    N = 10
    b = CreateSquare(dimensions = [N,N],origin=[-0.1,-0.1],spacing=[1.2/(N-1),1.2/(N-1)])

    from BasicTools.FE.Fields.FEField import FEField
    inputFEField = FEField(name="2DTo2D",mesh=inputmesh,space=space,numbering=numberings[0])
    x = inputmesh.nodes[:,0]
    y = inputmesh.nodes[:,1]
    data = (x -0.5)**2-y*0.5+x*y*0.25
    RunTransfer(inputFEField,data,b)
    return "ok"

def CheckIntegrity1DTo2D(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare


    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars
    inputmesh = CreateUniformMeshOfBars(0,1,10)
    inputmesh.nodes[:,1] = inputmesh.nodes[:,0]**2
    inputmesh.nodes = inputmesh.nodes[:,0:2].copy(order='C')

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, _offset, _NGauss = PrepareFEComputation(inputmesh)
    N = 10
    b = CreateSquare(dimensions = [N,N],origin=[-0.5,-0.5],spacing=[2/(N-1),2/(N-1)])

    from BasicTools.FE.Fields.FEField import FEField
    inputFEField = FEField(name="1DTo2D",mesh=inputmesh,space=space,numbering=numberings[0])

    x = inputmesh.nodes[:,0]
    y = inputmesh.nodes[:,1]
    data = (x -0.5)**2-y*0.5+x*y*0.25
    RunTransfer(inputFEField,data,b)
    return "ok"

def CheckIntegrity2DTo3D(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare, CreateCube

    N = 10
    inputmesh = CreateSquare(dimensions = [N,N],origin=[0,0],spacing=[1/(N-1),1/(N-1)])
    n  = inputmesh.nodes
    inputmesh.nodes = np.zeros((n.shape[0],3) )
    inputmesh.nodes[:,0:2] = n
    inputmesh.nodes[:,2] = n[:,0]**3

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, _offset, _NGauss = PrepareFEComputation(inputmesh)

    from BasicTools.FE.Fields.FEField import FEField
    inputFEField = FEField(name="2DTo3D",mesh=inputmesh,space=space,numbering=numberings[0])
    N = 10
    b = CreateCube(dimensions = [N,N,N],origin=[-0.5]*3,spacing=[2/(N-1),2/(N-1),2/(N-1)])

    x = inputmesh.nodes[:,0]
    y = inputmesh.nodes[:,1]
    data = (x -0.5)**2-y*0.5+x*y*0.25
    RunTransfer(inputFEField,data,b)
    return "ok"


def CheckIntegrity1DSecondOrderTo2D(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars
    inputmesh = CreateUniformMeshOfBars(0,1,11,secondOrder=True)
    inputmesh.nodes[:,1] = inputmesh.nodes[:,0]**2
    inputmesh.nodes = inputmesh.nodes[:,0:2].copy(order='C')

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, _offset, _NGauss = PrepareFEComputation(inputmesh)
    N = 10
    #b = CreateSquare(dimensions = [N,N],origin=[-0.5,-0.5],spacing=[2/(N-1),2/(N-1)])
    #b = CreateSquare(dimensions = [N,N],origin=[-0.1,0.0],spacing=[1./(N-1),0.7/(N-1)])
    b = CreateSquare(dimensions = [N,N],origin=[-0.5,-0.5],spacing=[2/(N-1),2/(N-1)])

    from BasicTools.FE.Fields.FEField import FEField
    inputFEField = FEField(name="1DSecondTo2D",mesh=inputmesh,space=space,numbering=numberings[0])

    x = inputmesh.nodes[:,0]
    y = inputmesh.nodes[:,1]
    data = (x -0.5)**2-y*0.5+x*y*0.25
    RunTransfer(inputFEField,data,b)
    return "ok"


def CheckIntegrity_GetValueAtPosLinearSymplecticMesh(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    points = [[-0.5,-0.5,-0.5],[2.5,-0.5,-0.5],[-0.5,2.5,-0.5],[-0.5,-0.5,2.5],[2.5,2.5,2.5]]
    tets = [[0,1,2,3]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    recMesh = ConstantRectilinearMesh()
    recMesh.SetDimensions([5,5,5])
    recMesh.SetSpacing([1, 1, 1])
    recMesh.SetOrigin([-1, -1, -1])

    #from BasicTools.IO.GeofWriter import WriteMeshToGeof
    #WriteMeshToGeof("mesh.geof", mesh)
    #WriteMeshToGeof("recMesh.geof", recMesh)

    res = GetValueAtPosLinearSymplecticMesh(np.array([np.arange(mesh.GetNumberOfNodes())]),mesh,recMesh)
    """import matplotlib
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6, 3.2))
    plt.pcolor(res[1,:,:].transpose(), cmap=None)
    plt.colorbar(orientation='vertical')
    plt.show()"""

    return "OK"

def CheckIntegrity_PointToCellData(GUI = False):
    myMesh = UnstructuredMesh()
    myMesh.nodes = np.array([[0,0,0],[1,0,0],[2,0,0]] ,dtype=PBasicFloatType)
    myMesh.originalIDNodes = np.array([0,1,2] ,dtype=int)
    tag = myMesh.GetNodalTag("linPoints")
    tag.AddToTag(0)
    tag.AddToTag(1)
    tag.AddToTag(2)
    import BasicTools.Containers.ElementNames as ElementNames
    elements = myMesh.GetElementsOfType(ElementNames.Bar_2)
    elements.AddNewElement([0,1],3)
    elements.AddNewElement([1,2],4)

    myMesh.AddElementToTagUsingOriginalId(3,'LinElements')
    myMesh.AddElementToTagUsingOriginalId(4,'LinElements')
    myMesh.PrepareForOutput()
    print(myMesh)
    res = PointToCellData(myMesh,np.array([[-2,2,4]]).T)
    ExactData = np.array([[0,3]], dtype=float).T
    print (res - ExactData)
    if (res - ExactData).any() :
        return ("Error CheckIntegrity_PointToCellData")
    return "ok"

def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_GetValueAtPosLinearSymplecticMesh,
    CheckIntegrity_PointToCellData,
    CheckIntegrity1DSecondOrderTo2D,
    CheckIntegrity1DTo2D,
    CheckIntegrity1D,
    CheckIntegrity2D,
    CheckIntegrity2DTo3D,
    CheckIntegrityApplyRotationMatrixTensorField
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
