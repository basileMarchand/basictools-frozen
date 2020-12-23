# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshCreationTools import QuadToLin
from BasicTools.Containers.UnstructuredMeshInspectionTools import  ExtractElementByDimensionalityNoCopy

def GetFieldTransferOp(inputField,targetPoints,method=None,verbose=False):
    """
        op : sparcematrix to do the transfer
        status: vector with the status of the transfer (1:interp, 2: extrap, 3:clamp, 4:ZeroFill, 0:nearest)
        return op, status
    """
    possibleMethods =["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"]

    if method is None:
        method = possibleMethods[0]
    elif method not in possibleMethods:
        raise(Exception("Method for transfert operator not know '{}' possible options are : {}".format(method,possibleMethods) ))

    imeshdim = 0
    for i in range(4):
        if inputField.mesh.GetNumberOfElements(i):
            imeshdim = i

    imesh = ExtractElementByDimensionalityNoCopy(inputField.mesh,imeshdim)
    numbering = inputField.numbering
    space = inputField.space
    inodes = imesh.nodes

    kdt = cKDTree(inodes)

    nbtp = targetPoints.shape[0]

    if method == "Nearest/Nearest" :
        dist, ids = kdt.query(targetPoints)

        if numbering is None or numbering["fromConnectivity"]:
            cols = ids
        else:
            cols = [ numbering["almanac"][('P',pid,None)] for pid in ids ]
        row = np.arange(nbtp)
        data = np.ones(nbtp)
        return coo_matrix((data, (row, cols)), shape=(nbtp , inodes.shape[0])), np.zeros(nbtp)

    # we build de Dual Coonectivity
    from BasicTools.Containers.UnstructuredMeshInspectionTools import GetDualGraphNodeToElement
    dualGraph,nused = GetDualGraphNodeToElement(imesh)

    from BasicTools.Containers.MeshTools import  GetElementsCenters

    centers = GetElementsCenters(imesh)

    # 30 to be sure to hold exa27 coefficients
    cols = np.empty(nbtp*30, dtype=int )
    rows = np.empty(nbtp*30, dtype=int )
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
                return data, enb
            else:
                enb -= data.GetNumberOfElements()

        raise(Exception("Element not found"))

    distTP, idsTP = kdt.query(targetPoints)

    if verbose:
        from BasicTools.Helpers.ProgressBar import printProgressBar
        printProgressBar(0, nbtp, prefix = 'Building Transfer '+method+':', suffix = 'Complete', length = 50)

    status = np.zeros(nbtp)
    ones = np.ones(50)
    for p in range(nbtp):
        if verbose:
            printProgressBar(p+1, nbtp, prefix = 'Building Transfer '+method+':', suffix = 'Complete', length = 50)

        TP = targetPoints[p,:]
        CP =  idsTP[p]
        LPE = nused[CP]
        potentialElements = dualGraph[CP,0:LPE]

        # compute distance to elements
        # for the moment we use the distance to the center, this gives a good estimate
        # of the order to check the elements
        dist = np.empty(LPE)
        for cpt,e in enumerate(potentialElements):
            data, lenb = GetElement(imesh,e)
            diff = centers[e,:]-TP
            dist[cpt] = diff.dot(diff)

        # order the element to test, closest element first
        index = np.argsort(dist)
        dist = dist[index]
        potentialElements = potentialElements[index]
        distmem = 1e10


        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
        for cpt,e in enumerate(potentialElements):
            data, lenb = GetElement(imesh,e)
            localnumbering = numbering[data.elementType]
            localspace = space[data.elementType]
            localspace.Create()

            posnumbering = data.connectivity
            posspace = LagrangeSpaceGeo[data.elementType]
            localspace.Create()
            coordAtDofs = imesh.nodes[posnumbering[lenb,:],:]

            #inside, shapeFunc, shapeFuncClamped = ComputeShapeFunctionsOnElement(coordAtDofs,localspace ,localnumbering,TP)

            inside, bary, baryClamped = ComputeBarycentricCoordinateOnElement(coordAtDofs,posspace ,posnumbering,TP,ElementNames.linear[data.elementType])

            if inside is None:
                continue

            posShapeFunc = posspace.GetShapeFunc(bary)
            posShapeFuncClamped = posspace.GetShapeFunc(baryClamped)

            #update the distance**2 with a *exact* distance
            distElemToTP =  posShapeFuncClamped.dot(coordAtDofs) -TP
            dist[cpt] = distElemToTP.dot(distElemToTP)

            #compute shape function of the incomming space using the xi eta phi
            shapeFunc = localspace.GetShapeFunc(bary)
            shapeFuncClamped = localspace.GetShapeFunc(baryClamped)

            if inside:
                sF =  shapeFunc
                status[p] = 1
                break

            # store the best elements (closest)
            if dist[cpt] < distmem:
                distmem = dist[cpt]
                memshapeFunc = shapeFunc
                memshapeFuncClamped =  shapeFuncClamped
                memdata = data
                memlenb = lenb
                memlocalnumbering = localnumbering
        else:
            if distmem == 1e10 or method == "Interp/Nearest":
                col = [CP]
                row = [p]
                dat = [1.]
                AddToOutput(len(col),col,row,dat,cood)
                #status[p] = 0
                continue
            elif method == "Interp/Clamp":
                sF = memshapeFuncClamped
                status[p] = 3
            elif method == "Interp/Extrap":
                # we use the shapeFunction in extrapolation
                sF = memshapeFunc
                status[p] = 2
            else:
                status[p] = 4
                continue

            data = memdata
            lenb = memlenb
            localnumbering = memlocalnumbering

        col = localnumbering[lenb,:]
        l = len(col)
        row = p*ones[0:l]
        dat = sF
        AddToOutput(l,col,row,dat,cood)

    return coo_matrix((cood[2][0:cood[3]], (cood[1][0:cood[3]], cood[0][0:cood[3]])), shape=(nbtp , inputField.numbering["size"])), status





def ComputeShapeFunctionsOnElement(coordAtDofs,localspace,localnumbering,point):
    inside, xietaphi, xietaphiClamped = ComputeBarycentricCoordinateOnElement(coordAtDofs,localspace,localnumbering,point)
    N = localspace.GetShapeFunc(xietaphi)
    NClamped = localspace.GetShapeFunc(xietaphiClamped)
    return inside, N , NClamped

def ddf(f,xietaphi,point,dN,GetShapeFuncDerDer,coordAtDofs,linear):
    dNX = dN.dot(coordAtDofs)
    res = 2*dNX.dot(dNX.T)

    if linear :
        for i in range(len(res) ):
            if res[i,i] == 0:
                res[i,i] = 1
        return res

    # TODO need Check this part not sure if is correct
    ddN = GetShapeFuncDerDer(xietaphi)
    for ccpt in range(coordAtDofs.shape[1]):
        for cpt,fct in enumerate(ddN):
           c = coordAtDofs[cpt,ccpt]
           res -= 2*f[ccpt]*(c*fct)

    for i in range(len(res) ):
        if res[i,i] == 0:
            res[i,i] = 1
    return res

def df(f,dN,coordAtDofs):
    dNX = dN.dot(coordAtDofs)
    res = 2*f.dot(-dNX.T)
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
    invA[0,0] = d/(a*d-b*c);
    invA[0,1] = -b/(a*d-b*c);
    invA[1,0] = -c/(a*d-b*c);
    invA[1,1] = a/(a*d-b*c);
    return invA

def ComputeBarycentricCoordinateOnElement(coordAtDofs,localspace,localnumbering,point,linear):

    xietaphi = np.array([0.5]*localspace.GetDimensionality())

    for x in range(10):
        N = localspace.GetShapeFunc(xietaphi)
        #print("x {} : {} ".format(x,xietaphi) )
        dN = localspace.GetShapeFuncDer(xietaphi)
        f = point - N.dot(coordAtDofs)
        df_num = df(f,dN,coordAtDofs)

        if df_num.dot(df_num) < 1e-10:
            break
        H = ddf(f,xietaphi,point,dN,localspace.GetShapeFuncDerDer,coordAtDofs,linear)
        if H.shape[0] == 2:
            xietaphi -=inv22(H).dot(df_num)
        elif H.shape[0] == 3:
            xietaphi -=hdinv(H).dot(df_num)
        else:
            xietaphi -= np.linalg.inv(H).dot(df_num)

    else:
        return None, xietaphi,localspace.ClampParamCoorninates(xietaphi)

    xichietaClamped = localspace.ClampParamCoorninates(xietaphi)
    inside = not np.any(xichietaClamped-xietaphi  )
    return inside, xietaphi, xichietaClamped


def TransportPos(imesh,tmesh,tspace,tnumbering):
    """
    Function to transport the position from the input mesh (imesh)
    to a target FEField target mesh, target space, tnumbering
    """
    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, offset, NGauss = PrepareFEComputation(imesh,numberOfComponents=1)
    field = FEField("",mesh=imesh,space=space,numbering=numberings[0])
    op,status = GetFieldTransferOp(field,tmesh.nodes,method="Interp/Clamp")
    pos = op.dot(field.mesh.nodes)
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

def GetValueAtPosLinearSymplecticMesh(fields,mesh,constantRectilinearMesh):
        """
        Works only for linear symplectic meshes
        """
        import math
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
        from BasicTools.FE.DofNumbering import ComputeDofNumbering

        numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity =True)

        mesh.ComputeBoundingBox()

        origin = constantRectilinearMesh.GetOrigin()
        spacing = constantRectilinearMesh.GetSpacing()
        dimensions = constantRectilinearMesh.GetDimensions()

        #print("origin =", origin)
        #print("spacing =", spacing)
        #print("dimensions =", dimensions)

        kmin, kmax = 0, 1

        shapeRes = [fields.shape[0]]
        for d in dimensions:
            shapeRes.append(d)
        result = np.zeros(tuple(shapeRes))
        mask = np.zeros(tuple(dimensions))
        for name, data in mesh.elements.items():
            #print("name =", name)
            #print("ElementNames.dimension[name] =", ElementNames.dimension[name])
            #print("mesh.GetDimensionality() =", mesh.GetDimensionality())
            #print("ElementNames.linear[name] =", ElementNames.linear[name])
            if (ElementNames.dimension[name] == mesh.GetDimensionality() and ElementNames.linear[name] == True):

                for el in range(data.GetNumberOfElements()):

                    localNumbering = numbering[name][el,:]

                    localNodes = mesh.nodes[data.connectivity[el,:]]
                    nodesCoords = localNodes - mesh.boundingMin
                    localBoundingMin = np.amin(localNodes, axis=0)
                    localBoundingMax = np.amax(localNodes, axis=0)
                    #print("nodesCoords =", nodesCoords)
                    #print("localBoundingMin =", localBoundingMin)
                    #print("localBoundingMax =", localBoundingMax)

                    numbering = ComputeDofNumbering(mesh, LagrangeSpaceGeo,fromConnectivity = True)

                    imin, imax = max(int(math.floor((localBoundingMin[0]-origin[0])/spacing[0])),0),min(int(math.floor((localBoundingMax[0]-origin[0])/spacing[0])+1),dimensions[0])
                    jmin, jmax = max(int(math.floor((localBoundingMin[1]-origin[1])/spacing[1])),0),min(int(math.floor((localBoundingMax[1]-origin[1])/spacing[1])+1),dimensions[1])
                    #print("imin, imax =", imin, imax)
                    #print("jmin, jmax =", jmin, jmax)

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
                                #print(point, rhs, qcoord)
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
    inputmesh.nodes = inputmesh.nodes[:,0:2]

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
    inputmesh.nodes = inputmesh.nodes[:,0:2]

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
    myMesh.nodes = np.array([[0,0,0],[1,0,0],[2,0,0]] ,dtype=np.float)
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
    CheckIntegrity2DTo3D
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
