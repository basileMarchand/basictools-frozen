# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshCreationTools import QuadToLin

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

    filt = ElementFilter(mesh,dimensionality=dim)
    cpt = 0
    for name,data,ids in filt:
        if len(pointfield.shape) == 1:
            valAtCenter = (np.sum(pointfield[data.connectivity],axis=1)/data.connectivity.shape[1]).flatten()
            res[cpt:cpt+data.GetNumberOfElements()] = valAtCenter
        else:
            for i in range(ncols):
                valAtCenter = (np.sum(pointfield[data.connectivity,i],axis=1)/data.connectivity.shape[1]).flatten()
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
        from BasicTools.FE.IntegrationsRules import Lagrange as Lagrange
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
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
