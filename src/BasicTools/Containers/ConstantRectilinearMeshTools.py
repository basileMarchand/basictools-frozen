# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
from scipy.sparse import coo_matrix

from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
from BasicTools.NumpyDefs import PBasicIndexType

def GetSubSuperMesh(inputmesh, _newDimensions):
    newDimensions = np.array(_newDimensions,dtype=PBasicIndexType)
    ## to generate meshes with more or less elements in each directions
    ## return the mesh
    newSpac = ((inputmesh.GetDimensions()-1)*inputmesh.GetSpacing()).astype(float)/(newDimensions-1)

    res = type(inputmesh)(dim=inputmesh.GetDimensionality())
    res.SetSpacing(newSpac)
    res.SetDimensions(newDimensions)
    res.SetOrigin(inputmesh.GetOrigin())

    return res

def GetNodeTrasfertMatrix(inputmesh,destinationmesh):
    # newVector   = oldToNew * oldVector
    nbNodes = 2**inputmesh.GetDimensionality()
    oldToNewVals = np.zeros((destinationmesh.GetNumberOfNodes(),nbNodes))
    oldToNewIK = np.zeros((destinationmesh.GetNumberOfNodes(),nbNodes), dtype=np.int_)
    oldToNewJK = np.zeros((destinationmesh.GetNumberOfNodes(),nbNodes), dtype=np.int_)

    for i in range(destinationmesh.GetNumberOfNodes()):

        pos= destinationmesh.GetPosOfNode(i)
        el = inputmesh.GetElementAtPos(pos)
        coon = inputmesh.GetConnectivityForElement(el)
        xiChiEta = inputmesh.GetElementShapeFunctionsAtPos(el,pos)
        oldToNewVals[i,:] = xiChiEta
        oldToNewIK[i,:] = i
        oldToNewJK[i,:] = coon

    oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destinationmesh.GetNumberOfNodes(), inputmesh.GetNumberOfNodes())).tocsc()

    return oldToNew

def GetElementTrasfertMatrix(inputmesh, destinationmesh):

    if not isinstance(inputmesh, ConstantRectilinearMesh):
        raise Exception("First argument must be a ConstantRectilinearMesh")
    nps = 3
    nps3 = nps**inputmesh.GetDimensionality()
    oldToNewVals = np.zeros((destinationmesh.GetNumberOfNodes(),nps3))
    oldToNewIK = np.zeros((destinationmesh.GetNumberOfNodes(),nps3), dtype=np.int_)
    oldToNewJK = np.zeros((destinationmesh.GetNumberOfNodes(),nps3), dtype=np.int_)

    numberOfElementsDest = destinationmesh.GetNumberOfElements(dim = destinationmesh.GetElementsDimensionality() )
    numberOfElementsInp = inputmesh.GetNumberOfElements(dim = inputmesh.GetElementsDimensionality() )

    for i in  range(numberOfElementsDest):
        coon = destinationmesh.GetConnectivityForElement(i)

        n0pos = destinationmesh.GetPosOfNode(coon[0])
        cpt =0
        for cx in range(0,nps):
            for cy in range(0,nps):
                if inputmesh.GetDimensionality() == 3:
                    for cz in range(0,nps):
                        pos = n0pos + destinationmesh.GetSpacing()*([cx+0.5,cy+0.5,cz+0.5])/nps
                        el = inputmesh.GetElementAtPos(pos)
                        oldToNewVals[i,cpt] += 1./nps3
                        oldToNewIK[i,cpt] += i
                        oldToNewJK[i,cpt] += el
                        cpt +=1
                else:
                    pos = n0pos + destinationmesh.GetSpacing()*([cx+0.5,cy+0.5])/nps
                    el = inputmesh.GetElementAtPos(pos)
                    oldToNewVals[i,cpt] += 1./nps3
                    oldToNewIK[i,cpt] += i
                    oldToNewJK[i,cpt] += el
                    cpt +=1

    oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(numberOfElementsDest, numberOfElementsInp)).tocsc()

    return oldToNew


#------------------------------------------------------------------------------
def CreateSquare(dimensions=[2,2], origin=[-1.0,-1.0], spacing=[1.,1.]):
    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshModificationTools import ComputeSkin
    import BasicTools.Containers.ElementNames as EN

    myMesh = ConstantRectilinearMesh(dim=2)
    myMesh.SetDimensions(dimensions);
    myMesh.SetOrigin(origin);
    myMesh.SetSpacing(spacing);

    # coorners
    d = np.array(dimensions)-1
    s = spacing
    indexs = [[   0,   0,   0],
              [d[0],   0,   0],
              [   0,d[1],   0],
              [d[0],d[1],   0]]

    for n in indexs:
        idx = myMesh.GetMonoIndexOfNode(n)
        name = "x"  + ("0" if n[0]== 0 else "1" )
        name += "y" + ("0" if n[1]== 0 else "1" )
        myMesh.nodesTags.CreateTag(name,False).SetIds([idx])


    skin = ComputeSkin(myMesh)
    for name,data in skin.elements.items():
        myMesh.GetElementsOfType(name).Merge(data)
    #print(skin)

    quads = myMesh.GetElementsOfType(EN.Quadrangle_4)
    quads.GetTag("2D").SetIds(range(quads.GetNumberOfElements()))

    skin = myMesh.GetElementsOfType(EN.Bar_2)
    #face tags

    x = myMesh.GetPosOfNodes()[skin.connectivity,0]
    y = myMesh.GetPosOfNodes()[skin.connectivity,1]
    tol = np.min(spacing)/10

    skin.GetTag("X0").SetIds( np.where(np.sum(np.abs(x - origin[0]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("X1").SetIds( np.where(np.sum(np.abs(x - (origin[0]+d[0]*s[0]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Y0").SetIds( np.where(np.sum(np.abs(y - origin[1]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Y1").SetIds( np.where(np.sum(np.abs(y - (origin[1]+d[1]*s[1]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])


    myMesh.PrepareForOutput()
    return myMesh


def CreateMesh(dim):
    myMesh = ConstantRectilinearMesh(dim)
    myMesh.SetDimensions([2,]*dim);
    myMesh.SetSpacing([1, ]*dim);

    return myMesh
def CheckIntegrity_GetSubSuperMesh(dim):

    newmesh = GetSubSuperMesh(CreateMesh(dim),[3,]*dim)

def CheckIntegrity_GetNodeTrasfertMatrix(dim):
    mesh1 = CreateMesh(dim)
    mesh2 = GetSubSuperMesh(mesh1,[3,]*dim)

    TMatrix = GetNodeTrasfertMatrix(mesh1,mesh2)

def CheckIntegrity_GetElementTrasfertMatrix(dim):
    mesh1 = CreateMesh(dim)
    mesh2 = GetSubSuperMesh(mesh1,[3,]*dim)

    TMatrix = GetElementTrasfertMatrix(mesh1,mesh2)

def CheckIntegrity(GUI=False):
    for dim in [2,3]:
        CheckIntegrity_GetSubSuperMesh(dim)
        CheckIntegrity_GetNodeTrasfertMatrix(dim)
        CheckIntegrity_GetElementTrasfertMatrix(dim)
    CreateSquare()
    return  "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(True))
    print("done")
