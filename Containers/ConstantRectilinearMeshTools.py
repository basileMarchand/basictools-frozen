# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import coo_matrix

from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

def GetSubSuperMesh(inputmesh, _newDimensions):
    newDimensions = np.array(_newDimensions,dtype=np.int)
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

    nps = 3
    nps3 = nps**inputmesh.GetDimensionality()
    oldToNewVals = np.zeros((destinationmesh.GetNumberOfNodes(),nps3))
    oldToNewIK = np.zeros((destinationmesh.GetNumberOfNodes(),nps3), dtype=np.int_)
    oldToNewJK = np.zeros((destinationmesh.GetNumberOfNodes(),nps3), dtype=np.int_)

    for i in  range(destinationmesh.GetNumberOfElements()):
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

    oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destinationmesh.GetNumberOfElements(), inputmesh.GetNumberOfElements())).tocsc()

    return oldToNew


#------------------------------------------------------------------------------
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
    return  "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(True))
    print("done")
