# -*- coding: utf-8 -*-

import numpy as np
import BasicTools.FE.UnstructuredMeshTools as UMT
import BasicTools.FE.ConstantRectilinearMesh as CRM

def GetElementsCenters(mesh):

    if mesh.IsConstantRectilinear():
        mesh.GenerateFullConnectivity()

    res = np.empty((mesh.GetNumberOfElements(),3) )

    cpt= 0
    pos = mesh.GetPosOfNodes()
    for elementName in mesh.elements:
        elements = mesh.elements[elementName]
        connectivity = elements.connectivity
        localRes = np.zeros((elements.GetNumberOfElements(),3) )

        for i in range(mesh.GetDimensionality()):
            localRes[:,i] += np.sum(pos[connectivity,i],axis=1)
        localRes /= connectivity.shape[1]

        res[cpt:cpt+elements.GetNumberOfElements()] = localRes
        cpt +=elements.GetNumberOfElements()
    return res


def CheckIntegrity_GetCellCenters():

    mesh1 = UMT.CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res = GetElementsCenters(mesh1)
    print(res)

    mesh2 = CRM.ConstantRectilinearMesh()
    mesh2.SetDimensions([2,3,2]);
    mesh2.SetSpacing([1, 1, 1]);
    res = GetElementsCenters(mesh2)
    print(res)

    return "ok"


def CheckIntegrity():

    CheckIntegrity_GetCellCenters()

    return "OK"



if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
