# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

import BasicTools.Containers.UnstructuredMeshTools as UMT
import BasicTools.Containers.ConstantRectilinearMesh as CRM
import BasicTools.Containers.ElementNames as EN

def GetElementsCenters(mesh, dim=None):

    if mesh.IsConstantRectilinear():
        mesh.GenerateFullConnectivity()

    if dim is None:
        res = np.empty((mesh.GetNumberOfElements(),3) )
    else:
        res = np.empty((mesh.GetNumberOfElements(dim),3) )

    cpt= 0
    pos = mesh.GetPosOfNodes()
    for elementName in mesh.elements:
        if dim is not None and EN.dimension[elementName] != dim:
            continue
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
