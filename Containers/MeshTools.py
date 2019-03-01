# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

import BasicTools.Containers.UnstructuredMeshTools as UMT
import BasicTools.Containers.ConstantRectilinearMesh as CRM
import BasicTools.Containers.ElementNames as EN



def GetElementsCenters(mesh=None,nodes=None,elements=None, dim=None):
    # this function is used in the filters implementation
    # no Filter can appear in this implementation

    if mesh is not None and elements is not None:
        raise(Exception("Cant trat mesh and element at the same time" ) )


    def traiteElements(nod,els):

        connectivity = els.connectivity
        localRes = np.zeros((els.GetNumberOfElements(),3) )
        for i in range(nod.shape[1]):
            localRes[:,i] += np.sum(nod[connectivity,i],axis=1)
        localRes /= connectivity.shape[1]
        return localRes

    if mesh is not None:
        if dim is None:
            res = np.empty((mesh.GetNumberOfElements(),3) )
        else:
            res = np.empty((mesh.GetNumberOfElements(dim),3) )

        cpt= 0
        for elementName,data in mesh.elements.items():
            if dim is not None and EN.dimension[elementName] != dim:
                continue
            res[cpt:cpt+data.GetNumberOfElements()] = traiteElements(mesh.nodes,data)
            cpt += data.GetNumberOfElements()
    else:
        return traiteElements(nodes,elements)



    pos = mesh.GetPosOfNodes()


    return res


def CheckIntegrity_GetCellCenters():

    mesh1 = UMT.CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res = GetElementsCenters(mesh1)
    print(res)

    mesh2 = CRM.ConstantRectilinearMesh()
    mesh2.SetDimensions([2,3,2]);
    mesh2.SetSpacing([1, 1, 1]);
    mesh2.GetPosOfNodes()
    res = GetElementsCenters(mesh=mesh2)
    print(res)

    return "ok"


def CheckIntegrity():

    CheckIntegrity_GetCellCenters()

    return "OK"



if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
