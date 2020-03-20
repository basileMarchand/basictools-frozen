# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       
import numpy as np

from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
import BasicTools.Containers.ConstantRectilinearMesh as CRM
import BasicTools.Containers.ElementNames as EN

def IsClose(mesh1,mesh2):
    if not np.all(np.isclose(mesh1.nodes,mesh2.nodes)):
        print("nodes not equal")
        return False

    for tag1 in mesh1.nodesTags:
        tag2 = mesh2.nodesTags[tag1.name]
        if not np.all(np.isclose(tag1.GetIds(),tag2.GetIds())):
            print("Nodal tag  "+ str(tag1.name) + " not equal")
            return False

    for name,data1 in mesh1.nodeFields.items():
        data2 = mesh2.nodeFields[name]
        if not np.all(np.isclose(data1,data2)):
            print("Field "+ str(name) + " not equal")
            return False

    for name, data1 in mesh1.elements.items():
        data2 = mesh2.elements[name]
        if not np.all(np.isclose(data1.connectivity,data2.connectivity)):
            print("Connectivity for  "+ str(name) + " not equal")
            return False
        for tag1 in data1.tags:
            tag2 = data2.tags[tag1.name]
            if not np.all(np.isclose(tag1.GetIds(),tag2.GetIds())):
                print("Tag " + str(tag1.name) + " is not equal for element" + str(name))
                return False

    for name,data1 in mesh1.elemFields.items():
        data2 = mesh2.elemFields[name]
        if not np.all(np.isclose(data1,data2)):
            print("Field "+ str(name) + " not equal")
            return False

    return True

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

    mesh1 = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
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
