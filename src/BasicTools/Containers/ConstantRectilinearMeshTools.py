# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
from scipy.sparse import coo_matrix

from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType, ArrayLike

def GetSubSuperMesh(inputmesh:ConstantRectilinearMesh, newDimensions:ArrayLike)-> ConstantRectilinearMesh:
    """Generate a new ConstantRectilinearMesh with the same volume (origin and size )
    but with different number of elements

    Parameters
    ----------
    inputmesh : ConstantRectilinearMesh
        the input mesh to recover the origin and size
    newDimensions : ArrayLike
        the dimensions of the new mesh

    Returns
    -------
    ConstantRectilinearMesh
        the new mesh
    """
    newDimensions = np.array(newDimensions, dtype=PBasicIndexType)
    ## to generate meshes with more or less elements in each directions
    ## return the mesh
    newSpac = ((inputmesh.GetDimensions()-1)*inputmesh.GetSpacing()).astype(float)/(newDimensions-1)

    res = type(inputmesh)(dim=inputmesh.GetDimensionality())
    res.SetSpacing(newSpac)
    res.SetDimensions(newDimensions)
    res.SetOrigin(inputmesh.GetOrigin())

    return res

def GetNodeTransferMatrix(inputMesh:ConstantRectilinearMesh, destinationMesh:ConstantRectilinearMesh)-> coo_matrix:
    """return the transfer operator to transfer a field defined on the nodes of the input mesh (inputMesh)
    to the output mesh. Data is extrapolated outsize if a point of the destination mesh lies outside the input mesh

    Parameters
    ----------
    inputMesh : ConstantRectilinearMesh
        the mesh defining the support of the point field to be transferred
    destinationMesh : ConstantRectilinearMesh
        destination mesh

    Returns
    -------
    np.ndarray
        OP: The transfer operation (sparse matrix in coo)
        field in the new mesh = OP.dot(field defined in the input Mesh)
    """
    # newVector   = oldToNew * oldVector
    nbNodes = 2**inputMesh.GetDimensionality()
    oldToNewVals = np.zeros((destinationMesh.GetNumberOfNodes(),nbNodes))
    oldToNewIK = np.zeros((destinationMesh.GetNumberOfNodes(),nbNodes), dtype=np.int_)
    oldToNewJK = np.zeros((destinationMesh.GetNumberOfNodes(),nbNodes), dtype=np.int_)

    for i in range(destinationMesh.GetNumberOfNodes()):

        pos= destinationMesh.GetPosOfNode(i)
        el = inputMesh.GetElementAtPos(pos)
        coon = inputMesh.GetConnectivityForElement(el)
        xiChiEta = inputMesh.GetElementShapeFunctionsAtPos(el,pos)
        oldToNewVals[i,:] = xiChiEta
        oldToNewIK[i,:] = i
        oldToNewJK[i,:] = coon

    oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destinationMesh.GetNumberOfNodes(), inputMesh.GetNumberOfNodes())).tocsc()

    return oldToNew

def GetElementTransferMatrix(inputMesh:ConstantRectilinearMesh, destinationMesh:ConstantRectilinearMesh)-> coo_matrix:
    """return the transfer operator to transfer a field defined on the elements of the input mesh (inputMesh)
    to the output mesh.

    Parameters
    ----------
    inputMesh : ConstantRectilinearMesh
        the mesh defining the support of the element field to be transferred
    destinationMesh : ConstantRectilinearMesh
        destination mesh

    Returns
    -------
    np.ndarray
        OP: The transfer operation (sparse matrix in coo)
        field in the new mesh = OP.dot(field defined in the input Mesh)
    """
    if not isinstance(inputMesh, ConstantRectilinearMesh):
        raise Exception("First argument must be a ConstantRectilinearMesh")
    nps = 3
    nps3 = nps**inputMesh.GetDimensionality()
    oldToNewVals = np.zeros((destinationMesh.GetNumberOfNodes(),nps3))
    oldToNewIK = np.zeros((destinationMesh.GetNumberOfNodes(),nps3), dtype=np.int_)
    oldToNewJK = np.zeros((destinationMesh.GetNumberOfNodes(),nps3), dtype=np.int_)

    numberOfElementsDest = destinationMesh.GetNumberOfElements(dim = destinationMesh.GetElementsDimensionality() )
    numberOfElementsInp = inputMesh.GetNumberOfElements(dim = inputMesh.GetElementsDimensionality() )

    for i in  range(numberOfElementsDest):
        coon = destinationMesh.GetConnectivityForElement(i)

        n0pos = destinationMesh.GetPosOfNode(coon[0])
        cpt =0
        for cx in range(0,nps):
            for cy in range(0,nps):
                if inputMesh.GetDimensionality() == 3:
                    for cz in range(0,nps):
                        pos = n0pos + destinationMesh.GetSpacing()*([cx+0.5,cy+0.5,cz+0.5])/nps
                        el = inputMesh.GetElementAtPos(pos)
                        oldToNewVals[i,cpt] += 1./nps3
                        oldToNewIK[i,cpt] += i
                        oldToNewJK[i,cpt] += el
                        cpt +=1
                else:
                    pos = n0pos + destinationMesh.GetSpacing()*([cx+0.5,cy+0.5])/nps
                    el = inputMesh.GetElementAtPos(pos)
                    oldToNewVals[i,cpt] += 1./nps3
                    oldToNewIK[i,cpt] += i
                    oldToNewJK[i,cpt] += el
                    cpt +=1

    oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(numberOfElementsDest, numberOfElementsInp)).tocsc()

    return oldToNew


#------------------------------------------------------------------------------
def CreateSquare(dimensions:ArrayLike=[2,2], origin:ArrayLike=[-1.0,-1.0], spacing:ArrayLike=[1.,1.]) ->ConstantRectilinearMesh:
    """Create ConstantRectilinearMesh using the dimension, origin and spacing.
    Extra data is added:
        the bulk element with tag "2D"
        the border elements with tags "X0", "X1", "Y0", "Y1",
        the nodes of the corner with tags "x0y0", "x1y0", "x1y1", "x0y1"

    Parameters
    ----------
    dimensions : ArrayLike, optional
        Number of node in each dimension, by default [2,2]
    origin : ArrayLike, optional
        position of the first corner, by default [-1.0,-1.0]
    spacing : ArrayLike, optional
        size of the elements, by default [1.,1.]

    Returns
    -------
    ConstantRectilinearMesh
        the
    """
    spacing = np.asarray(spacing,dtype=PBasicFloatType)
    origin = np.asarray(origin,dtype=PBasicFloatType)
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshModificationTools import ComputeSkin
    import BasicTools.Containers.ElementNames as EN

    myMesh = ConstantRectilinearMesh(dim=2)
    myMesh.SetDimensions(dimensions)
    myMesh.SetOrigin(origin)
    myMesh.SetSpacing(spacing)

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


def CreateMesh(dim:int):
    """Helper class to create a ConstantRectilinearMesh of dimension dim and with only one element of size 1"""

    myMesh = ConstantRectilinearMesh(dim)
    myMesh.SetDimensions([2,]*dim)
    myMesh.SetSpacing([1, ]*dim)
    return myMesh

def CheckIntegrity_GetSubSuperMesh(dim):
    """CheckIntegrity function

    Parameters
    ----------
    dim : int
        dimensionality of the mesh to test
    """
    newmesh = GetSubSuperMesh(CreateMesh(dim),[3,]*dim)

def CheckIntegrity_GetNodeTransferMatrix(dim:int):
    """CheckIntegrity function

    Parameters
    ----------
    dim : int
        dimensionality of the mesh to test
    """

    mesh1 = CreateMesh(dim)
    mesh2 = GetSubSuperMesh(mesh1,[3,]*dim)

    TMatrix = GetNodeTransferMatrix(mesh1,mesh2)

def CheckIntegrity_GetElementTransferMatrix(dim:int):
    """CheckIntegrity function

    Parameters
    ----------
    dim : int
        dimensionality of the mesh to test
    """
    mesh1 = CreateMesh(dim)
    mesh2 = GetSubSuperMesh(mesh1,[3,]*dim)

    TMatrix = GetElementTransferMatrix(mesh1,mesh2)

def CheckIntegrity(GUI=False)-> str:
    """CheckIntegrity function. Tests

    Parameters
    ----------
    GUI : bool, optional
        if True, generate (in some case) an output on a new window, by default False

    Returns
    -------
    str
        ok if all ok
    """
    for dim in [2,3]:
        CheckIntegrity_GetSubSuperMesh(dim)
        CheckIntegrity_GetNodeTransferMatrix(dim)
        CheckIntegrity_GetElementTransferMatrix(dim)
    CreateSquare()
    return  "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(True))
    print("done")
