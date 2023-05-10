# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Union, Optional
import copy

import numpy as np

from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType, ArrayLike
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh, ElementsContainer
from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
from BasicTools.Containers.UnstructuredMeshModificationTools import ComputeSkin
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes

def CreateUniformMeshOfBars(startPoint: Union[PBasicFloatType,ArrayLike], stopPoint: Union[PBasicFloatType,ArrayLike], nbPoints: PBasicIndexType=50, secondOrder: bool=False) -> UnstructuredMesh:
    """Create a uniform mesh of bars.
    In the case the starPoint and stopPoint are scalar the nodes of the mesh are 3D (with zeros for the other columns)
    The user can give a list with only one component to generate nodes with only 1 coordinate
    0D elements are created at the start and stop points with tags "L", "H"
    nodal tags are created at the start and stop points with tags "L", "H"


    Parameters
    ----------
    startPoint : Union[PBasicFloatType,ArrayLike]
        Start point of the linear mesh
    stopPoint : Union[PBasicFloatType,ArrayLike]
        Stop point of the linear mesh
    nbPoints : PBasicIndexType, optional
        Number of point in the mesh, by default 50
    secondOrder : bool, optional
        if true second order bars are generated. In this case the number of points (nbPoints) must be odd

    Returns
    -------
    UnstructuredMesh
        the mesh in UnstructuredMesh format

    Raises
    ------
    Exception
        if the startPoint and stopPoint are not compatible
    """

    if not hasattr(startPoint,"__len__"):
        startPoint = [startPoint,0.,0.]

    if not hasattr(stopPoint,"__len__"):
        stopPoint = [stopPoint,0.,0.]

    if len(startPoint) != len(stopPoint):# pragma: no cover
        raise Exception(f"size of startPoint ({len(startPoint)}) and stopPoint ({len(stopPoint)}) not compatible")

    points = np.linspace(startPoint, stopPoint, nbPoints, dtype=PBasicFloatType)

    if secondOrder :
        if nbPoints % 2 == 0:# pragma: no cover
            raise(Exception("the number of point must be odd in secondOrder"))

        bars = np.empty(((nbPoints-1)//2 ,3))
        bars[:,0] = np.arange(0,nbPoints-2,2)
        bars[:,1] = np.arange(2,nbPoints,2)
        bars[:,2] = np.arange(1,nbPoints-1,2)
        res = CreateMeshOf(points,bars,elemName = ElementNames.Bar_3 )
        #print(bars)
        #raise
    else:
        bars = np.empty((nbPoints-1,2))
        bars[:,0] = np.arange(nbPoints-1)
        bars[:,1] = np.arange(1,nbPoints)
        res = CreateMeshOf(points,bars,elemName = ElementNames.Bar_2 )

    res.nodesTags.CreateTag("L",).SetIds([0])
    res.nodesTags.CreateTag("H",).SetIds([nbPoints-1])
    elements = res.GetElementsOfType(ElementNames.Point_1)
    elements.connectivity = np.array([[0],[nbPoints-1]],dtype=PBasicIndexType)
    elements.originalIds = np.arange(2,dtype=PBasicIndexType)
    elements.cpt = elements.connectivity.shape[0]
    elements.tags.CreateTag("L",).SetIds([0])
    elements.tags.CreateTag("H",).SetIds([1])
    res.PrepareForOutput()
    return res

def CreateMeshOfTriangles(points:ArrayLike, triangles:ArrayLike) -> UnstructuredMesh:
    """Helper function to create a Unstructured mesh using only points
    and the connectivity matrix for the triangles.

    Parameters
    ----------
    points : ArrayLike
        The positions of the point (nb Points,3)
    triangles : ArrayLike
        The connectivity of the triangles (nb Triangles, 3)

    Returns
    -------
    UnstructuredMesh
        A instance of UnstructuredMesh with the triangles
    """
    return CreateMeshOf(points, triangles, elemName = ElementNames.Triangle_3 )

def CreateMeshOf(points: ArrayLike, connectivity: ArrayLike, elemName: str) -> UnstructuredMesh:
    """Helper function to create a mesh of homogeneous elements.

    Parameters
    ----------
    points : ArrayLike
        The positions of the point (nb Points,3)
    connectivity : ArrayLike
        The connectivity of the elements (nb of elements, number of point per element)
    elemName : str, optional
        the name of the elements to create. , by default None

    Returns
    -------
    UnstructuredMesh
        A instance of UnstructuredMesh with the given point and elements

    """

    res = UnstructuredMesh()

    res.nodes = np.array(points, dtype=np.double)
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=PBasicIndexType)

    elements = res.GetElementsOfType(elemName)
    elements.connectivity = np.array(connectivity,dtype=PBasicIndexType)
    elements.originalIds = np.arange(0,elements.connectivity.shape[0],dtype=PBasicIndexType)
    elements.cpt = elements.connectivity.shape[0]
    elements.tags.CreateTag(str(ElementNames.dimension[elemName])+"D").SetIds(np.arange(elements.GetNumberOfElements() ) )
    res.PrepareForOutput()

    return res

def CreateSquare(dimensions: ArrayLike=[2,2], origin: ArrayLike=[-1.0,-1.0], spacing: ArrayLike=[1.,1.], ofTriangles: bool=False)-> UnstructuredMesh:
    """Create a mesh of a square

    Parameters
    ----------
    dimensions : ArrayLike, optional
        Number of points in every direction, by default [2,2]
    origin : ArrayLike, optional
        origin of the mesh , by default [-1.0,-1.0]
    spacing : ArrayLike, optional
        x and y size of every element, by default [1.,1.]
    ofTriangles : bool, optional
        , by default False

    Returns
    -------
    UnstructuredMesh
        _description_
    """
    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)

    myMesh = ConstantRectilinearMesh(dim=2)
    myMesh.SetDimensions(dimensions)
    myMesh.SetOrigin(origin)
    myMesh.SetSpacing(spacing)

    # corners
    d = np.array(dimensions)-1
    s = spacing
    indices = [[   0,   0,   0],
                [d[0],   0,   0],
                [   0,d[1],   0],
                [d[0],d[1],   0]]

    for n in indices:
        idx = myMesh.GetMonoIndexOfNode(n)
        name = "x"  + ("0" if n[0]== 0 else "1" )
        name += "y" + ("0" if n[1]== 0 else "1" )
        myMesh.nodesTags.CreateTag(name,False).SetIds([idx])


    mesh = CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex=ofTriangles)
    skin = ComputeSkin(mesh)
    for name,data in skin.elements.items():
        mesh.GetElementsOfType(name).Merge(data)
    #print(skin)

    if ofTriangles:
        triangles = mesh.GetElementsOfType(ElementNames.Triangle_3)
        triangles.GetTag("2D").SetIds(np.arange(triangles.GetNumberOfElements()))
    else:
        quads = mesh.GetElementsOfType(ElementNames.Quadrangle_4)
        quads.GetTag("2D").SetIds(np.arange(quads.GetNumberOfElements()))
    skin = mesh.GetElementsOfType(ElementNames.Bar_2)
    #face tags

    x = mesh.GetPosOfNodes()[skin.connectivity,0]
    y = mesh.GetPosOfNodes()[skin.connectivity,1]
    tol = np.min(spacing)/10

    skin.GetTag("X0").SetIds( np.where(np.sum(np.abs(x - origin[0]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("X1").SetIds( np.where(np.sum(np.abs(x - (origin[0]+d[0]*s[0]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Y0").SetIds( np.where(np.sum(np.abs(y - origin[1]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Y1").SetIds( np.where(np.sum(np.abs(y - (origin[1]+d[1]*s[1]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])


    mesh.PrepareForOutput()
    return mesh

def CreateDisk(nr: PBasicIndexType=10, nTheta: PBasicIndexType=10, r0: PBasicFloatType=0.5, r1: PBasicFloatType=1., theta0: PBasicFloatType=0, theta1: PBasicFloatType=np.pi/2 ) -> UnstructuredMesh:
    """Function to create a disk section

    Parameters
    ----------
    nr : PBasicIndexType, optional
        number of points in the radial direction, by default 10
    nTheta : PBasicIndexType, optional
        number of point in the angular direction, by default 10
    r0 : PBasicFloatType, optional
        internal radius, by default 0.5
    r1 : PBasicFloatType, optional
        external radius, by default 1.
    theta0 : PBasicFloatType, optional
        start angle, by default 0
    theta1 : PBasicFloatType, optional
        end angle, by default np.pi/2

    Returns
    -------
    UnstructuredMesh
        An UnstructuredMesh of a disk sector
    """

    myMesh = CreateSquare(dimensions=[nr,nTheta],origin=[r0,theta0],spacing=[(r1-r0)/(nr-1),(theta1-theta0)/(nTheta-1)])

    r = myMesh.nodes[:,0].copy()
    theta = myMesh.nodes[:,1].copy()

    myMesh.nodes[:,0] = r*np.cos(theta)
    myMesh.nodes[:,1] = r*np.sin(theta)

    myMesh.elements[ElementNames.Bar_2].tags["X0"].name = "R0"
    myMesh.elements[ElementNames.Bar_2].tags["X1"].name = "R1"
    myMesh.elements[ElementNames.Bar_2].tags["Y0"].name = "Theta0"
    myMesh.elements[ElementNames.Bar_2].tags["Y1"].name = "Theta1"

    return myMesh

def CreateMeshFromConstantRectilinearMesh(CRM: ConstantRectilinearMesh, ofSimplex: bool=False)-> UnstructuredMesh:
    """Function to convert a ConstantRectilinear mesh to a UnstructuredMesh

    Parameters
    ----------
    CRM : ConstantRectilinearMesh
        The input mesh
    ofSimplex : bool, optional
        True to make simplex (tetras and triangles), by default False

    Returns
    -------
    UnstructuredMesh
        An UnstructuredMesh instance of the input mesh
    """
    res = UnstructuredMesh()


    res.CopyProperties(CRM)

    res.nodes = np.copy(CRM.GetPosOfNodes())
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=PBasicIndexType)
    res.nodesTags = CRM.nodesTags.Copy()

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearElementContainer
    from BasicTools.Containers.UnstructuredMesh import ElementsContainer

    for name, data in CRM.elements.items():
        if isinstance(data,ConstantRectilinearElementContainer):
            elements = ElementsContainer(name)
            elements.connectivity = data.connectivity
            elements.tags = data.tags.Copy()
            elements.cpt = data.GetNumberOfElements()
            elements.originalIds = np.arange(0,data.GetNumberOfElements(),dtype=PBasicIndexType)
        else:
            elements = data.Copy()

        res.elements[name] = elements

    if ofSimplex:
        MeshToSimplex(res)
    return res

def CreateCube(dimensions: ArrayLike=[2,2,2], origin: ArrayLike=[-1.0,-1.0,-1.0], spacing: ArrayLike=[1.,1.,1.], ofTetras: bool=False)-> UnstructuredMesh:
    """Create a UnstructuredMesh of a cube

    Parameters
    ----------
    dimensions : ArrayLike, optional
        Number of point in every direction, by default [2,2,2]
    origin : ArrayLike, optional
        origin of the lowest point, by default [-1.0,-1.0,-1.0]
    spacing : ArrayLike, optional
        size in the 3 direction of every element, by default [1.,1.,1.]
    ofTetras : bool, optional
        if False: mesh of hexahedrons
        if True: mesh of tetrahedrons
        by default False

    Returns
    -------
    UnstructuredMesh
        A mesh of the cube
    """

    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions(dimensions)
    myMesh.SetOrigin(origin)
    myMesh.SetSpacing(spacing)

    # corners
    d = np.array(dimensions)-1
    s = spacing
    indices = [[   0,   0,   0],
                [d[0],   0,   0],
                [   0,d[1],   0],
                [d[0],d[1],   0],
                [   0,   0,d[2]],
                [d[0],   0,d[2]],
                [   0,d[1],d[2]],
                [d[0],d[1],d[2]]]

    for n in indices:
        idx = myMesh.GetMonoIndexOfNode(n)
        name = "x"  + ("0" if n[0]== 0 else "1" )
        name += "y" + ("0" if n[1]== 0 else "1" )
        name += "z" + ("0" if n[2]== 0 else "1" )
        myMesh.nodesTags.CreateTag(name,False).SetIds([idx])

    mesh = CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex=ofTetras)
    skin = ComputeSkin(mesh)
    for name,data in skin.elements.items():
        mesh.GetElementsOfType(name).Merge(data)

    if ofTetras:
        tets = mesh.GetElementsOfType(ElementNames.Tetrahedron_4)
        tets.GetTag("3D").SetIds(np.arange(tets.GetNumberOfElements()))
        skin = mesh.GetElementsOfType(ElementNames.Triangle_3)
    else:
        hexahedrons = mesh.GetElementsOfType(ElementNames.Hexaedron_8)
        hexahedrons.GetTag("3D").SetIds(np.arange(hexahedrons.GetNumberOfElements()))
        skin = mesh.GetElementsOfType(ElementNames.Quadrangle_4)
    #face tags

    x = mesh.GetPosOfNodes()[skin.connectivity,0]
    y = mesh.GetPosOfNodes()[skin.connectivity,1]
    z = mesh.GetPosOfNodes()[skin.connectivity,2]
    tol = np.min(spacing)/10

    skin.GetTag("X0").SetIds( np.where(np.sum(np.abs(x - origin[0]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("X1").SetIds( np.where(np.sum(np.abs(x - (origin[0]+d[0]*s[0]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Y0").SetIds( np.where(np.sum(np.abs(y - origin[1]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Y1").SetIds( np.where(np.sum(np.abs(y - (origin[1]+d[1]*s[1]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Z0").SetIds( np.where(np.sum(np.abs(z - origin[2]          )<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])
    skin.GetTag("Z1").SetIds( np.where(np.sum(np.abs(z - (origin[2]+d[2]*s[2]))<tol,axis=1) == skin.GetNumberOfNodesPerElement())[0])

    mesh.PrepareForOutput()
    return mesh

def MeshToSimplex(mesh: UnstructuredMesh) -> UnstructuredMesh:
    """(EXPERIMENTAL) Convert mesh to only tetrahedron/triangle/bars/points

    Warning!!! This function is dangerous, we don't check the compatibility at each interface between
    the elements.

    Parameters
    ----------
    mesh : UnstructuredMesh
        A mesh

    Returns
    -------
    UnstructuredMesh
        A mesh composed only by simplex
    """


    from BasicTools.Containers.UnstructuredMesh import ElementsContainer,AllElements

    ae = AllElements()


    for elementName, data in mesh.elements.items():
        res = data

        if elementName == ElementNames.Hexaedron_8:

            res = ElementsContainer(ElementNames.Tetrahedron_4)
            nbElements = data.GetNumberOfElements()
            res.Allocate(nbElements*6)
            conn = data.connectivity
            res.connectivity[0:nbElements*6:6] = conn[:,[0, 6, 2, 3]]
            res.connectivity[1:nbElements*6:6] = conn[:,[0, 6, 3, 7]]
            res.connectivity[2:nbElements*6:6] = conn[:,[0, 6, 7, 4]]
            res.connectivity[3:nbElements*6:6] = conn[:,[0, 6, 4, 5]]
            res.connectivity[4:nbElements*6:6] = conn[:,[0, 6, 5, 1]]
            res.connectivity[5:nbElements*6:6] = conn[:,[0, 6, 1, 2]]

            res.originalIds =  np.repeat(data.originalIds,6)
            for tagName in data.tags.keys():
                ids = data.tags[tagName].GetIds()
                res.tags.CreateTag(tagName).SetIds(np.repeat(ids,6)*6+np.tile(range(6),len(ids)) )

        elif elementName == ElementNames.Quadrangle_4:

            res = ElementsContainer(ElementNames.Triangle_3)
            nbElements = data.GetNumberOfElements()
            res.Allocate(nbElements*2)
            conn = data.connectivity
            res.connectivity[0:nbElements*2:2] = conn[:,[0, 1, 2]]
            res.connectivity[1:nbElements*2:2] = conn[:,[0, 2, 3]]

            res.originalIds =  np.repeat(data.originalIds,2)
            for tagName in data.tags:
                ids = data.tags[tagName].GetIds()
                res.tags.CreateTag(tagName).SetIds(np.repeat(ids,2)*2+np.tile(range(2),len(ids)) )
        elif elementName in [ElementNames.Triangle_3,ElementNames.Triangle_6,ElementNames.Tetrahedron_4,ElementNames.Tetrahedron_10,ElementNames.Bar_2,ElementNames.Bar_3,ElementNames.Point_1]  :
            pass
        else:
            raise(Exception("Don't know how to convert {} to simplices".format(elementName)))

        if elementName in ae:
            ae[res.elementType].merge(res)
        else:
            ae[res.elementType] = res

    mesh.elements = ae



def ToQuadraticMesh(inputMesh: UnstructuredMesh) -> UnstructuredMesh:
    """Function to convert any mesh to a quadratic mesh.
    Nodes fields and element fields are lost.

    Parameters
    ----------
    inputMesh : UnstructuredMesh
        the in put mesh

    Returns
    -------
    UnstructuredMesh
        A mesh composed only by quadratic elements
    """

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP2
    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.Containers.Filters import ElementFilter
    from BasicTools.FE.Fields.FieldTools import FillFEField

    dim = inputMesh.GetDimensionality()

    numbering = ComputeDofNumbering(inputMesh, LagrangeSpaceP2)

    res = UnstructuredMesh()
    res.nodes = np.empty((numbering.size,dim), dtype=PBasicFloatType)

    def GetPos(j):
        return lambda x: x[j]

    for i in range(dim):
        newPos_ =  FEField("newPos_", inputMesh, LagrangeSpaceP2, numbering)
        newPos_.data = res.nodes[:,i]
        FillFEField(newPos_,[(ElementFilter(inputMesh),GetPos(i))])

    res.originalIDNodes = np.zeros(numbering.size, dtype=PBasicIndexType)-1
    newIds = np.empty(inputMesh.GetNumberOfNodes(), dtype=PBasicIndexType)
    for i in range(inputMesh.GetNumberOfNodes()):
        dof = numbering.GetDofOfPoint(i)
        newIds[i] = dof

    res.originalIDNodes[newIds] = range(inputMesh.GetNumberOfNodes())

    for tag in inputMesh.nodesTags:
        res.nodesTags.CreateTag(tag.name).SetIds(newIds[tag.GetIds()])

    import BasicTools.Containers.ElementNames as EN
    toQuad = {}
    toQuad[ElementNames.Point_1] = ElementNames.Point_1
    toQuad[ElementNames.Bar_2] = ElementNames.Bar_3
    toQuad[ElementNames.Bar_3] = ElementNames.Bar_3
    toQuad[ElementNames.Triangle_3] = ElementNames.Triangle_6
    toQuad[ElementNames.Triangle_6] = ElementNames.Triangle_6
    toQuad[ElementNames.Tetrahedron_4]  = ElementNames.Tetrahedron_10
    toQuad[ElementNames.Tetrahedron_10] = ElementNames.Tetrahedron_10
    toQuad[ElementNames.Quadrangle_4] = ElementNames.Quadrangle_9
    toQuad[ElementNames.Quadrangle_8] = ElementNames.Quadrangle_9
    toQuad[ElementNames.Quadrangle_9] = ElementNames.Quadrangle_9
    toQuad[ElementNames.Hexaedron_8]  = ElementNames.Hexaedron_27
    toQuad[ElementNames.Hexaedron_20] = ElementNames.Hexaedron_27
    toQuad[ElementNames.Hexaedron_27] = ElementNames.Hexaedron_27
    toQuad[ElementNames.Wedge_6] = ElementNames.Wedge_15
    toQuad[ElementNames.Wedge_15] = ElementNames.Wedge_15

    for name, data in inputMesh.elements.items():
        if data.GetNumberOfElements() == 0:
            continue
        elements = res.elements.GetElementsOfType(toQuad[name])
        elements.AddNewElements(numbering[name],data.originalIds)
        for tag in data.tags:
            elements.tags.CreateTag(tag.name,False).SetIds(tag.GetIds())

    res.PrepareForOutput()

    return res

def QuadToLin(inputMesh: UnstructuredMesh, divideQuadElements: bool=True, linearizedMiddlePoints: bool=False)-> UnstructuredMesh:
    """Convert a quadratic mesh to a linear mesh

    Parameters
    ----------
    inputMesh : UnstructuredMesh
        the input mesh
    divideQuadElements : bool, optional
        if the quadratic element must be divided, by default True
    linearizedMiddlePoints : bool, optional
        if the middle point of the quadratic element must line in the segment forme by the extreme values , by default False

    Returns
    -------
    UnstructuredMesh
        a linear mesh

    """

    from BasicTools.Containers.UnstructuredMeshFieldOperations import QuadFieldToLinField

    res = type(inputMesh)()
    res.CopyProperties(inputMesh)

    res.nodes = np.copy(inputMesh.GetPosOfNodes())
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=PBasicIndexType)

    res.nodesTags = inputMesh.nodesTags.Copy()

    for elementName in inputMesh.elements:
        quadElement = inputMesh.elements[elementName]
        if elementName == ElementNames.Tetrahedron_10:

            linearElements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
            initNbElem = linearElements.GetNumberOfElements()
            if divideQuadElements:
                nbOfNewElements = 8
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*8)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*8
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i*8+0,:] = quadConn[[0,4,6,7]]
                    linearElements.connectivity[initNbElem+i*8+1,:] = quadConn[[1,5,4,8]]
                    linearElements.connectivity[initNbElem+i*8+2,:] = quadConn[[2,6,5,9]]
                    linearElements.connectivity[initNbElem+i*8+3,:] = quadConn[[7,8,9,3]]
                    linearElements.connectivity[initNbElem+i*8+4,:] = quadConn[[4,5,6,7]]
                    linearElements.connectivity[initNbElem+i*8+5,:] = quadConn[[4,5,7,8]]
                    linearElements.connectivity[initNbElem+i*8+6,:] = quadConn[[5,6,7,9]]
                    linearElements.connectivity[initNbElem+i*8+7,:] = quadConn[[5,7,8,9]]
                    if linearizedMiddlePoints :
                        res.nodes[quadConn[4],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                        res.nodes[quadConn[5],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[6],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[0],:] )/2
                        res.nodes[quadConn[7],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[3],:] )/2
                        res.nodes[quadConn[8],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[3],:] )/2
                        res.nodes[quadConn[9],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[3],:] )/2
            else:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*1)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*1
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3]]


        elif elementName == ElementNames.Triangle_6:

            linearElements = res.GetElementsOfType(ElementNames.Triangle_3)
            initNbElem = linearElements.GetNumberOfElements()
            if divideQuadElements:
                nbOfNewElements = 4
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*4)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*4
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i*4+0,:] = quadConn[[0,3,5]]
                    linearElements.connectivity[initNbElem+i*4+1,:] = quadConn[[1,4,3]]
                    linearElements.connectivity[initNbElem+i*4+2,:] = quadConn[[2,5,4]]
                    linearElements.connectivity[initNbElem+i*4+3,:] = quadConn[[3,4,5]]
                    if linearizedMiddlePoints :
                        res.nodes[quadConn[3],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                        res.nodes[quadConn[4],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[5],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[0],:] )/2
            else:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i,:] = quadConn[[0,1,2]]

        elif elementName == ElementNames.Quadrangle_8:

            linearElements = res.GetElementsOfType(ElementNames.Quadrangle_4)
            initNbElem = linearElements.GetNumberOfElements()
            if divideQuadElements:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*1)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*1
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i*1+0,:] = quadConn[[0,1,2,3]]
                    #if linearizedMiddlePoints :
                    #    res.nodes[quadConn[4],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                    #    res.nodes[quadConn[5],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                    #    res.nodes[quadConn[6],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[3],:] )/2
                    #    res.nodes[quadConn[7],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[0],:] )/2
            else:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3]]

        elif elementName == ElementNames.Hexaedron_20:

            linearElements = res.GetElementsOfType(ElementNames.Hexaedron_8)
            initNbElem = linearElements.GetNumberOfElements()
            if divideQuadElements:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*1)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*1
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i*1+0,:] = quadConn[[0,1,2,3,4,5,6,7]]
                    if linearizedMiddlePoints :
                        res.nodes[quadConn[8],:]  = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                        res.nodes[quadConn[9],:]  = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[10],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[3],:] )/2
                        res.nodes[quadConn[11],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[0],:] )/2
                        res.nodes[quadConn[12],:] = (res.nodes[quadConn[4],:] + res.nodes[quadConn[5],:] )/2
                        res.nodes[quadConn[13],:] = (res.nodes[quadConn[5],:] + res.nodes[quadConn[6],:] )/2
                        res.nodes[quadConn[14],:] = (res.nodes[quadConn[6],:] + res.nodes[quadConn[7],:] )/2
                        res.nodes[quadConn[15],:] = (res.nodes[quadConn[7],:] + res.nodes[quadConn[4],:] )/2
                        res.nodes[quadConn[16],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[4],:] )/2
                        res.nodes[quadConn[17],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[5],:] )/2
                        res.nodes[quadConn[18],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[6],:] )/2
                        res.nodes[quadConn[19],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[7],:] )/2
            else:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3,4,5,6,7]]


        elif elementName == ElementNames.Hexaedron_27:

            linearElements = res.GetElementsOfType(ElementNames.Hexaedron_8)
            initNbElem = linearElements.GetNumberOfElements()
            if divideQuadElements:
                nbOfNewElements = 8
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*8)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*8
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i*8+0,:] = quadConn[[0,8,24,11,16,22,26,20]]
                    linearElements.connectivity[initNbElem+i*8+1,:] = quadConn[[8,0,9,24,22,17,21,26]]
                    linearElements.connectivity[initNbElem+i*8+2,:] = quadConn[[11,24,10,3,20,26,23,19]]
                    linearElements.connectivity[initNbElem+i*8+3,:] = quadConn[[24,9,2,10,26,21,18,23]]
                    linearElements.connectivity[initNbElem+i*8+4,:] = quadConn[[16,22,26,20,4,12,25,15]]
                    linearElements.connectivity[initNbElem+i*8+5,:] = quadConn[[22,17,21,26,12,5,13,25]]
                    linearElements.connectivity[initNbElem+i*8+6,:] = quadConn[[20,26,23,19,15,25,14,7]]
                    linearElements.connectivity[initNbElem+i*8+7,:] = quadConn[[26,21,18,23,25,13,6,14]]

                    if linearizedMiddlePoints :
                        res.nodes[quadConn[8],:]  = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                        res.nodes[quadConn[9],:]  = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[10],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[3],:] )/2
                        res.nodes[quadConn[11],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[0],:] )/2
                        res.nodes[quadConn[12],:] = (res.nodes[quadConn[4],:] + res.nodes[quadConn[5],:] )/2
                        res.nodes[quadConn[13],:] = (res.nodes[quadConn[5],:] + res.nodes[quadConn[6],:] )/2
                        res.nodes[quadConn[14],:] = (res.nodes[quadConn[6],:] + res.nodes[quadConn[7],:] )/2
                        res.nodes[quadConn[15],:] = (res.nodes[quadConn[7],:] + res.nodes[quadConn[4],:] )/2
                        res.nodes[quadConn[16],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[4],:] )/2
                        res.nodes[quadConn[17],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[5],:] )/2
                        res.nodes[quadConn[18],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[6],:] )/2
                        res.nodes[quadConn[19],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[7],:] )/2
                        res.nodes[quadConn[20],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[4],:] )/2
                        res.nodes[quadConn[21],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[6],:] )/2
                        res.nodes[quadConn[22],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[5],:] )/2
                        res.nodes[quadConn[23],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[7],:] )/2
                        res.nodes[quadConn[24],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[25],:] = (res.nodes[quadConn[4],:] + res.nodes[quadConn[6],:] )/2
                        res.nodes[quadConn[26],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[6],:] )/2


            else:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3,4,5,6,7]]


        elif elementName == ElementNames.Bar_3:

            linearElements = res.GetElementsOfType(ElementNames.Bar_2)
            initNbElem = linearElements.GetNumberOfElements()
            if divideQuadElements:
                nbOfNewElements = 2
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements()*2)
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()*2
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i*2+0,:] = quadConn[[0,2]]
                    linearElements.connectivity[initNbElem+i*2+1,:] = quadConn[[2,1]]
                    if linearizedMiddlePoints :
                        res.nodes[quadConn[2],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
            else:
                nbOfNewElements = 1
                linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:]
                    linearElements.connectivity[initNbElem+i,:] = quadConn[[0,1]]
        elif ElementNames.linear[elementName] :
            linearElements = res.GetElementsOfType(elementName)
            initNbElem = linearElements.GetNumberOfElements()

            linearElements.Reserve(initNbElem+quadElement.GetNumberOfElements())
            nbOfNewElements = 1
            linearElements.connectivity[initNbElem:initNbElem+quadElement.GetNumberOfElements(),:] = quadElement.connectivity
            linearElements.cpt = initNbElem+quadElement.GetNumberOfElements()

        else:
            raise Exception('Error : not coded yet for this type of elements ' + str(elementName))# pragma: no cover
        #copy of tags
        for originalTag in quadElement.tags :
            destinationTag = linearElements.GetTag(originalTag.name)
            ids = originalTag.GetIds()
            for i in range(originalTag.cpt):
                for t in range(nbOfNewElements):
                    destinationTag.AddToTag(initNbElem+ids[i]*nbOfNewElements+t)

            destinationTag.Tighten()

        res.ComputeGlobalOffset()

    res.PrepareForOutput()

    if divideQuadElements == False:
        CleanLonelyNodes(res)

    return res

def MirrorMesh(inmesh :UnstructuredMesh, x: Optional[PBasicFloatType]=None, y: Optional[PBasicFloatType]=None, z: Optional[PBasicFloatType]=None) -> UnstructuredMesh:
    """Compute a new symmetric mesh from the inmesh
    All inMesh element with be copied to the output
    The user is responsible to call RemoveDoubles to eliminate the possible double nodes on the symmetry plane

    Parameters
    ----------
    inmesh : UnstructuredMesh
        input mesh
    x : Optional[PBasicFloatType], optional
        the position of the yz plane to build the symmetric part, if Non, no symmetry to the plane yz
        , by default None
    y : Optional[PBasicFloatType], optional
        the position of the xz plane to build the symmetric part, if Non, no symmetry to the plane xz
        , by default None
    z : Optional[PBasicFloatType], optional
        the position of the xy plane to build the symmetric part, if Non, no symmetry to the plane xy
        _description_, by default None

    Returns
    -------
    UnstructuredMesh
        the output mesh
    """
    nbPoints = inmesh.GetNumberOfNodes()

    outMesh = type(inmesh)()
    outMesh.CopyProperties(inmesh)

    d = 0
    if x is not None:
        d += 1
    if y is not None:
        d += 1
    if z is not None:
        d += 1

    outMesh.nodes = np.empty((nbPoints*(2**d),inmesh.GetDimensionality()), dtype=PBasicFloatType)
    outMesh.originalIDNodes = np.empty((nbPoints*(2**d),), dtype=PBasicIndexType)

    #copy of points:
    outMesh.nodes[0:nbPoints,:] = inmesh.nodes

    #copy of points:
    outMesh.originalIDNodes[0:nbPoints] = inmesh.originalIDNodes


    outMesh.nodesTags = inmesh.nodesTags.Copy()
    cpt = nbPoints

    def increaseTags(tags,oldSize):
        for tag in tags:
            ids = tag.GetIds()[:]  # make a copy
            tag.SetIds(np.hstack((ids,ids+oldSize)) )

    if x is not None:
        vec = np.array([ [  -1,1,1], ],dtype=float)
        outMesh.nodes[cpt:(2*cpt),:] = (outMesh.nodes[0:cpt,:]-[x,0,0])*vec+ [x,0,0]
        outMesh.originalIDNodes[cpt:(2*cpt)] = outMesh.originalIDNodes[0:cpt]
        increaseTags(outMesh.nodesTags,cpt)
        cpt = cpt*2

    if y is not None:
        vec = np.array([ [  1,-1,1], ],dtype=float)
        outMesh.nodes[cpt:(2*cpt),:] = (outMesh.nodes[0:cpt,:]-[0,y,0])*vec+ [0,y,0]
        outMesh.originalIDNodes[cpt:(2*cpt)] = outMesh.originalIDNodes[0:cpt]
        increaseTags(outMesh.nodesTags,cpt)
        cpt = cpt*2

    if z is not None:
        vec = np.array([ [  1,1,-1], ],dtype=float)
        outMesh.nodes[cpt:(2*cpt),:] = (outMesh.nodes[0:cpt,:]-[0,0,z])*vec+ [0,0,z]
        outMesh.originalIDNodes[cpt:(2*cpt)] = outMesh.originalIDNodes[0:cpt]
        increaseTags(outMesh.nodesTags,cpt)
        cpt = cpt*2

    for name,vals in inmesh.elements.items():
        nbElements = vals.GetNumberOfElements()
        outElements = outMesh.GetElementsOfType(name)
        outElements.Reserve(nbElements*(2**d))
        outElements.connectivity[0:nbElements,:] = vals.connectivity
        outElements.tags = copy.deepcopy(vals.tags)
        cpt = nbElements
        pointCpt = nbPoints
        permutation = ElementNames.mirrorPermutation[name]
        if x is not None:
            outElements.connectivity[cpt:(2*cpt),:] = (outElements.connectivity[0:cpt,:]+pointCpt)[:,permutation]
            pointCpt = pointCpt*2
            increaseTags(outElements.tags,cpt)
            cpt = cpt*2

        if y is not None:
            outElements.connectivity[cpt:(2*cpt),:] = (outElements.connectivity[0:cpt,:]+pointCpt)[:,permutation]
            pointCpt = pointCpt*2
            increaseTags(outElements.tags,cpt)
            cpt = cpt*2

        if z is not None:
            outElements.connectivity[cpt:(2*cpt),:] = (outElements.connectivity[0:cpt,:]+pointCpt)[:,permutation]
            pointCpt = pointCpt*2
            increaseTags(outElements.tags,cpt)
            cpt = cpt*2

        outElements.cpt = cpt

    outMesh.PrepareForOutput()
    return outMesh

def Create0DElementContainerForEveryPoint(mesh: UnstructuredMesh) -> ElementsContainer:
    """Create a ElementContainer with 0D elements based on the point of the mesh
    The user is responsible to put this ElementContainer back to a mesh (using merge for example)
    the Nodal tags are transferred to the new 0D elements
    the originalIds is constructed based on the nodes (range(nb Nodes))

    Parameters
    ----------
    mesh : UnstructuredMesh
        the input mesh

    Returns
    -------
    ElementsContainer
        An ElementsContainer with the new 0D elements
    """
    nbNodes = mesh.GetNumberOfNodes()
    elements1D = ElementsContainer(ElementNames.Point_1)
    elements1D.tags = mesh.nodesTags.Copy()
    elements1D.connectivity = np.arange(nbNodes)
    elements1D.connectivity.shape  = (nbNodes,1)
    elements1D.originalIds = np.arange(nbNodes)
    elements1D.cpt = mesh.GetNumberOfNodes()
    return elements1D

def SubDivideMesh(mesh: UnstructuredMesh, level: PBasicIndexType=1) -> UnstructuredMesh :
    """Subdivide a mesh

    Parameters
    ----------
    mesh : UnstructuredMesh
        The input mesh
    level : PBasicIndexType, optional
        Number of times the mesh must be divided, by default 1

    Returns
    -------
    UnstructuredMesh
        A new mesh divided level times
    """
    if level == 0:
        return mesh

    res = UnstructuredMesh()

    subdivisionAlmanac = {}
    subdivisionAlmanac[ElementNames.Point_1] = [(ElementNames.Point_1, [0]) ]
    subdivisionAlmanac[ElementNames.Bar_2] = [(ElementNames.Bar_2, [0,2]),
                                                (ElementNames.Bar_2, [2,1])]

    subdivisionAlmanac[ElementNames.Quadrangle_4] = [(ElementNames.Quadrangle_4, [0,4,8,7]),
                                                    (ElementNames.Quadrangle_4, [4,1,5,8]),
                                                    (ElementNames.Quadrangle_4, [8,5,2,6]),
                                                    (ElementNames.Quadrangle_4, [7,8,6,3])]

    subdivisionAlmanac[ElementNames.Hexaedron_8] = [(ElementNames.Hexaedron_8, [ 0, 8,24,11,16,22,26,20]),
                                                    (ElementNames.Hexaedron_8, [ 8, 1, 9,24,22,17,21,26]),
                                                    (ElementNames.Hexaedron_8, [24, 9, 2,10,26,21,18,23]),
                                                    (ElementNames.Hexaedron_8, [11,24,10, 3,20,26,23,19]),
                                                    (ElementNames.Hexaedron_8, [16,22,26,20, 4,12,25,15]),
                                                    (ElementNames.Hexaedron_8, [22,17,21,26,12, 5,13,25]),
                                                    (ElementNames.Hexaedron_8, [26,21,18,23,25,13, 6,14]),
                                                    (ElementNames.Hexaedron_8, [20,26,23,19,15,25,14, 7])]

    subdivisionAlmanac[ElementNames.Triangle_3] =  [(ElementNames.Triangle_3, [ 0, 3, 5]),
                                                    (ElementNames.Triangle_3, [ 3, 1, 4]),
                                                    (ElementNames.Triangle_3, [ 3, 4, 5]),
                                                    (ElementNames.Triangle_3, [ 5, 4, 2])]

    subdivisionAlmanac[ElementNames.Tetrahedron_4] =  [(ElementNames.Tetrahedron_4, [ 0, 4, 6, 7]),
                                                    (ElementNames.Tetrahedron_4, [ 1, 5, 4, 8]),
                                                    (ElementNames.Tetrahedron_4, [ 2, 6, 5, 9]),
                                                    (ElementNames.Tetrahedron_4, [ 7, 8, 9, 3]),

                                                    (ElementNames.Tetrahedron_4, [ 5, 6, 7, 9]),
                                                    (ElementNames.Tetrahedron_4, [ 5, 6, 4, 7]),
                                                    (ElementNames.Tetrahedron_4, [ 5, 7, 4, 8]),
                                                    (ElementNames.Tetrahedron_4, [ 5, 9, 7, 8]),

                                                    ]


    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo, LagrangeSpaceP2
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.IntegrationsRules import NodalEvaluationP2

    numberingGeo = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)
    numberingP2 = ComputeDofNumbering(mesh,LagrangeSpaceP2)

    ## Generation of nodes
    res.nodes = np.empty((numberingP2.size,3), dtype=PBasicFloatType)
    res.originalIDNodes = np.zeros(res.nodes.shape[0],dtype=PBasicIndexType)-1

    oldToNewDofs =  np.zeros(mesh.GetNumberOfNodes(),dtype=PBasicIndexType)

    oldToNewDofs[numberingP2["doftopointLeft"]] = numberingP2["doftopointRight"]

    res.originalIDNodes[oldToNewDofs] = mesh.originalIDNodes

    for tag in mesh.nodesTags.keys():
        name = mesh.nodesTags[tag].name
        ids = mesh.nodesTags[tag].GetIds()
        res.nodesTags.CreateTag(name).SetIds(oldToNewDofs[ids])

    for elemType, data in mesh.elements.items():
        spaceGeo = LagrangeSpaceGeo[elemType]
        spaceGeo.Create()
        p,w = NodalEvaluationP2[elemType]
        sGeoAtIp = spaceGeo.SetIntegrationRule(p,w)

        nGeo = numberingGeo[elemType]
        nP2 = numberingP2[elemType]

        # generation of nodes
        for sf in range(len(p)):
            geoNs = sGeoAtIp.valN[sf]
            for c in range(3):
                res.nodes[nP2[:,sf], c] = np.sum(mesh.nodes[:,c][nGeo]*geoNs,axis=1)

        #generation of elements
        for t, nn in subdivisionAlmanac[elemType]:
            newElements = res.GetElementsOfType(t)
            offset = newElements.GetNumberOfElements()
            tne =  newElements.AddNewElements(nP2[:,nn], data.originalIds)
            for tag in data.tags.keys():
                name = data.tags[tag].name
                ids = data.tags[tag].GetIds()
                if len(ids) == 0:
                    continue
                newElements.tags.CreateTag(name,False).AddToTag(ids+offset)

    mesh.PrepareForOutput()
    return SubDivideMesh(res,level-1)
#------------------------- CheckIntegrity ------------------------
def CheckIntegrity_CreateDisk(GUI=False):
    """ CheckIntegrity_CreateDisk """
    a = CreateDisk()
    return "ok"

def CheckIntegrity_SubDivideMesh(GUI=False):
    points = [[0,0,0],
                [1,0,0],
                [1,1,0],
                [0,1,0],
                [0,0,1],
                [1,0,1.5],
                [1,1,1],
                [0,1,1.5]]
    hexahedrons = [[0,1,2,3,4,5,6,7],]
    mesh = CreateMeshOf(points, hexahedrons, ElementNames.Hexaedron_8)
    mesh.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    mesh.GetElementsOfType(ElementNames.Hexaedron_8)
    mesh.GetElementsOfType(ElementNames.Hexaedron_8).tags
    mesh.GetElementsOfType(ElementNames.Hexaedron_8).tags.CreateTag("OnlyHex")
    mesh.GetElementsOfType(ElementNames.Hexaedron_8).tags.CreateTag("OnlyHex",False).AddToTag(0)

    outMesh = SubDivideMesh(mesh,1)

    print(mesh)
    print(outMesh)
    if GUI:
        from BasicTools.Bridges.vtkBridge import PlotMesh
        PlotMesh(outMesh)

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    WriteMeshToXdmf(tempdir+"CheckIntegrity_SubDivideMesh.xdmf", outMesh, PointFields=[outMesh.originalIDNodes],PointFieldsNames=["originalIDNodes"] )
    print(tempdir)

    if outMesh.GetNumberOfNodes() != 27:
        raise# pragma: no cover

    if outMesh.GetNumberOfElements() != 8:
        raise# pragma: no cover
    return "ok"


def CheckIntegrity_Create0DElementContainerForEveryPoint(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("FirstPoint").AddToTag(0)

    print(mesh)
    zeroDElements = Create0DElementContainerForEveryPoint(mesh)
    print(zeroDElements)
    if zeroDElements.GetNumberOfElements() != mesh.GetNumberOfNodes():
        raise(Exception("Error "))
    return "ok"

def CheckIntegrity_MirrorMesh(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    mesh.GetElementsOfType(ElementNames.Tetrahedron_4).tags.CreateTag("OnlyTet").AddToTag(0)
    outMesh = MirrorMesh(mesh,x=0)

    if outMesh.GetNumberOfNodes() != 8:
        raise# pragma: no cover

    if outMesh.GetNumberOfElements() != 2:
        raise# pragma: no cover
    return "ok"


def CheckIntegrity_QuadToLin(GUI=False):
    myMesh = UnstructuredMesh()
    myMesh.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0.5,0,0],[0.5,0.5,0],[0,0.5,0],[0,0,0.5],[0.5,0,0.5],[0,0.5,0.5]] ,dtype=float)
    tag = myMesh.GetNodalTag("linPoints")
    tag.AddToTag(0)
    tag.AddToTag(1)
    tag.AddToTag(2)
    tag.AddToTag(3)
    import BasicTools.Containers.ElementNames as ElementNames

    elements = myMesh.GetElementsOfType(ElementNames.Tetrahedron_10)
    elements.AddNewElement([0,1,2,3,4,5,6,7,8,9],0)
    elements = myMesh.GetElementsOfType(ElementNames.Triangle_6)
    elements.AddNewElement([0,1,2,4,5,6],1)
    elements = myMesh.GetElementsOfType(ElementNames.Bar_3)
    elements.AddNewElement([0,1,4],2)
    elements = myMesh.GetElementsOfType(ElementNames.Bar_2)
    elements.AddNewElement([0,1],3)

    myMesh.AddElementToTagUsingOriginalId(3,'LinElements')

    print(myMesh)
    linMesh = QuadToLin(myMesh,divideQuadElements=False)
    print(linMesh)
    print(QuadToLin(myMesh,divideQuadElements=True))
    print(QuadToLin(myMesh,divideQuadElements=True,linearizedMiddlePoints=True))
    from BasicTools.Containers.UnstructuredMeshFieldOperations import QuadFieldToLinField
    QuadFieldToLinField(myMesh, np.arange(myMesh.GetNumberOfNodes()))
    QuadFieldToLinField(myMesh, np.arange(myMesh.GetNumberOfNodes()), linMesh)
    return "ok"


def CheckIntegrity_CreateMeshFromConstantRectilinearMesh(GUI=False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(2)
    myMesh.SetDimensions([3,3])
    myMesh.SetSpacing([1, 1])

    print(myMesh)
    print(CreateMeshFromConstantRectilinearMesh(myMesh))

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([3,3,3])
    myMesh.SetSpacing([1, 1,1])
    for name,data in myMesh.elements.items():
        data.tags.CreateTag("First Element per Type").SetIds([0])
    print(myMesh)
    print(CreateMeshFromConstantRectilinearMesh(myMesh))

    res2 = CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex=True)
    print(res2.GetNumberOfElements())

    return "OK"

def CheckIntegrity_CreateMeshOfTriangles(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    print(res)
    return "OK"


def CheckIntegrity_CreateCube(GUI = False):
    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=False)
    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)
    return "ok"


def CheckIntegrity_CreateSquare(GUI = False):
    mesh = CreateSquare(dimensions=[20,21],spacing=[2.,2.],ofTriangles=False)
    mesh = CreateSquare(dimensions=[20,21],spacing=[2.,2.],ofTriangles=True)
    return "ok"

def CheckIntegrity_CreateUniformMeshOfBars(GUI=False):
    print(CreateUniformMeshOfBars(0,8,10))
    print(CreateUniformMeshOfBars(0,8,11,secondOrder=True))
    return "ok"

def CheckIntegrity_ToQuadraticMesh(GUI=False):
    inMesh =CreateCube()
    outMesh = ToQuadraticMesh(inMesh)

    return "ok"

def CheckIntegrity(GUI=False):
    toTest= [
    CheckIntegrity_ToQuadraticMesh,
    CheckIntegrity_CreateUniformMeshOfBars,
    CheckIntegrity_CreateCube,
    CheckIntegrity_CreateSquare,
    CheckIntegrity_CreateMeshOfTriangles,
    CheckIntegrity_CreateMeshFromConstantRectilinearMesh,
    CheckIntegrity_QuadToLin,
    CheckIntegrity_MirrorMesh,
    CheckIntegrity_Create0DElementContainerForEveryPoint,
    CheckIntegrity_SubDivideMesh,
    CheckIntegrity_CreateDisk,
    ]
    for f in toTest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
