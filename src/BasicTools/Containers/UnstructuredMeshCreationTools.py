# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
from BasicTools.Containers.UnstructuredMeshModificationTools import ComputeSkin
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes

def CreateUniformMeshOfBars(pmin,pmax,npoints,secondOrder=False):

    points = np.zeros((npoints,3))
    points[:,0] = np.linspace(pmin,pmax,npoints)

    if secondOrder :
        if npoints % 2 == 0:
            raise(Exception("the number of point must be odd in secondOrder"))

        bars = np.empty(((npoints-1)//2 ,3))
        bars[:,0] = np.arange(0,npoints-2,2)
        bars[:,1] = np.arange(2,npoints,2)
        bars[:,2] = np.arange(1,npoints-1,2)
        res = CreateMeshOf(points,bars,elemName = ElementNames.Bar_3 )
        #print(bars)
        #raise
    else:
        bars = np.empty((npoints-1,2))
        bars[:,0] = np.arange(npoints-1)
        bars[:,1] = np.arange(1,npoints)
        res = CreateMeshOf(points,bars,elemName = ElementNames.Bar_2 )

    elements = res.GetElementsOfType(ElementNames.Point_1)
    elements.connectivity = np.array([[0],[npoints-1]],dtype=PBasicIndexType)
    elements.originalIds = np.arange(2,dtype=PBasicIndexType)
    elements.cpt = elements.connectivity.shape[0]
    elements.tags.CreateTag("L",).SetIds([0])
    elements.tags.CreateTag("H",).SetIds([1])
    res.PrepareForOutput()
    return res

def CreateMeshOfTriangles(points,tris):
    return CreateMeshOf(points,tris,elemName = ElementNames.Triangle_3 )

def CreateMeshOf(points,connectivity,elemName = None,out=None):

    if elemName is None:
        raise Exception("Need a element name ")# pragma: no cover

    if out is None:
        res = UnstructuredMesh()
    else:
        res = out # pragma: no cover

    res.nodes = np.array(points, dtype=np.double)
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=PBasicIndexType)

    elements = res.GetElementsOfType(elemName)
    elements.connectivity = np.array(connectivity,dtype=PBasicIndexType)
    elements.originalIds = np.arange(0,elements.connectivity.shape[0],dtype=PBasicIndexType)
    elements.cpt = elements.connectivity.shape[0]
    elements.tags.CreateTag(str(ElementNames.dimension[elemName])+"D").SetIds(np.arange(elements.GetNumberOfElements() ) )
    res.PrepareForOutput()

    return res

def CreateSquare(dimensions=[2,2], origin=[-1.0,-1.0], spacing=[1.,1.], ofTris=False):
    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)
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


    mesh = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=ofTris)
    skin = ComputeSkin(mesh)
    for name,data in skin.elements.items():
        mesh.GetElementsOfType(name).Merge(data)
    #print(skin)

    if ofTris:
        tris = mesh.GetElementsOfType(EN.Triangle_3)
        tris.GetTag("2D").SetIds(np.arange(tris.GetNumberOfElements()))
    else:
        quads = mesh.GetElementsOfType(EN.Quadrangle_4)
        quads.GetTag("2D").SetIds(np.arange(quads.GetNumberOfElements()))
    skin = mesh.GetElementsOfType(EN.Bar_2)
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

def CreateDisk(nr=10, ntheta=10, r0=0.5, r1=1, theta0=0, theta1=np.pi/2 ):
    """Functionto create a disk section
    nr: number of points in the radial direction
    ntheta: number of point in the angular direction
    r0: internal radius
    r1: external radius
    theta0: start angle
    theta1: end angle

    return : Unstructured mesh of a disk sector
    """
    myMesh = CreateSquare(dimensions=[nr,ntheta],origin=[r0,theta0],spacing=[(r1-r0)/(nr-1),(theta1-theta0)/(ntheta-1)])

    r = myMesh.nodes[:,0].copy()
    theta = myMesh.nodes[:,1].copy()

    myMesh.nodes[:,0] = r*np.cos(theta)
    myMesh.nodes[:,1] = r*np.sin(theta)

    import BasicTools.Containers.ElementNames as EN
    myMesh.elements[EN.Bar_2].tags["X0"].name = "R0"
    myMesh.elements[EN.Bar_2].tags["X1"].name = "R1"
    myMesh.elements[EN.Bar_2].tags["Y0"].name = "Theta0"
    myMesh.elements[EN.Bar_2].tags["Y1"].name = "Theta1"

    return myMesh


def CreateCube(dimensions=[2,2,2], origin=[-1.0,-1.0,-1.0], spacing=[1.,1.,1.], ofTetras=False):
    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshModificationTools import ComputeSkin
    import BasicTools.Containers.ElementNames as EN

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions(dimensions);
    myMesh.SetOrigin(origin);
    myMesh.SetSpacing(spacing);

    # coorners
    d = np.array(dimensions)-1
    s = spacing
    indexs = [[   0,   0,   0],
              [d[0],   0,   0],
              [   0,d[1],   0],
              [d[0],d[1],   0],
              [   0,   0,d[2]],
              [d[0],   0,d[2]],
              [   0,d[1],d[2]],
              [d[0],d[1],d[2]]]
    for n in indexs:
        idx = myMesh.GetMonoIndexOfNode(n)
        name = "x"  + ("0" if n[0]== 0 else "1" )
        name += "y" + ("0" if n[1]== 0 else "1" )
        name += "z" + ("0" if n[2]== 0 else "1" )
        myMesh.nodesTags.CreateTag(name,False).SetIds([idx])


    mesh = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=ofTetras)
    skin = ComputeSkin(mesh)
    for name,data in skin.elements.items():
        mesh.GetElementsOfType(name).Merge(data)


    if ofTetras:
        tets = mesh.GetElementsOfType(EN.Tetrahedron_4)
        tets.GetTag("3D").SetIds(range(tets.GetNumberOfElements()))
        skin = mesh.GetElementsOfType(EN.Triangle_3)
    else:
        hexs = mesh.GetElementsOfType(EN.Hexaedron_8)
        hexs.GetTag("3D").SetIds(range(hexs.GetNumberOfElements()))
        skin = mesh.GetElementsOfType(EN.Quadrangle_4)
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

def MeshToSimplex(mesh, legacy=False):
    """
    Convert mesh to only tetrahedron/triangle/bars/points
    """

    from BasicTools.Containers.UnstructuredMesh import ElementsContainer,AllElements

    ae = AllElements()


    for elemtype, data in mesh.elements.items():
        res = data

        if elemtype == ElementNames.Hexaedron_8:

            res = ElementsContainer(ElementNames.Tetrahedron_4)
            nbelems = data.GetNumberOfElements()
            res.Allocate(nbelems*6)
            conn = data.connectivity
            res.connectivity[0:nbelems*6:6] = conn[:,[0, 6, 2, 3]]
            res.connectivity[1:nbelems*6:6] = conn[:,[0, 6, 3, 7]]
            res.connectivity[2:nbelems*6:6] = conn[:,[0, 6, 7, 4]]
            res.connectivity[3:nbelems*6:6] = conn[:,[0, 6, 4, 5]]
            res.connectivity[4:nbelems*6:6] = conn[:,[0, 6, 5, 1]]
            res.connectivity[5:nbelems*6:6] = conn[:,[0, 6, 1, 2]]

            res.originalIds =  np.repeat(data.originalIds,6)
            for tname in data.tags:
                ids = data.tags[tname].GetIds()
                res.tags.CreateTag(tname).SetIds(np.repeat(ids,6)*6+np.tile(range(6),len(ids)) )

        elif elemtype == ElementNames.Quadrangle_4:

            res = ElementsContainer(ElementNames.Triangle_3)
            nbelems = data.GetNumberOfElements()
            res.Allocate(nbelems*2)
            conn = data.connectivity
            res.connectivity[0:nbelems*2:2] = conn[:,[0, 1, 2]]
            res.connectivity[1:nbelems*2:2] = conn[:,[0, 2, 3]]

            res.originalIds =  np.repeat(data.originalIds,2)
            for tname in data.tags:
                ids = data.tags[tname].GetIds()
                res.tags.CreateTag(tname).SetIds(np.repeat(ids,2)*2+np.tile(range(2),len(ids)) )
        elif elemtype in [ElementNames.Triangle_3,ElementNames.Triangle_6,ElementNames.Tetrahedron_4,ElementNames.Tetrahedron_10,ElementNames.Bar_2,ElementNames.Bar_3,ElementNames.Point_1]  :
            pass
        else:
            raise(Exception("Dont know how to convert {} to simplices".format(elemtype)))

        if elemtype in ae:
            ae[res.elementType].merge(res)
        else:
            ae[res.elementType] = res

    mesh.elements = ae

def CreateMeshFromConstantRectilinearMesh(CRM, ofTetras= False,out=None, legacy=False):
    if out is None:
        res = UnstructuredMesh()
    else:
        res = out # pragma: no cover

    res.CopyProperties(CRM)

    res.nodes = CRM.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=PBasicIndexType)
    res.nodesTags = CRM.nodesTags

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearElementContainer
    from BasicTools.Containers.UnstructuredMesh import ElementsContainer

    for elemtype, data in CRM.elements.items():
        if isinstance(data,ConstantRectilinearElementContainer):
            eres = ElementsContainer(elemtype)
            eres.connectivity = data.connectivity
            eres.tags = data.tags
            eres.cpt = data.GetNumberOfElements()
            eres.originalIds = np.arange(0,data.GetNumberOfElements(),dtype=PBasicIndexType)
        else:
            eres = data

        res.elements[elemtype] = eres

    if ofTetras:
        MeshToSimplex(res)
    return res

def QuadToLin(inputmesh, divideQuadElements=True,lineariseMiddlePoints=False):
    from BasicTools.Containers.UnstructuredMeshFieldOperations import QuadFieldToLinField

    res = type(inputmesh)()
    res.CopyProperties(inputmesh)

    res.nodes = inputmesh.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=PBasicIndexType)
    import copy
    res.nodesTags = copy.deepcopy(inputmesh.nodesTags)

    for elementName in inputmesh.elements:
        quadElement = inputmesh.elements[elementName]
        if elementName == ElementNames.Tetrahedron_10:

            lineelements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 8
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*8)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*8
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*8+0,:] = quadConn[[0,4,6,7]];
                    lineelements.connectivity[initNbElem+i*8+1,:] = quadConn[[1,5,4,8]];
                    lineelements.connectivity[initNbElem+i*8+2,:] = quadConn[[2,6,5,9,]];
                    lineelements.connectivity[initNbElem+i*8+3,:] = quadConn[[7,8,9,3]];
                    lineelements.connectivity[initNbElem+i*8+4,:] = quadConn[[4,5,6,7]];
                    lineelements.connectivity[initNbElem+i*8+5,:] = quadConn[[4,5,7,8]];
                    lineelements.connectivity[initNbElem+i*8+6,:] = quadConn[[5,6,7,9]];
                    lineelements.connectivity[initNbElem+i*8+7,:] = quadConn[[5,7,8,9]];
                    if lineariseMiddlePoints :
                        res.nodes[quadConn[4],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                        res.nodes[quadConn[5],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[6],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[0],:] )/2
                        res.nodes[quadConn[7],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[3],:] )/2
                        res.nodes[quadConn[8],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[3],:] )/2
                        res.nodes[quadConn[9],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[3],:] )/2
            else:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*1)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*1
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3]];


        elif elementName == ElementNames.Triangle_6:

            lineelements = res.GetElementsOfType(ElementNames.Triangle_3)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 4
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*4)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*4
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*4+0,:] = quadConn[[0,3,5]];
                    lineelements.connectivity[initNbElem+i*4+1,:] = quadConn[[1,4,3]];
                    lineelements.connectivity[initNbElem+i*4+2,:] = quadConn[[2,5,4]];
                    lineelements.connectivity[initNbElem+i*4+3,:] = quadConn[[3,4,5]];
                    if lineariseMiddlePoints :
                        res.nodes[quadConn[3],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                        res.nodes[quadConn[4],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                        res.nodes[quadConn[5],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[0],:] )/2
            else:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2]];

        elif elementName == ElementNames.Quadrangle_8:

            lineelements = res.GetElementsOfType(ElementNames.Quadrangle_4)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*1)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*1
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*1+0,:] = quadConn[[0,1,2,3]];
                    #if lineariseMiddlePoints :
                    #    res.nodes[quadConn[4],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
                    #    res.nodes[quadConn[5],:] = (res.nodes[quadConn[1],:] + res.nodes[quadConn[2],:] )/2
                    #    res.nodes[quadConn[6],:] = (res.nodes[quadConn[2],:] + res.nodes[quadConn[3],:] )/2
                    #    res.nodes[quadConn[7],:] = (res.nodes[quadConn[3],:] + res.nodes[quadConn[0],:] )/2
            else:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3]];

        elif elementName == ElementNames.Hexaedron_20:

            lineelements = res.GetElementsOfType(ElementNames.Hexaedron_8)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*1)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*1
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*1+0,:] = quadConn[[0,1,2,3,4,5,6,7]];
                    if lineariseMiddlePoints :
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
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3,4,5,6,7]];


        elif elementName == ElementNames.Hexaedron_27:

            lineelements = res.GetElementsOfType(ElementNames.Hexaedron_8)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 8
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*8)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*8
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*8+0,:] = quadConn[[0,8,24,11,16,22,26,20]];
                    lineelements.connectivity[initNbElem+i*8+1,:] = quadConn[[8,0,9,24,22,17,21,26]];
                    lineelements.connectivity[initNbElem+i*8+2,:] = quadConn[[11,24,10,3,20,26,23,19]];
                    lineelements.connectivity[initNbElem+i*8+3,:] = quadConn[[24,9,2,10,26,21,18,23]];
                    lineelements.connectivity[initNbElem+i*8+4,:] = quadConn[[16,22,26,20,4,12,25,15]];
                    lineelements.connectivity[initNbElem+i*8+5,:] = quadConn[[22,17,21,26,12,5,13,25]];
                    lineelements.connectivity[initNbElem+i*8+6,:] = quadConn[[20,26,23,19,15,25,14,7]];
                    lineelements.connectivity[initNbElem+i*8+7,:] = quadConn[[26,21,18,23,25,13,6,14]];

                    if lineariseMiddlePoints :
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
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3,4,5,6,7]];


        elif elementName == ElementNames.Bar_3:

            lineelements = res.GetElementsOfType(ElementNames.Bar_2)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 2
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*2)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*2
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*2+0,:] = quadConn[[0,2]];
                    lineelements.connectivity[initNbElem+i*2+1,:] = quadConn[[2,1]];
                    if lineariseMiddlePoints :
                        res.nodes[quadConn[2],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
            else:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in range(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1]];
        elif ElementNames.linear[elementName] :
            lineelements = res.GetElementsOfType(elementName)
            initNbElem = lineelements.GetNumberOfElements();

            lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
            nbOfNewElements = 1
            lineelements.connectivity[initNbElem:initNbElem+quadElement.GetNumberOfElements(),:] = quadElement.connectivity
            lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()

        else:
            raise Exception('Error : not coded yet for this type of elements ' + str(elementName))# pragma: no cover
        #copy of tags
        for originaltag in quadElement.tags :
            destinationtag = lineelements.GetTag(originaltag.name)
            ids = originaltag.GetIds()
            for i in range(originaltag.cpt):
                for t in range(nbOfNewElements):
                    destinationtag.AddToTag(initNbElem+ids[i]*nbOfNewElements+t)

            destinationtag.Tighten()

        res.ComputeGlobalOffset()

    res.PrepareForOutput()

    if divideQuadElements == False:
        CleanLonelyNodes(res)

    return res

def MirrorMesh(inmesh,x=None,y=None,z=None) :
    nbpoints = inmesh.GetNumberOfNodes()

    outmesh = type(inmesh)()
    outmesh.CopyProperties(inmesh)

    d = 0
    if x is not None:
        d += 1
    if y is not None:
        d += 1
    if z is not None:
        d += 1

    outmesh.nodes = np.empty((nbpoints*(2**d),inmesh.GetDimensionality()), dtype=PBasicFloatType)
    outmesh.originalIDNodes = np.empty((nbpoints*(2**d),), dtype=PBasicIndexType)

    #copy of points:
    outmesh.nodes[0:nbpoints,:] = inmesh.nodes
    #copy of points:
    outmesh.originalIDNodes[0:nbpoints] = inmesh.originalIDNodes
    import copy
    outmesh.nodesTags = copy.deepcopy(inmesh.nodesTags)
    cpt = nbpoints

    def increaseTags(tags,oldSize):
        for tag in tags:
            ids = tag.GetIds()[:]  # make a copy
            tag.SetIds(np.hstack((ids,ids+oldSize)) )

    if x is not None:
        vec = np.array([ [  -1,1,1], ],dtype=float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[x,0,0])*vec+ [x,0,0]
        outmesh.originalIDNodes[cpt:(2*cpt)] = outmesh.originalIDNodes[0:cpt]
        increaseTags(outmesh.nodesTags,cpt)
        cpt = cpt*2

    if y is not None:
        vec = np.array([ [  1,-1,1], ],dtype=float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[0,y,0])*vec+ [0,y,0]
        outmesh.originalIDNodes[cpt:(2*cpt)] = outmesh.originalIDNodes[0:cpt]
        increaseTags(outmesh.nodesTags,cpt)
        cpt = cpt*2

    if z is not None:
        vec = np.array([ [  1,1,-1], ],dtype=float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[0,0,z])*vec+ [0,0,z]
        outmesh.originalIDNodes[cpt:(2*cpt)] = outmesh.originalIDNodes[0:cpt]
        increaseTags(outmesh.nodesTags,cpt)
        cpt = cpt*2

    for name,vals in inmesh.elements.items():
        nbelements = vals.GetNumberOfElements()
        outelements = outmesh.GetElementsOfType(name)
        outelements.Reserve(nbelements*(2**d))
        outelements.connectivity[0:nbelements,:] = vals.connectivity
        outelements.tags = copy.deepcopy(vals.tags)
        cpt = nbelements
        pcpt = nbpoints
        permutation = ElementNames.mirrorPermutation[name]
        if x is not None:
            outelements.connectivity[cpt:(2*cpt),:] = (outelements.connectivity[0:cpt,:]+pcpt)[:,permutation]
            pcpt = pcpt *2
            increaseTags(outelements.tags,cpt)
            cpt = cpt*2

        if y is not None:
            outelements.connectivity[cpt:(2*cpt),:] = (outelements.connectivity[0:cpt,:]+pcpt)[:,permutation]
            pcpt = pcpt *2
            increaseTags(outelements.tags,cpt)
            cpt = cpt*2

        if z is not None:
            outelements.connectivity[cpt:(2*cpt),:] = (outelements.connectivity[0:cpt,:]+pcpt)[:,permutation]
            pcpt = pcpt *2
            increaseTags(outelements.tags,cpt)
            cpt = cpt*2

        outelements.cpt = cpt

    outmesh.PrepareForOutput()
    return outmesh

def Creat0DElementAtEveryPoint(mesh):
    res   = UnstructuredMesh()
    res.nodes = mesh.nodes
    res.originalIDNodes = mesh.originalIDNodes
    res.nodesTags = mesh.nodesTags

    elems1D = res.elements.GetElementsOfType(ElementNames.Point_1)
    elems1D.tags = mesh.nodesTags
    elems1D.connectivity = np.arange(mesh.GetNumberOfNodes())
    elems1D.connectivity.shape  = (mesh.GetNumberOfNodes(),1)
    elems1D.originalIds = mesh.originalIDNodes
    elems1D.cpt = mesh.GetNumberOfNodes()
    return res

def SubDivideMesh(mesh,level=1):
    if level == 0:
        return mesh

    res = UnstructuredMesh()

    subdivitionAlmanac = {}
    subdivitionAlmanac[ElementNames.Point_1] = [(ElementNames.Point_1, [0]) ]
    subdivitionAlmanac[ElementNames.Bar_2] = [(ElementNames.Bar_2, [0,2]),
                                              (ElementNames.Bar_2, [2,1])]

    subdivitionAlmanac[ElementNames.Quadrangle_4] = [(ElementNames.Quadrangle_4, [0,4,8,7]),
                                                     (ElementNames.Quadrangle_4, [4,1,5,8]),
                                                     (ElementNames.Quadrangle_4, [8,5,2,6]),
                                                     (ElementNames.Quadrangle_4, [7,8,6,3])]

    subdivitionAlmanac[ElementNames.Hexaedron_8] = [(ElementNames.Hexaedron_8, [ 0, 8,24,11,16,22,26,20]),
                                                    (ElementNames.Hexaedron_8, [ 8, 1, 9,24,22,17,21,26]),
                                                    (ElementNames.Hexaedron_8, [24, 9, 2,10,26,21,18,23]),
                                                    (ElementNames.Hexaedron_8, [11,24,10, 3,20,26,23,19]),
                                                    (ElementNames.Hexaedron_8, [16,22,26,20, 4,12,25,15]),
                                                    (ElementNames.Hexaedron_8, [22,17,21,26,12, 5,13,25]),
                                                    (ElementNames.Hexaedron_8, [26,21,18,23,25,13, 6,14]),
                                                    (ElementNames.Hexaedron_8, [20,26,23,19,15,25,14, 7])]

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo, LagrangeSpaceP2
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.IntegrationsRules import NodalEvaluationP2
    numberingGeo = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)
    numberingP2 = ComputeDofNumbering(mesh,LagrangeSpaceP2)

    ## Generation of nodes
    res.nodes = np.empty((numberingP2.size,3), dtype=PBasicFloatType)
    res.originalIDNodes = np.zeros(res.nodes.shape[0],dtype=PBasicIndexType)-1

    oldToNewDofs = []

    for i in range(mesh.GetNumberOfNodes()):
        oldToNewDofs.append(numberingP2.GetDofOfPoint(i))

    oldToNewDofs = np.array(oldToNewDofs,dtype=PBasicIndexType)
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
        for t,nn in subdivitionAlmanac[elemType]:
            #t = elementype
            #nn = new numbering
            nelems = res.GetElementsOfType(t)
            offset = nelems.GetNumberOfElements()
            tne =  nelems.AddNewElements(nP2[:,nn], data.originalIds)
            for tag in data.tags.keys():
                name = data.tags[tag].name
                ids = data.tags[tag].GetIds()
                if len(ids) == 00:
                    continue
                nelems.tags.CreateTag(name,False).AddToTag(ids+offset)
                #nelems.tags.CreateTag(name,False).AddToTag(np.arange(offset,tne))

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
    hexa= [[0,1,2,3,4,5,6,7],]
    mesh = CreateMeshOf(points,hexa,ElementNames.Hexaedron_8)
    mesh.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    mesh.GetElementsOfType(ElementNames.Hexaedron_8)
    mesh.GetElementsOfType(ElementNames.Hexaedron_8).tags
    mesh.GetElementsOfType(ElementNames.Hexaedron_8).tags.CreateTag("OnlyHex")
    mesh.GetElementsOfType(ElementNames.Hexaedron_8).tags.CreateTag("OnlyHex",False).AddToTag(0)

    outmesh = SubDivideMesh(mesh,1)

    print(mesh)
    print(outmesh)
    if GUI:
        from BasicTools.Containers.vtkBridge import PlotMesh
        PlotMesh(outmesh)

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    WriteMeshToXdmf(tempdir+"CheckIntegrity_SubDivideMesh.xdmf",outmesh,PointFields=[outmesh.originalIDNodes],PointFieldsNames=["originalIDNodes"] )
    print(tempdir)

    if outmesh.GetNumberOfNodes() != 27:
        raise# pragma: no cover

    if outmesh.GetNumberOfElements() != 8:
        raise# pragma: no cover
    return "ok"


def CheckIntegrity_Creat0DElementAtEveryPoint(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("FirstPoint").AddToTag(0)

    print(mesh)
    outputmesh = Creat0DElementAtEveryPoint(mesh)
    print(outputmesh)
    if outputmesh.elements.GetElementsOfType(ElementNames.Point_1).GetNumberOfElements() != mesh.GetNumberOfNodes():
        raise(Exception("Error "))
    return "ok"

def CheckIntegrity_MirrorMesh(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    mesh.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    mesh.GetElementsOfType(ElementNames.Tetrahedron_4).tags.CreateTag("OnlyTet").AddToTag(0)
    outmesh = MirrorMesh(mesh,x=0)

    if outmesh.GetNumberOfNodes() != 8:
        raise# pragma: no cover

    if outmesh.GetNumberOfElements() != 2:
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
    print(QuadToLin(myMesh,divideQuadElements=True,lineariseMiddlePoints=True))
    from BasicTools.Containers.UnstructuredMeshFieldOperations import QuadFieldToLinField
    QuadFieldToLinField(myMesh, np.arange(myMesh.GetNumberOfNodes()))
    QuadFieldToLinField(myMesh, np.arange(myMesh.GetNumberOfNodes()), linMesh)
    return "ok"


def CheckIntegrity_CreateMeshFromConstantRectilinearMesh(GUI=False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(2)
    myMesh.SetDimensions([3,3]);
    myMesh.SetSpacing([1, 1]);
    print(myMesh)
    print(CreateMeshFromConstantRectilinearMesh(myMesh))

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([3,3,3]);
    myMesh.SetSpacing([1, 1,1]);
    print(myMesh)
    print(CreateMeshFromConstantRectilinearMesh(myMesh))

    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=True)
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
    mesh = CreateSquare(dimensions=[20,21],spacing=[2.,2.],ofTris=False)
    mesh = CreateSquare(dimensions=[20,21],spacing=[2.,2.],ofTris=True)
    return "ok"

def CheckIntegrity_CreateUniformMeshOfBars(GUI=False):
    print(CreateUniformMeshOfBars(0,8,10))
    return "ok"

def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_CreateUniformMeshOfBars,
    CheckIntegrity_CreateCube,
    CheckIntegrity_CreateSquare,
    CheckIntegrity_CreateMeshOfTriangles,
    CheckIntegrity_CreateMeshFromConstantRectilinearMesh,
    CheckIntegrity_QuadToLin,
    CheckIntegrity_MirrorMesh,
    CheckIntegrity_Creat0DElementAtEveryPoint,
    CheckIntegrity_SubDivideMesh,
    CheckIntegrity_CreateDisk,
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
