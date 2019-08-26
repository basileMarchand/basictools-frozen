# -*- coding: utf-8 -*-
import numpy as np
import math

import BasicTools.Helpers.BaseOutputObject as BaseOutputObject

from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
import BasicTools.Containers.ElementNames as ElementNames

def CreateUniformMeshOfBars(pmin,pmax,npoints):
    points = np.zeros((npoints,3))
    points[:,0] = np.linspace(pmin,pmax,npoints)
    bars = np.empty((npoints-1,2))
    bars[:,0] = np.arange(npoints-1)
    bars[:,1] = np.arange(1,npoints)
    return CreateMeshOf(points,bars,elemName = ElementNames.Bar_2 )

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
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int)

    elements = res.GetElementsOfType(elemName)
    elements.connectivity = np.array(connectivity,dtype=np.int)
    elements.originalIds = np.arange(0,elements.connectivity.shape[0],dtype=np.int)
    elements.cpt = elements.connectivity.shape[0]
    elements.tags.CreateTag(str(ElementNames.dimension[elemName])+"D").SetIds(np.arange(elements.GetNumberOfElements() ) )
    res.PrepareForOutput()

    return res

def CreateSquare(dimensions=[2,2], origin=[-1.0,-1.0], spacing=[1.,1.], ofTris=False):
    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshFromConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshTools import ComputeSkin
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
        tris.GetTag("2D").SetIds(range(tris.GetNumberOfElements()))
    else:
        quads = mesh.GetElementsOfType(EN.Quadrangle_4)
        quads.GetTag("2D").SetIds(range(quads.GetNumberOfElements()))
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


def CreateCube(dimensions=[2,2,2], origin=[-1.0,-1.0,-1.0], spacing=[1.,1.,1.], ofTetras=False):
    spacing = np.array(spacing,dtype=float)
    origin = np.array(origin,dtype=float)
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshFromConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshTools import ComputeSkin
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


def CreateMeshFromConstantRectilinearMesh(CRM, ofTetras= False,out=None):
    if out is None:
        res = UnstructuredMesh()
    else:
        res = out # pragma: no cover

    res.CopyProperties(CRM)

    res.nodes = CRM.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int);

    res.nodesTags = CRM.nodesTags

    nbelements = CRM.GetNumberOfElements()

    elementtype = ElementNames.Tetrahedron_4
    if(CRM.GetDimensionality() == 3):
        if ofTetras:
            elementtype = ElementNames.Tetrahedron_4
            nbelements = CRM.GetNumberOfElements()*6
        else:
            elementtype = ElementNames.Hexaedron_8
    else:
        if ofTetras:
            elementtype = ElementNames.Triangle_3
            nbelements = CRM.GetNumberOfElements()*2
        else:
            elementtype = ElementNames.Quadrangle_4

    elements = res.GetElementsOfType(elementtype)
    elements.connectivity = np.zeros((nbelements,ElementNames.numberOfNodes[elementtype]),dtype=np.int)
    elements.cpt = nbelements

    if ofTetras:
        if(CRM.GetDimensionality() == 3):
            p0 = np.array([0,1,2,3,4,5,6,7])
            p1=  np.array([1,2,3,0,5,6,7,4])
            p2=  np.array([3,0,1,2,7,4,5,6])
            p3=  np.array([2,3,0,1,6,7,4,5])

            for elem in range(CRM.GetNumberOfElements()):

                index = CRM.GetMultiIndexOfElement(elem)
                idx = index[0]%2+ 2*(index[1]%2)+4*(index[2]%2)
                if idx == 0:
                    per = p0
                elif idx == 1:
                    per = p1
                elif idx == 2:
                    per = p2
                elif idx == 3:
                    per = p3
                elif idx == 4:
                    per = p3
                elif idx == 5:
                    per = p2
                elif idx == 6:
                    per = p1
                elif idx == 7:
                    per = p0
                else:
                    raise # pragma: no cover

                conn = CRM.GetConnectivityForElement(elem)
                elements.connectivity[elem*6+0,:] = conn[per[[0,6,5,1]]];
                elements.connectivity[elem*6+1,:] = conn[per[[0,6,1,2]]];
                elements.connectivity[elem*6+2,:] = conn[per[[0,6,2,3]]];
                elements.connectivity[elem*6+3,:] = conn[per[[0,6,3,7]]];
                elements.connectivity[elem*6+4,:] = conn[per[[0,6,7,4]]];
                elements.connectivity[elem*6+5,:] = conn[per[[0,6,4,5]]];
        else:
            for elem in range(CRM.GetNumberOfElements()):
                conn = CRM.GetConnectivityForElement(elem)
                elements.connectivity[elem*2+0,:] = conn[[0, 1, 2]];
                elements.connectivity[elem*2+1,:] = conn[[0, 2, 3]];

    else:
        elements.connectivity = CRM.GenerateFullConnectivity()

    elements.originalIds = np.arange(0,elements.GetNumberOfElements(),dtype=np.int)
    res.PrepareForOutput()
    return res

def QuadToLin(inputmesh, divideQuadElements=True,lineariseMiddlePoints=False):

    res = type(inputmesh)()
    res.CopyProperties(inputmesh)

    res.nodes = inputmesh.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int);
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


    if divideQuadElements == False:
        CleanLonelyNodes(res)

    res.PrepareForOutput()
    return res

def CleanDoubleNodes(res, tol = None, nodesToTestMask= None):

    BaseOutputObject.BaseOutputObject().PrintDebug("in CleanDoubleNodes")

    res.ComputeBoundingBox()
    if tol is None:
        tol = np.linalg.norm(res.boundingMax - res.boundingMin)*1e-7

    nbnodes = res.GetNumberOfNodes()
    toKeep = np.zeros(nbnodes, dtype=np.bool )
    newindex = np.zeros(nbnodes, dtype=np.int )
    tol2 = tol**2
    def dist2(array,value):
        d = array-value
        return np.inner(d,d)

    if nodesToTestMask is None and tol == 0:
        # optimized version for tol = 0.0
        database = {}
        cpt  = 0
        for i in range(nbnodes):
            point = tuple(res.nodes[i,:])
            ni = database.get(point)
            if ni is None:
                database[point] = cpt
                newindex[i] = cpt
                cpt +=1
                toKeep[i] = True;
            else:
                newindex[i] = ni

    elif nodesToTestMask is None :

        from BasicTools.Containers.Octree import Octree
        cpt  = 0

        Ma = res.boundingMax
        Mi = res.boundingMin
        diag2 = np.linalg.norm(res.boundingMax-res.boundingMin)**2
        tree = Octree(Ma[0],Ma[1],Ma[2], Mi[0], Mi[1], Mi[2])

        for i in range(nbnodes):
            #if float(i)/1000 == int(i/1000):
                #print(100.*i/nbnodes)
            point = tuple(res.nodes[i,:])
            entries = tree.find_within_range_cube(point, tol)

            odist = diag2
            if len(entries):
                for entry in entries:
                    dist = (point[0] - entry[0][0])**2   +(point[1] - entry[0][1])**2+(point[2] - entry[0][2])**2
                    if dist < odist:
                        odist = dist
                        index = entry[1]

            if odist <= tol2:
                newindex[i] = index
            else:
                tree.add_item(cpt, point )
                newindex[i] = cpt
                cpt += 1
                toKeep[i] = True;


    else:
        cpt =0
        for i in range(nbnodes):
            if not nodesToTestMask[i]:
                  newindex[i] = cpt
                  cpt += 1
                  toKeep[i] = True;
                  continue

            posi = res.nodes[i,:]

            for j in range(i):
                if not nodesToTestMask[j]:
                    continue
                #dist =np.linalg.norm(res.nodes[i,:]-res.nodes[j,:])
                if toKeep[j]:
                    dist = dist2(posi,res.nodes[j,:] )
                    if dist < tol2 :
                        newindex[i] = newindex[j]
                        break
            else:
                newindex[i] = cpt
                cpt += 1
                toKeep[i] = True;


    res.nodes = res.nodes[toKeep,:]
    res.originalIDNodes = np.where(toKeep)[0]

    for tag in res.nodesTags :
        tag.SetIds(np.unique(newindex[tag.GetIds()]) )


    for elementName in res.elements:
        elements = res.elements[elementName]
        elements.connectivity = newindex[elements.connectivity]

    BaseOutputObject.BaseOutputObject().PrintDebug("CleanDoubleNodes Done")

def CleanLonelyNodes(res,out=None):

    usedNodes = np.zeros(res.GetNumberOfNodes(),dtype=np.bool )
    for elementName in res.elements:
        elements = res.elements[elementName]
        usedNodes[elements.connectivity.ravel()] = True;

    cpt = 0 ;
    NewIndex =  np.zeros(res.GetNumberOfNodes(),dtype=np.int )-1
    originalIDNodes = np.zeros(res.GetNumberOfNodes(),dtype=np.int)
    for n in range(res.GetNumberOfNodes()):
        if usedNodes[n]:
            NewIndex[n] = cpt
            originalIDNodes[cpt] = n
            cpt += 1


    #filter the nodes

    #inplace
    if out is None:
        res.nodes = res.nodes[usedNodes ,:]
        res.originalIDNodes = np.where(usedNodes)[0]
        #node tags
        for tag in res.nodesTags :
            tag.SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])

                #renumbering the connectivity matrix
        for elementName in res.elements:
            elements = res.elements[elementName]
            elements.connectivity = NewIndex[elements.connectivity]

    else:
        out.nodes = res.nodes[usedNodes ,:]
        out.originalIDNodes = np.where(usedNodes)[0]
        #node tags
        for tag in res.nodesTags :
            outtag = out.nodesTags.CreateTag(tag.name)
            outtag.SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])

        import copy
        for elementName in res.elements:
            elements = res.elements[elementName]
            outelements = out.GetElementsOfType(elementName)
            outelements.connectivity = NewIndex[elements.connectivity]

            outelements.tags = copy.deepcopy(elements.tags)

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

    outmesh.nodes = np.empty((nbpoints*(2**d),inmesh.GetDimensionality()), dtype=np.float)

    #copy of points:
    outmesh.nodes[0:nbpoints,:] = inmesh.nodes
    import copy
    outmesh.nodesTags = copy.deepcopy(inmesh.nodesTags)
    cpt = nbpoints

    def increaseTags(tags,oldSize):
        for tag in tags:
            ids = tag.GetIds()[:]  # make a copy
            tag.SetIds(np.hstack((ids,ids+oldSize)) )


    if x is not None:
        vec = np.array([ [  -1,1,1], ],dtype=np.float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[x,0,0])*vec+ [x,0,0]
        increaseTags(outmesh.nodesTags,cpt)
        cpt = cpt*2

    if y is not None:
        vec = np.array([ [  1,-1,1], ],dtype=np.float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[0,y,0])*vec+ [0,y,0]
        increaseTags(outmesh.nodesTags,cpt)
        cpt = cpt*2

    if z is not None:
        vec = np.array([ [  1,1,-1], ],dtype=np.float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[0,0,z])*vec+ [0,0,z]
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


def ExtractElementsByMask(inelems, _mask):
    outelems = type(inelems)(inelems.elementType)

    newIndex = np.empty(inelems.GetNumberOfElements(),dtype=np.int)


    if _mask.dtype == np.bool:

        nbels =0;
        for i in range(inelems.GetNumberOfElements()):
           newIndex[i] = nbels
           nbels += 1 if _mask[i] else 0
        mask = _mask
    else:
        nbels = len(_mask)
        mask = np.zeros(inelems.GetNumberOfElements(),dtype=np.bool)
        cpt =0;
        for index in _mask:
           newIndex[index ] = cpt
           mask[index] = True
           cpt += 1

    outelems.Allocate(nbels)
    outelems.connectivity = inelems.connectivity[mask,:]
    outelems.originalIds = np.where(mask)[0]
    outelems.originalOffset = inelems.globaloffset

    for tag in inelems.tags  :
       temp = np.extract(mask[tag.GetIds()],tag.GetIds())
       newid = newIndex[temp]
       outelems.tags.CreateTag(tag.name).SetIds(newid)

    return outelems

def ExtractElementByDimensionalityNoCopy(inmesh,dimensionalityFilter):
    outmesh = type(inmesh)()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = inmesh.nodes
    outmesh.originalIDNodes = inmesh.originalIDNodes
    outmesh.nodesTags = inmesh.nodesTags

    if dimensionalityFilter >0 :
        for name,elems in inmesh.elements.items():
            if ElementNames.dimension[name] != dimensionalityFilter:
                continue
            else:
                outmesh.elements[name] = elems
    else:
        for name,elems in inmesh.elements.items():
            if ElementNames.dimension[name] == -dimensionalityFilter:
                continue
            else:
                outmesh.elements[name] = elems

    outmesh.PrepareForOutput()
    return outmesh

def ExtractElementByTags(inmesh,tagsToKeep, allNodes=False,dimensionalityFilter= None, cleanLonelyNodes=True):

    outmesh = UnstructuredMesh()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = np.copy(inmesh.nodes)
    outmesh.originalIDNodes = None

    import copy
    for tag in inmesh.nodesTags:
        outmesh.nodesTags.AddTag(copy.deepcopy(tag) )


    nodalMask = np.zeros(inmesh.GetNumberOfNodes(),dtype = np.bool)
    for name,elems in inmesh.elements.items():

       #if dimensionalityFilter is not None:
       #    if dimensionalityFilter !=  ElementNames.dimension[name]:
       #        continue

       if (np.any([x in elems.tags for x in tagsToKeep] ) == False) and (dimensionalityFilter is None) :
           if np.any([x in inmesh.nodesTags for x in tagsToKeep]) == False:
               continue# pragma: no cover


       toKeep = np.zeros(elems.GetNumberOfElements(), dtype=np.bool)
       # check elements tags
       for tagToKeep in tagsToKeep:
           if tagToKeep in elems.tags:
               toKeep[elems.tags[tagToKeep].GetIds()] = True

       # check for nodes tags
       for tagToKeep in tagsToKeep:
           if tagToKeep in inmesh.nodesTags:
             nodalMask.fill(False)
             tag = inmesh.GetNodalTag(tagToKeep)
             nodalMask[tag.GetIds()] = True
             elemMask = np.sum(nodalMask[elems.connectivity],axis=1)
             if allNodes :
                 toKeep[elemMask == elems.GetNumberOfNodesPerElement()] = True
             else:
                 toKeep[elemMask > 0] = True

       # if dimensionality is ok and no tag to keep, we keep all elements
       if dimensionalityFilter is not None:
           if dimensionalityFilter ==  ElementNames.dimension[name]:
               toKeep[:] = True
           #if len(tagsToKeep)  == 0:
           #    toKeep[:] = True

       newIndex = np.empty(elems.GetNumberOfElements(), dtype=np.int )
       cpt =0;
       for i in range(elems.GetNumberOfElements()):
           newIndex[i] = cpt
           cpt += 1 if toKeep[i] else 0





       outelem = outmesh.GetElementsOfType(name)
       nbTokeep = np.sum(toKeep)
       outelem.Allocate(nbTokeep)
       outelem.connectivity = elems.connectivity[toKeep,:]
       outelem.originalIds = np.where(toKeep)[0]

       for tag in elems.tags  :
           temp = np.extract(toKeep[tag.GetIds()],tag.GetIds())
           newid = newIndex[temp]
           outelem.tags.CreateTag(tag.name,errorIfAlreadyCreated=False).SetIds(newid)

    if cleanLonelyNodes:
        CleanLonelyNodes(outmesh)
    outmesh.PrepareForOutput()
    return outmesh

def VolumeOfTetrahedrons(inmesh):

    elems =inmesh.GetElementsOfType(ElementNames.Tetrahedron_4)
    conn = elems.connectivity
    a = inmesh.nodes[conn[:,0],:]
    b = inmesh.nodes[conn[:,1],:]
    c = inmesh.nodes[conn[:,2],:]
    d = inmesh.nodes[conn[:,3],:]
    e = np.cross(b-d,c-d)
    f = (a-d)
    res = np.empty(elems.GetNumberOfElements(),dtype=np.float)
    for n in range(elems.GetNumberOfElements()):
        res[n] = np.abs( np.dot(f[n,:],e[n,:])  )

    return res*(1./6.)

def VolumeOfHexaedrons(inmesh):

    elems =inmesh.GetElementsOfType(ElementNames.Hexaedron_8)
    conn = elems.connectivity

    def VolumeInternal(a,b,c,d):
        e = np.cross(b-d,c-d)
        f = (a-d)
        res = np.empty(elems.GetNumberOfElements(),dtype=np.float)
        for n in range(elems.GetNumberOfElements()):
            res[n] = np.abs( np.dot(f[n,:],e[n,:])  )
        return res*(1./6.)

    res = np.zeros(elems.GetNumberOfElements(),dtype=np.float)

    p0 = inmesh.nodes[conn[:,0],:]
    p1 = inmesh.nodes[conn[:,1],:]
    p2 = inmesh.nodes[conn[:,2],:]
    p3 = inmesh.nodes[conn[:,3],:]
    p4 = inmesh.nodes[conn[:,4],:]
    p5 = inmesh.nodes[conn[:,5],:]
    p6 = inmesh.nodes[conn[:,6],:]
    p7 = inmesh.nodes[conn[:,7],:]

    res += VolumeInternal(p0,p6,p5,p1)
    res += VolumeInternal(p0,p6,p1,p2)
    res += VolumeInternal(p0,p6,p2,p3)
    res += VolumeInternal(p0,p6,p3,p7)
    res += VolumeInternal(p0,p6,p7,p4)
    res += VolumeInternal(p0,p6,p4,p5)

    return res


def GetVolume(inmesh) :

    vol = 0;
    for name,elems in inmesh.elements.items():
        if ElementNames.dimension[name] != 3:
            continue# pragma: no cover
        elif name == ElementNames.Tetrahedron_4:
            vol += np.sum(VolumeOfTetrahedrons(inmesh))
        elif name == ElementNames.Hexaedron_8:
            vol += np.sum(VolumeOfHexaedrons(inmesh))
        else:
            raise Exception('Not implemented for elements of type "'+str(name)+'" code me please...')# pragma: no cover
    return vol


def CleanEmptyTags(inmesh):
    for name,elems in inmesh.elements.items():
        elems.tags.RemoveEmptyTags()
    inmesh.nodesTags.RemoveEmptyTags()

def GetDualGraphNodeToElement(inmesh, maxNumConnections=200):
    # generation of the dual graph
    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=int )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=int );

    cpt =0

    for name,elems in inmesh.elements.items():
        for i in range(elems.GetNumberOfElements()):
            coon = elems.connectivity[i,:]
            for j in coon:
                dualGraph[j,usedPoints[j]] =  cpt
                usedPoints[j] += 1
            cpt += 1

    #we crop the output data
    maxsize = np.max(np.sum(dualGraph>=0,axis=1))
    dualGraph = dualGraph[:,0:maxsize]
    return dualGraph,usedPoints

def GetDualGraph(inmesh, maxNumConnections=200):

    # generation of the dual graph
    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=int )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=int );

    for name,elems in inmesh.elements.items():
        size = elems.GetNumberOfNodesPerElement()
        for i in range(elems.GetNumberOfElements()):
            coon = elems.connectivity[i,:]
            for j in range(size):
                myIndex = coon[j]
                for k in range(size):
                    if k == j:
                        continue
                    dualGraph[myIndex,usedPoints[myIndex]] =  coon[k]
                    usedPoints[myIndex] += 1
                    # we reached the maximun number of connection
                    # normally we have some data duplicated
                    # we try to shrink the vector
                    if usedPoints[myIndex] == maxNumConnections :# pragma: no cover
                        # normaly we shouldn't pas here
                        c = np.unique(dualGraph[myIndex,:])
                        dualGraph[myIndex,0:len(c)] = c
                        usedPoints[myIndex] = len(c)

    maxsize = 0
    # we finish now we try to compact the structure
    for i in range(inmesh.GetNumberOfNodes()):
        c = np.unique(dualGraph[i,0:usedPoints[i]])
        dualGraph[i,0:len(c)] = c
        usedPoints[i] = len(c)
        maxsize = max(len(c),maxsize)

    #we crop the output data
    dualGraph = dualGraph[:,0:maxsize]
    return dualGraph,usedPoints

# to generate one tag per body
# a body is defines by all the nodes connected by the elements
def AddTagPerBody(inmesh):

    dualGraph,usedPoints = GetDualGraph(inmesh)

    # Connectivity walk
    nbOfNodes = inmesh.GetNumberOfNodes()
    treated = np.zeros(inmesh.GetNumberOfNodes(),dtype=np.bool)
    nextpoint = np.zeros(inmesh.GetNumberOfNodes(),dtype=np.int_)

    # we start from the first point and body number 0
    nextpointcpt = 0
    bodyCpt = 0

    cpt = 0
    pointsPerBody = []
    while(True):
        # we finish all the points
        if cpt == nbOfNodes :
            break

        # we already explored all the point in this body
        if cpt == nextpointcpt:
            initialPoint = np.argmax(treated == False)
            #(edge case) in the case we have only Trues in the treated
            if treated[initialPoint]:# pragma: no cover
                break
            treated[initialPoint] = True
            nextpoint[nextpointcpt] = initialPoint
            nextpointcpt +=1

            tagName = "Body_"+str(bodyCpt)
            tag = inmesh.GetNodalTag(tagName)
            tag.AddToTag(initialPoint)
            bodyCpt += 1
            pointsInThisBody = 1


        else:
            raise Exception ("Error in this function")# pragma: no cover


        while cpt < nextpointcpt:
            workingIndex = nextpoint[cpt]
            indexes = dualGraph[workingIndex,0:usedPoints[workingIndex]]

            for index in indexes:
                if not treated[index] :
                    treated[index] = True
                    nextpoint[nextpointcpt] = index
                    nextpointcpt +=1
                    tag.AddToTag(index)
                    pointsInThisBody +=1

            cpt += 1
        pointsPerBody.append(pointsInThisBody)
    return pointsPerBody






def DeleteElements(mesh,mask):
    OriginalNumberOfElements = mesh.GetNumberOfElements()
    datatochange = {}

    mesh.ComputeGlobalOffset()

    for name,data in mesh.elements.items():
        offset = data.globaloffset
        localmask = mask[offset:offset+data.GetNumberOfElements()]
        datatochange[name] = ExtractElementsByMask(data, np.logical_not(localmask))

    for name, data in datatochange.items():
        mesh.elements[name] = data

    if OriginalNumberOfElements != mesh.GetNumberOfElements():
        print("Number Of Elements Changed: Droping elemFields")
        mesh.elemFields = {}

    mesh.PrepareForOutput()

def DeleteInternalFaces(mesh):
    OriginalNumberOfElements = mesh.GetNumberOfElements()
    skin = ComputeSkin(mesh)
    datatochange = {}

    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] != 2:
            continue
        ne = data.GetNumberOfElements()
        mask = np.zeros(ne, dtype=np.bool)

        data2 = skin.elements[name]
        ne2 =data2.GetNumberOfElements()
        surf2 = {}
        key = np.array([ne**x for x in range(ElementNames.numberOfNodes[name]) ])

        for i in range(ne2):
            cc = data2.connectivity[i,:]
            lc = np.sort(cc)
            ehash = np.sum(lc*key)
            surf2[ehash] = [1,cc]

        for i in range(ne):
            cc = data.connectivity[i,:]
            lc = np.sort(cc)

            ehash = np.sum(lc*key)
            if ehash in surf2:
                mask[i] = True

        datatochange[name] = ExtractElementsByMask(data,mask)



    for name, data in datatochange.items():
        mesh.elements[name] = data

    if OriginalNumberOfElements != mesh.GetNumberOfElements():
        print("Number Of Elements Changed: Droping elemFields")
        mesh.elemFields = {}

    mesh.PrepareForOutput()
    return mesh







    elems2D = {}
    elemsMask = {}
    usedPointsBy2DElements = []
    mesh.PrepareForOutput()
    #generate a list of all the 2D elements and an empty mask
    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] == 2:
            conn = np.sort(data.connectivity,axis=1)

            ind = np.lexsort((conn[:,1],conn[:,0], conn[:,2]))
            #print data.connectivity
            #print ind
            elems2D[name] = conn[ind,:]
            #print elems2D[name]
            elemsMask[name] = np.zeros(data.GetNumberOfElements())
            usedPointsBy2DElements = np.unique(np.hstack( (data.connectivity.ravel(),usedPointsBy2DElements) )  )
            #print elemsMask[name]
            #print usedPointsBy2DElements

    #print(usedPointsBy2DElements)
    elems3D = {}
    # generate a potential list of 3D elements touched by the 2D elements
    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] == 3:
            elems3D[name] = data
            elemsMask[name] =  np.sum(np.in1d(data.connectivity,usedPointsBy2DElements).reshape(data.connectivity.shape),axis=1) >=3
            #np.zeros(data.GetNumberOfElements())
            #print(elemsMask[name])


    for name3D,data3D in mesh.elements.items():
        if ElementNames.dimension[name3D] == 3:
            for facedata in ElementNames.faces[name3D]:
                mask = elemsMask[facedata[0]]
                conn = facedata[1]
                #print(conn)
                usedElements = data3D.connectivity[elemsMask[name3D],:]
                faces = np.sort(usedElements[:,conn],axis=1)

                ind = np.lexsort((faces [:,1],faces [:,0], faces [:,2]))
                faces  = faces [ind]
                #print(faces)
                dt  = [ ('col'+str(x),np.int) for x in range(faces.shape[1]) ]

                assert faces.flags['C_CONTIGUOUS']
                faces3D = faces.ravel().view(dt)
                #print(faces3D)

                assert elems2D[facedata[0]].flags['C_CONTIGUOUS']
                faces2D = elems2D[facedata[0]].ravel().view(dt)
                #print(elems2D[facedata[0]])
                #print(faces2D)

                mask += np.in1d(faces2D,faces3D)


    # keep only the faces with one or zero volumes attached
    newElements = {}
    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] == 2:
            mask = elemsMask[name] <=1
            newElements[name] = ExtractElementsByMask(data,mask)

    for name,data in newElements.items():
        mesh.elements[name] = data
    mesh.PrepareForOutput()


def ExtractElementsByImplicitZone(inmesh,op,allNodes=True,cellCenter=False):

    mask = op(inmesh.nodes) <= 0.

    outmesh = type(inmesh)()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = np.copy(inmesh.nodes)

    # keep only the faces with one or zero volumes attached
    for name,data in inmesh.elements.items():
        if allNodes:
            ElemMask = np.sum(mask[data.connectivity],axis=1) == data.GetNumberOfNodesPerElement()
        else:
            ElemMask = np.sum(mask[data.connectivity],axis=1) >= 1
        outmesh.elements[name] = ExtractElementsByMask(data,ElemMask)

    CleanLonelyNodes(outmesh)
    outmesh.PrepareForOutput()
    return outmesh

def ComputeSkin(mesh, md=None ,inplace=False):

    #dualGraph,usedPoints = GetDualGraph(mesh)

    if md is None:
        md = mesh.GetDimensionality()

    # first we add the 2D element to the skin almanac
    from BasicTools.Containers.Filters import ElementFilter

    surf = {}

    class Operator():
        def __init__(self,surf):
            self.surf = surf

        def __call__(self,name,data,ids):
            if name not in surf:
                self.surf[name] = {}
            surf2 = self.surf.get(name)
            #print(ids)
            for i in ids:
                cc = data.connectivity[i,:]
                lc = np.sort(cc)
                ehash = (name,tuple(lc))
                if ehash in surf2:
                    surf2[ehash][0] +=1
                    if  surf2[ehash][0] >= 2:
                        del surf2[ehash]
                else:
                    surf2[ehash] = [1,cc]

    op = Operator(surf)
    ElementFilter(mesh,dimensionality = md-1 ).ApplyOnElements(op)

    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] < md:
            continue
        faces = ElementNames.faces[name]
        ne = data.GetNumberOfElements()
        for faceType,localFaceConnectivity in faces:
            globalFaceConnectivity = data.connectivity[:,localFaceConnectivity]
            if not faceType in surf:
                surf[faceType] = {}
            surf2 = surf[faceType]
            for i in range(ne):
                cc = globalFaceConnectivity[i,:]
                lc = np.sort(cc)
                ehash = (faceType,tuple(lc))
                if ehash in surf2:
                    surf2[ehash][0] +=1
                    if  surf2[ehash][0] >= 2:
                        del surf2[ehash]
                else:
                    surf2[ehash] = [1,cc]
    if inplace:
        res = mesh
    else:
        res = UnstructuredMesh()
        res.nodes = mesh.nodes

    for name in surf:
        data = res.GetElementsOfType(name)
        surf2 = surf[name]
        for hashh,data2 in surf2.items():
            if data2[0] == 1:
                data.AddNewElement(data2[1],-1)
        data.tags.CreateTag("ExteriorSurf").SetIds(np.arange(data.GetNumberOfElements()))

    res.PrepareForOutput()

    return res


def ComputeFeatures(inputmesh,FeatureAngle=90,skin=None):

    from BasicTools.Containers.Filters import ElementFilter
    import copy

    if skin is None:
        skinmesh = ComputeSkin(inputmesh)
        skinmeshSave = copy.deepcopy(skinmesh)
        for name,data,ids in ElementFilter(inputmesh, dimensionality = 2):
            skinmesh.elements[name].Merge(data)
    else:
        skinmesh = skin


    # we have to merge all the 2D elements form the original mesh to the skinmesh

    md = skinmesh.GetDimensionality()

    nex = skinmesh.GetNumberOfElements()

    #we use the original id to count the number of time the faces is used
    surf = {}
    for name,data in skinmesh.elements.items():
        if ElementNames.dimension[name] != md-1:
            continue
        faces = ElementNames.faces[name]
        numberOfNodes = ElementNames.numberOfNodes[name]

        ne = data.GetNumberOfElements()

        for faceType,localFaceConnectivity in faces:
            globalFaceConnectivity = data.connectivity[:,localFaceConnectivity]
            if not faceType in surf:
                surf[faceType] = {}
            surf2 = surf[faceType]
            key = np.array([nex**x for x in range(ElementNames.numberOfNodes[faceType]) ])
            for i in range(ne):
                baricentre = np.sum(skinmesh.nodes[data.connectivity[i,:] ,:],axis=0)/numberOfNodes
                cc = globalFaceConnectivity[i,:]
                lc = np.sort(cc)

                ehash = np.sum(lc*key)

                edgeVector = skinmesh.nodes[lc[0],:] - skinmesh.nodes[lc[1],:]
                planeVector = baricentre - skinmesh.nodes[lc[1],:]
                normal = np.cross(edgeVector, planeVector)
                normal /= np.linalg.norm(normal)

                if ehash in surf2:
                    surf2[ehash][0] +=1
                    normal1 = surf2[ehash][1]
                    cross = np.cross(normal, normal1)
                    angle = np.arcsin(np.linalg.norm(cross))
                    surf2[ehash][2] = 180*angle/np.pi
                else:
                    #[number of of used, normal of the first insertion,angle,   connectivity
                    surf2[ehash] = [1,normal,None,cc]
                #print(lc),
                #print(ehash),
                #print(surf2[ehash])
            #print(surf2)


    edgemesh = type(inputmesh)()
    edgemesh.nodes = inputmesh.nodes

    for name in surf:
        data = edgemesh.GetElementsOfType(name)
        surf2 = surf[name]
        for hashh,data2 in surf2.items():
            if data2[0] == 1 or data2[0] > 2:
                #print(data2)
                data.AddNewElement(data2[3],-1)
            elif data2[0] == 2 and data2[2]  >= FeatureAngle:
                data.AddNewElement(data2[3],-1)


    for eltype in [ElementNames.Bar_2, ElementNames.Bar_3]:
        bars = edgemesh.GetElementsOfType(eltype)
        bars.tags.CreateTag("Ridges").SetIds(np.arange(bars.GetNumberOfElements()))
    skinmesh.PrepareForOutput()
    edgemesh.PrepareForOutput()

    """
    Now we compute the corners

    The corner are the points:
        where 3 Ridges meet
        touched by only one Ridge
        the angle of 2 ridges is bigger than the FeatureAngle
    """

    extractedBars = edgemesh.GetElementsOfType(ElementNames.Bar_2)
    extractedBars.tighten()

    originalBars = inputmesh.GetElementsOfType(ElementNames.Bar_2)
    originalBars.tighten()

    #corners
    mask = np.zeros((edgemesh.GetNumberOfNodes()), dtype=bool )

    almanac = {}
    for bars in [originalBars,extractedBars]:
        # working form to do the operation (Please read the documentation)
        # mask[bars.connectivity.ravel()] += 1

        intmask = np.zeros((edgemesh.GetNumberOfNodes()), dtype=int )

        np.add.at(intmask, bars.connectivity.ravel(),1)


        mask[intmask > 2 ] = True
        mask[intmask == 1 ] = True

        for bar in range(bars.GetNumberOfElements()):
            id0,id1 =  bars.connectivity[bar,:]
            p0 = inputmesh.nodes[id0,:]
            p1 = inputmesh.nodes[id1,:]
            vec1 = (p1-p0)
            vec1 /= np.linalg.norm(vec1)
            for p,sign in zip([id0,id1],[1,-1]):
                if mask[p]:
                    #already in the mask no need to treat this point
                    continue

                if p in almanac:
                    vec2 = almanac[p][1]
                    angle = (180/np.pi)*np.arccos(np.dot(vec1, vec2*sign) )
                    almanac[p][2] = angle
                    almanac[p][0] +=1
                    #print(angle)
                else:
                    almanac[p] = [1,vec1*sign,id0]
                    #print(vec1)


    for p,data in almanac.items():
        if (180-data[2]) >= FeatureAngle:
            mask[p] = True


    edgemesh.nodesTags.CreateTag("Corners").AddToTag(np.where(mask)[0])

    edgemesh.PrepareForOutput()
    skinmesh.PrepareForOutput()
    return (edgemesh,skinmeshSave)

def PointToCellData(mesh,pointfield):

    if len(pointfield.shape) == 2:
        ncols = pointfield.shape[1]
    else:
        ncols  = 1

    res = np.zeros((mesh.GetNumberOfElements(),ncols),dtype=float)
    mesh.ComputeGlobalOffset()

    for name,data in mesh.elements.items():
        for i in range(ncols):
            res[data.globaloffset:data.globaloffset+data.GetNumberOfElements(),i] = (np.sum(pointfield[data.connectivity,i],axis=1)/data.connectivity.shape[1]).flatten()

    return res

def Morphing(mesh,BC,bord_tot,rayon=None,tridim=False,IDW=False):
##### method for computing the deform mesh knowing displacement of some nodes #######################################
##### BC is the known displacement in a numpy array (shape [number_of_of_known_nodes,3]) ############################
##### bord_tot contains the ids of known nodes (list of ids or boolean array)  in same order as BC ##################
##### you can choose a radius by setting rayon to a value ###########################################################
##### put tridim to True if you want to do the morphing in one step #################################################
##### put IDW to True if you want to use another morphing method without system inversion ###########################

##### https://www.researchgate.net/publication/288624175_Mesh_deformation_based_on_radial_basis_function_interpolation_computers_and_structure
    def dphi(x):
      table=x>=1
      y=-20*x*(1-x)**3
      if len(np.argwhere(table))>1:
        y[table]=np.zeros(len(np.argwhere(table)))
      return y

    def phi(x):
      table=x>=1
      y=(1-x)**4*(4*x+1)
      if len(np.argwhere(table))>1:
        y[table]=np.zeros(len(np.argwhere(table)))
      return y


    grad_max=0.8
    max_step=10
    nb_step=1
    if tridim:
      nb_step=1
    step=0
  ##################################RBF###############################################
  ####################################################################################
    nb_nodes = mesh.GetNumberOfNodes()
    new_nodes =  np.copy(mesh.GetPosOfNodes())

    if rayon==None:
      mesh.ComputeBoundingBox()
      r=np.linalg.norm(mesh.boundingMax-mesh.boundingMin)/2
    else:
      r=rayon
    if r==0:
      r=1
    while step<nb_step:
      border_nodes=new_nodes[bord_tot,:]

      rhs=BC/nb_step
      M=np.eye(np.shape(border_nodes)[0])
      #print('Building RBF operator {}/{}(shape {})'.format(step+1,nb_step,np.shape(M)))
      for j in range(np.shape(M)[0]):
          d=np.linalg.norm(border_nodes-border_nodes[j],axis=1)
          M[:,j]=phi(d/r)
      op=M
      del(M)
      #print('Solving RBF problem (shape {})'.format(np.shape(op)))
      try:
          ab=np.linalg.solve(op,rhs)
      except np.linalg.LinAlgError:
          #print('Bad conditioning of RBF operator, using least square solution')
          ab=np.linalg.lstsq(op,rhs,cond=10**(9))
      del(op)
      alpha=ab
      ds=np.zeros((nb_nodes,mesh.GetDimensionality()))
      s=np.zeros((nb_nodes,mesh.GetDimensionality()))
    #        pbar=ProgressBar()
      #print('Building RBF displacement field (shape {})'.format((nb_nodes,np.shape(border_nodes)[0])))
      for j in range(np.shape(border_nodes)[0]):
          d=np.array([np.linalg.norm(new_nodes-border_nodes[j],axis=1)])
          s=s+(phi(d/r)).T*alpha[j]
          ds=ds+(dphi(d/r)).T*alpha[j]/r
      if step==0 and tridim:
          nb_step=min(max_step,int(np.floor(np.max(np.linalg.norm(ds,axis=1)/grad_max)))+1)
          s=s/nb_step


      new_nodes += s
      step+=1

    return new_nodes

def EnsureUniquenessElements(mesh):
    cpt = 0
    for name,data in mesh.elements.items():
        dd = dict()
        for el in range(data.GetNumberOfElements()):
            n = tuple(np.sort(data.connectivity[el,:]))
            if len(n) != len(np.unique(n)):
                raise(Exception("("+str(name)+") element " + str(el) +" (global["+str(el+cpt)+"])" + " use a point more than once ("+str(data.connectivity[el,:])+")"  ))
            if n in dd.keys():
                raise(Exception("("+str(name)+") element " + str(el) +" (global["+str(el+cpt)+"])" + " is a duplication of element " + str(dd[n]) +" (global["+str(dd[n]+cpt)+"])" ))
            dd[n] = el
        cpt = data.GetNumberOfElements()

def MeshQualityAspectRatioBeta(mesh):
    #https://cubit.sandia.gov/public/15.2/help_manual/WebHelp/mesh_generation/mesh_quality_assessment/tetrahedral_metrics.htm
    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] == 0:
            pass
        elif ElementNames.dimension[name] == 1:
            pass
        elif ElementNames.dimension[name] == 2:
            pass
        elif ElementNames.dimension[name] == 3:
            if name == ElementNames.Tetrahedron_4:
               mmax = 0
               for el in range(data.GetNumberOfElements()):
                   n = data.connectivity[el,:]
                   nodes = mesh.nodes[n,:]
                   p0 = nodes[0,:]
                   p1 = nodes[1,:]
                   p2 = nodes[2,:]
                   p3 = nodes[3,:]
                   #print("---")

                   a = np.linalg.norm(p1-p0)
                   b = np.linalg.norm(p2-p0)
                   #print(a)
                   #print(b)
                   normal = np.cross(p1-p0,p2-p0)
                   #https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
                   base_area = 0.5*a*b*math.sqrt(1-(np.dot(p1-p0,p2-p0)/(a*b))**2)

                   normal /= np.linalg.norm(normal)

                   # distance to the oposite point
                   d = np.dot(normal,p3-p0)

                   volume = base_area*d
                   #print(base_area)
                   #print(d)
                   inscribed_sphere_radius = d/3
                   #https://math.stackexchange.com/questions/2820212/circumradius-of-a-tetrahedron

                   c = np.linalg.norm(p3-p0)

                   A = np.linalg.norm(p3-p2)
                   B = np.linalg.norm(p1-p3)
                   C = np.linalg.norm(p2-p1)

                   #print("-")
                   #print((a*A+b*B+c*C)*
                   #                                      (a*A+b*B-c*C)*
                   #                                      (a*A-b*B+c*C)*
                   #                                      (-a*A+b*B+c*C))
                   circumradius_sphere_radius = math.sqrt((a*A+b*B+c*C)*
                                                         (a*A+b*B-c*C)*
                                                         (a*A-b*B+c*C)*
                                                         (-a*A+b*B+c*C))/(24*volume)


                   AspectRatioBeta = circumradius_sphere_radius/(3.0 * inscribed_sphere_radius)
                   #print(volume)
                   #print(circumradius_sphere_radius)
                   #print(inscribed_sphere_radius)
                   #print(circumradius_sphere_radius)
                   #print(AspectRatioBeta)
                   mmax = max([mmax,AspectRatioBeta])
                   if AspectRatioBeta > 1000:
                       raise(Exception("Element " +str(el) + " has quality of " +str(AspectRatioBeta)) )
               #print("mmax")
               #print(mmax)

        else:
            raise
###############################################################################
def CheckIntegrity_CreateUniformMeshOfBars(GUI=False):
    print(CreateUniformMeshOfBars(0,8,10))
    return "ok"


def CheckIntegrity_EnsureUniquenessElements(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    # this function must work
    EnsureUniquenessElements(mesh)

    tets = [[0,1,2,3],[3,1,0,2]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    # This function must not work
    try:
       EnsureUniquenessElements(mesh)
       return "No ok"
    except :
       pass
    return "ok"


def CheckIntegrity_CreateCube(GUI = False):
    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=False)
    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)
    return "ok"


def CheckIntegrity_CreateSquare(GUI = False):
    mesh = CreateSquare(dimensions=[20,21],spacing=[2.,2.],ofTris=False)
    mesh = CreateSquare(dimensions=[20,21],spacing=[2.,2.],ofTris=True)
    return "ok"


def CheckIntegrity_Morphing(GUI = False):

   mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)
   BC = np.empty((mesh.GetNumberOfNodes(),3),dtype=float)
   bord_tot = np.empty(mesh.GetNumberOfNodes(),dtype=int)
   cpt = 0
   print(mesh)
   for name,data in mesh.elements.items():

       if ElementNames.dimension[name] != 2:
           continue

       ids = data.GetNodesIdFor(data.GetTag("X0").GetIds())
       print(ids)
       bord_tot[cpt:cpt+len(ids)] = ids
       BC[cpt:cpt+len(ids),:] = 0
       cpt += len(ids)

       ids = data.GetNodesIdFor(data.GetTag("X1").GetIds())
       print(ids)
       bord_tot[cpt:cpt+len(ids)] = ids
       BC[cpt:cpt+len(ids),:] = [[0,0,10]]
       cpt += len(ids)

   BC = BC[0:cpt,:]
   bord_tot = bord_tot[0:cpt]

   new_p1 = Morphing(mesh, BC,bord_tot)
   new_p2 = Morphing(mesh, BC,bord_tot,rayon= 20. )

   new_p0 = np.copy(mesh.nodes)
   new_p0[bord_tot,:] += BC
   mesh.nodeFields["morph0"] = new_p0
   mesh.nodeFields["morph1"] = new_p1
   mesh.nodeFields["morph2"] = new_p2


   if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

   return "ok"

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


def CheckIntegrity_ComputeFeatures(GUI =False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([2,3,4]);
    myMesh.SetOrigin([-1.0,-1.0,-1.0]);
    myMesh.SetSpacing([2., 2.,2]/myMesh.GetDimensions());
    print(myMesh)
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=True)
    print(res2)


    edges,skin = ComputeFeatures(res2,FeatureAngle=80)
    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=edges,filename="edges.xmf")
        OpenInParaView(mesh=skin,filename="skin.xmf")

        for name,data in edges.elements.items():
              res2.GetElementsOfType(name).Merge(data)

        OpenInParaView(res2,filename="all+edges.xmf")
        print(res2)

    print(edges)
    return "ok"

def CheckIntegrity_ComputeSkin(GUI=False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([2,3,4]);
    myMesh.SetOrigin([-1.0,-1.0,-1.0]);
    myMesh.SetSpacing([2., 2.,2]/myMesh.GetDimensions());
    print(myMesh)
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=True)
    print(res2)

    skin = ComputeSkin(res2)
    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView

        OpenInParaView(skin,filename="skin.xmf")

    print(skin)
    return "ok"

def CheckIntegrity_ExtractElementsByImplicitZone(GUI=False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([20,30,40]);
    myMesh.SetOrigin([-1.0,-1.0,-1.0]);
    myMesh.SetSpacing([2., 2.,2]/myMesh.GetDimensions());
    print(myMesh)
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=False)
    print(res2)

    class OPSphere(object):
        def __init__(self):
            self.center = np.array([0.0,0.0,0.0],dtype=np.float)
            self.radius = 0.5
        def __call__(self,pos):
            res = np.sqrt(np.sum((pos-self.center)**2,axis=1))-self.radius
            return res

    myOp = OPSphere()

    res = ExtractElementsByImplicitZone(res2,myOp)
    print(res)

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    WriteMeshToXdmf(tempdir+"Test_ExtractElementsByImplicitZone.xdmf",res,PointFields=[res.originalIDNodes],PointFieldsNames=["originalIDNodes"] )
    print(tempdir)
    return "OK"

def CheckIntegrity_DeleteInternalFaces(GUI=False):

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],[3,0,2,4]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    ## add 2 tris
    tris = mesh.GetElementsOfType(ElementNames.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([3,0,2],1)
    tris.AddNewElement([3,0,2],2)
    #print(mesh)

    DeleteInternalFaces(mesh)
    print(mesh)
    return "ok"


def LowerNodesDimension(mesh):
    newDim = mesh.GetDimensionality() - 1
    newNodes = np.empty((mesh.GetNumberOfNodes(),newDim))
    for i in range(mesh.GetNumberOfNodes()):
        newNodes[i,:] = mesh.nodes[i,0:newDim]
    mesh.nodes = newNodes
    return mesh


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
                


def CheckIntegrity_GetValueAtPosLinearSymplecticMesh(GUI=False):
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOf
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

def CheckIntegrity_CreateMeshOfTriangles(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    print(res)
    return "OK"

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

def CheckIntegrity_QuadToLin(GUI=False):
    myMesh = UnstructuredMesh()
    myMesh.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0.5,0,0],[0.5,0.5,0],[0,0.5,0],[0,0,0.5],[0.5,0,0.5],[0,0.5,0.5]] ,dtype=np.float)
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

    print("-----")
    print(QuadToLin(myMesh,divideQuadElements=False))
    print(QuadToLin(myMesh,divideQuadElements=True))
    print(QuadToLin(myMesh,divideQuadElements=True,lineariseMiddlePoints=True))
    return "ok"

def CheckIntegrity_CleanDoubleNodes(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh)
    CleanDoubleNodes(mesh,tol =0)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh, nodesToTestMask= np.array([True,False,True,False,True],dtype=np.bool) )
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover


    return "ok"

def CheckIntegrity_CleanLonelyNodes(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanLonelyNodes(mesh)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_GetVolume(GUI=False):
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    vol = GetVolume(mesh)
    if vol != (1./6.):
        raise Exception('Error en the calculation of the volumen')# pragma: no cover

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([3,3,3]);
    myMesh.SetSpacing([0.5, 0.5,0.5]);
    print(myMesh)
    vol  = GetVolume(CreateMeshFromConstantRectilinearMesh(myMesh))
    if abs(vol-1.) > 1e-8 :
        raise Exception('Error en the calculation of the volumen')# pragma: no cover

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

def CheckIntegrity_ExtractElementByTags(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.AddElementToTagUsingOriginalId(0,"first")
    res.AddElementToTagUsingOriginalId(1,"second")
    res.AddElementToTagUsingOriginalId(0,"all")
    res.AddElementToTagUsingOriginalId(1,"all")
    res.GetNodalTag("Point3").AddToTag(3)

    if ExtractElementByTags(res,["first"] ).GetNumberOfElements() != 1:
        raise# pragma: no cover
    if ExtractElementByTags(res,["Point3"] ).GetNumberOfElements() != 1:
        raise# pragma: no cover
    if ExtractElementByTags(res,["Point3"],allNodes=True ).GetNumberOfElements() != 0:
        raise# pragma: no cover

    if ExtractElementByTags(res,[],dimensionalityFilter=2 ).GetNumberOfElements() != 2:
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_CleanEmptyTags(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetNodalTag("Point0").AddToTag(0)
    res.GetNodalTag("EmptyTagN")

    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("EmptyTagE")
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("Element0").AddToTag(0)

    CleanEmptyTags(res)
    if len(res.nodesTags) != 1 :
        raise# pragma: no cover

    if "EmptyTagE" in  res.GetNamesOfElemTags() :
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_AddTagPerBody(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    resII = MirrorMesh(res,x=0,y=0,z=0)
    print( resII.nodes)
    print( resII.GetElementsOfType(ElementNames.Triangle_3))
    print(AddTagPerBody(resII))

    if len(resII.nodesTags) != 8:
        raise # pragma: no cover
    print(resII.nodesTags)

    return "ok"

def CheckIntegrity_ExtractElementsByMask(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("tri1").AddToTag(0)

    #tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([False,True]))
    #print(tri.connectivity)
    print(res.GetElementsOfType(ElementNames.Triangle_3).connectivity)

    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0],dtype=np.int))
    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0,1],dtype=np.bool))

    print(tri.connectivity)
    print(tri.originalIds)
    return "ok"

def CheckIntegrity_GetDualGraph(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    dg, nused = GetDualGraph(res)

    return "ok"



def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_GetValueAtPosLinearSymplecticMesh,
    CheckIntegrity_CreateUniformMeshOfBars,
    CheckIntegrity_CreateCube,
    CheckIntegrity_CreateSquare,
    CheckIntegrity_EnsureUniquenessElements,
    CheckIntegrity_Morphing,
    CheckIntegrity_PointToCellData,
    CheckIntegrity_ComputeFeatures,
    CheckIntegrity_ComputeSkin,
    CheckIntegrity_ExtractElementsByImplicitZone,
    CheckIntegrity_DeleteInternalFaces,
    CheckIntegrity_ExtractElementsByMask,
    CheckIntegrity_GetVolume,
    CheckIntegrity_CreateMeshOfTriangles,
    CheckIntegrity_CreateMeshFromConstantRectilinearMesh,
    CheckIntegrity_QuadToLin,
    CheckIntegrity_CleanDoubleNodes,
    CheckIntegrity_CleanLonelyNodes,
    CheckIntegrity_MirrorMesh,
    CheckIntegrity_ExtractElementByTags,
    CheckIntegrity_CleanEmptyTags,
    CheckIntegrity_AddTagPerBody,
    CheckIntegrity_GetDualGraph,
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
