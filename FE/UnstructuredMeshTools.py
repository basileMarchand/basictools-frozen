# -*- coding: utf-8 -*-
import numpy as np

from BasicTools.FE.UnstructuredMesh import UnstructuredMesh
import BasicTools.FE.ElementNames as ElementNames
import BasicTools.Helpers.BaseOutputObject as BaseOutputObject


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
    return res

def CreateMeshFromConstantRectilinearMesh(CRM, ofTetras= False,out=None):
    if out is None:
        res = UnstructuredMesh()
    else:
        res = out # pragma: no cover

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
        if ofTetras:# pragma: no cover
            raise

        elementtype = ElementNames.Quadrangle_4

    elements = res.GetElementsOfType(elementtype)
    elements.connectivity = np.zeros((nbelements,ElementNames.numberOfNodes[elementtype]),dtype=np.int)
    elements.cpt = nbelements

    if ofTetras:
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
        CRM.GenerateFullConnectivity()
        elements.connectivity  =  CRM.connectivity

    elements.originalIds = np.arange(0,elements.GetNumberOfElements(),dtype=np.int)

    return res

def QuadToLin(inputmesh, divideQuadElements=True,lineariseMiddlePoints=False):

    res = type(inputmesh)()
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

    return res

def CleanDoubleNodes(res, tol = None, nodesToTestMask= None):

    BaseOutputObject.BaseOutputObject().PrintDebug("in CleanDoubleNodes")

    if tol is None:
        res.ComputeBoundingBox()
        tol = np.linalg.norm(res.boundingMax - res.boundingMin)*1e-7

    nbnodes = res.GetNumberOfNodes()
    toKeep = np.zeros(nbnodes, dtype=np.bool )
    newindex = np.zeros(nbnodes, dtype=np.int )
    tol = tol**2
    def dist2(array,value):
        x = array-value
        return x[0]**2 +x[1]**2+x[2]**2

    if nodesToTestMask is None:
        cpt =0
        for i in range(nbnodes):
            posi = res.nodes[i,:]
            for j in range(i):
                #dist =np.linalg.norm(res.nodes[i,:]-res.nodes[j,:])
                dist = dist2(posi,res.nodes[j,:] )
                if dist < tol and toKeep[j]:
                    newindex[i] = newindex[j]
                    break
            else:
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
                    if dist < tol :
                        newindex[i] = newindex[j]
                        break
            else:
                newindex[i] = cpt
                cpt += 1
                toKeep[i] = True;





    res.nodes = res.nodes[toKeep,:]

    for elementName in res.elements:
        elements = res.elements[elementName]
        elements.connectivity = newindex[elements.connectivity]

    BaseOutputObject.BaseOutputObject().PrintDebug("CleanDoubleNodes Done")

def CleanLonelyNodes(res):

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
    res.nodes = res.nodes[usedNodes ,:]
    res.originalIDNodes = originalIDNodes[0:cpt]

    #node tags
    for tag in res.nodesTags :
        tag.SetIds(NewIndex[np.extract(usedNodes[tag.GetIds()],tag.GetIds() )])

    #renumbering the connectivity matrix
    for elementName in res.elements:
        elements = res.elements[elementName]
        elements.connectivity = NewIndex[elements.connectivity]

def MirrorMesh(inmesh,x=None,y=None,z=None) :
    nbpoints = inmesh.GetNumberOfNodes()

    outmesh = type(inmesh)()
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


    for tag in inelems.tags  :
       temp = np.extract(mask[tag.GetIds()],tag.GetIds())
       newid = newIndex[temp]
       outelems.tags.CreateTag(tag.name).SetIds(newid)

    return outelems

def ExtractElementByTags(inmesh,tagsToKeep, allNodes=False,dimensionalityFilter= None):

    outmesh = type(inmesh)()

    outmesh.nodes = np.copy(inmesh.nodes)
    outmesh.originalIDNodes = np.copy(inmesh.originalIDNodes)

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

       for tag in elems.tags  :
           #print(tag.id)
           temp = np.extract(toKeep[tag.GetIds()],tag.GetIds())
           #print(temp)
           newid = newIndex[temp]
           #print(newid)
           outelem.tags.CreateTag(tag.name,errorIfAlreadyCreated=False).SetIds(newid)


    CleanLonelyNodes(outmesh)
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

def GetDualGraph(inmesh, maxNumConnections=100):

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


def MeshToVtk(mesh, vtkobject=None, TagsAsFields=False):


    # From www.vtk;org/wp-content/updloads/2015/04/file-formats.pdf

    vtknumbers = {}

    vtknumbers[ElementNames.Point_1] = 1

    vtknumbers[ElementNames.Bar_2] = 3

    vtknumbers[ElementNames.Triangle_3] = 5
    vtknumbers[ElementNames.Quadrangle_4] = 9
    vtknumbers[ElementNames.Tetrahedron_4] = 10

    vtknumbers[ElementNames.Hexaedron_8] = 12
    vtknumbers[ElementNames.Wedge_6] = 13
    vtknumbers[ElementNames.Pyramid_5] = 14

    vtknumbers[ElementNames.Bar_3] = 21
    vtknumbers[ElementNames.Triangle_6] = 22
    vtknumbers[ElementNames.Quadrangle_8] = 23
    vtknumbers[ElementNames.Tetrahedron_10] = 24
    vtknumbers[ElementNames.Hexaedron_20] = 25

    try:
        from paraview.vtk import vtkPolyData as vtkPolyData
        from paraview.vtk import vtkUnstructuredGrid as vtkUnstructuredGrid
        from paraview.vtk import vtkPoints
        from paraview.vtk import vtkFloatArray
        from paraview.vtk import vtkIntArray
        from paraview.vtk import vtkIdList
    except :
        from vtk import vtkPolyData
        from vtk import vtkUnstructuredGrid
        from vtk import vtkPoints
        from vtk import vtkFloatArray
        from vtk import vtkIntArray
        from vtk import vtkIdList

    if vtkobject is None:


        usePoly = True
        for  elementsname,elementContainer in mesh.elements.items():
            if ElementNames.dimension[elementsname] == 3:
                usePoly = False
                break
        if usePoly:
            output = vtkPolyData()
        else:
            output = vtkUnstructuredGrid()

    else:
        output = vtkobject # pragma: no cover



    output.Allocate(mesh.GetNumberOfElements())
    ##copy points
    pts = vtkPoints()
    pts.Allocate(mesh.GetNumberOfNodes())
    if mesh.nodes.shape[1] == 3 :
        for p in range(mesh.GetNumberOfNodes()):
            point = mesh.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],point[2])
    else:
        #2DCase
        for p in range(mesh.GetNumberOfNodes()):
            point = mesh.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],0.0)

    output.SetPoints(pts)

    if hasattr(mesh,"nodeFields"):
        for name,data in mesh.nodeFields.items():
            #VTK_data = numpy_support.numpy_to_vtk(num_array=np.swapaxes(phi,0,2).ravel(), deep=True, array_type=vtk.VTK_FLOAT)
            #VTK_data.SetName(name)

            pd = vtkFloatArray()
            pd.SetName(name)
            if len(data.shape) == 1:
                pd.SetNumberOfComponents(1)
            else:
                pd.SetNumberOfComponents(data.shape[1])
            pd.SetNumberOfTuples(mesh.GetNumberOfNodes())

            if len(data.shape) > 1:
              cpt = 0
              for i in range(mesh.GetNumberOfNodes()):
                 for j in range(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1
              output.GetPointData().AddArray(pd)
            else:
              cpt = 0
              for i in range(mesh.GetNumberOfNodes()):
                    pd.SetValue(cpt, data[i])
                    cpt +=1
              output.GetPointData().AddArray(pd)

    if TagsAsFields:
        for tag in mesh.nodesTags:
            pd = vtkIntArray()
            pd.SetName(tag.name)
            pd.SetNumberOfComponents(1)
            pd.SetNumberOfTuples(mesh.GetNumberOfNodes())
            pd.FillComponent(0,0);

            for i in tag.GetIds():
                pd.SetValue(i,1 )
            output.GetPointData().AddArray(pd)


    for elementsname,elementContainer in mesh.elements.items():
        pointIds = vtkIdList()
        npe = elementContainer.GetNumberOfNodesPerElement()
        pointIds.SetNumberOfIds(npe)
        vtknumber = vtknumbers[elementsname]
        for e in range(elementContainer.GetNumberOfElements()):
            for i in range(npe):
                pointIds.SetId(i,elementContainer.connectivity[e,i])
            output.InsertNextCell(vtknumber, pointIds)

    if hasattr(mesh,"elemFields"):
        for name,data in mesh.elemFields.items():

            pd = vtkFloatArray()
            pd.SetName(name)

            if len(data.shape) == 1:
                pd.SetNumberOfComponents(1)
            else:
                pd.SetNumberOfComponents(data.shape[1])

            pd.SetNumberOfTuples(mesh.GetNumberOfElements())

            if len(data.shape) > 1:
              cpt = 0
              for i in range(mesh.GetNumberOfElements()):
                 for j in range(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1
            else:
              cpt = 0
              for i in range(mesh.GetNumberOfElements()):
                    pd.SetValue(cpt, data[i])
                    cpt +=1

            output.GetCellData().AddArray(pd)
            #output.GetCellData().SetScalars(pd)

    if TagsAsFields:
        elementTags = mesh.GetNamesOfElemTags()
        for tagname in elementTags:
            ids = mesh.GetElementsInTag(tagname)
            pd = vtkIntArray()
            pd.SetName(tagname)
            pd.SetNumberOfComponents(1)
            pd.SetNumberOfTuples(mesh.GetNumberOfElements())
            pd.FillComponent(0,0);
            for i in ids:
                pd.SetValue(i,1 )
            output.GetCellData().AddArray(pd)

    return output

def VtkToMesh(vtkmesh, meshobject=None, TagsAsFields=False):

    if meshobject is None:
        out = UnstructuredMesh()
    else:

        out = meshobject

    from vtk.util import numpy_support
    data = vtkmesh.GetPoints().GetData()
    out.nodes = numpy_support.vtk_to_numpy(data)

    nc = vtkmesh.GetNumberOfCells()

    vtknumbers = {}
    vtknumbers[3 ] = ElementNames.Bar_2
    vtknumbers[4 ] = ElementNames.Bar_2
    vtknumbers[5 ] = ElementNames.Triangle_3
    vtknumbers[9 ] = ElementNames.Quadrangle_4
    vtknumbers[12] = ElementNames.Hexaedron_8
    vtknumbers[25] = ElementNames.Hexaedron_20
    vtknumbers[21] = ElementNames.Bar_3
    vtknumbers[22] = ElementNames.Triangle_6
    vtknumbers[10] = ElementNames.Tetrahedron_4
    vtknumbers[24] = ElementNames.Tetrahedron_10

    for i in range(nc):
        cell= vtkmesh.GetCell(i)
        ct = cell.GetCellType()
        et = vtknumbers[ct]
        np = cell.GetNumberOfPoints()
        #polyline case
        # we have to be careful because we potentialy change the number of
        # elements in the mesh if we have polylines
        if ct ==4:
            for i in range(np-1):
                out.GetElementsOfType(et).AddNewElement([cell.GetPointId(i),cell.GetPointId(i+1) ] ,i)
        else:
            out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j) for j in range(np)] ,i)
    return out

def DeleteInternalFaces(mesh):
    elems2D = {}
    elemsMask = {}
    usedPointsBy2DElements = []
    mesh.PrepareForOutput()
    #generate a list of all the 2D elements and a empty mask
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
                dt  = [ ('col'+str(x),np.int) for x in xrange(faces.shape[1]) ]

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



###############################################################################
def CheckIntegrity_DeleteInternalFaces():

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




def CheckIntegrity_VtkToMesh():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.nodeFields = {"x": res.nodes[:,0], "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"firstPoint": res.GetElementsOfType(ElementNames.Triangle_3).connectivity[:,0], "conn": res.GetElementsOfType(ElementNames.Triangle_3).connectivity }
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("FirstTriangle").AddToTag(0)
    sol = MeshToVtk(res,TagsAsFields= True)

    print(res)
    print(VtkToMesh(sol))
    return 'ok'

def CheckIntegrity_MeshToVtk():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.nodeFields = {"x": res.nodes[:,0], "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"firstPoint": res.GetElementsOfType(ElementNames.Triangle_3).connectivity[:,0], "conn": res.GetElementsOfType(ElementNames.Triangle_3).connectivity }
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("FirstTriangle").AddToTag(0)
    sol = MeshToVtk(res,TagsAsFields= True)

    res = CreateMeshOfTriangles([[0,0],[1,0],[0,1],[1,1] ], [[0,1,2],[2,1,3]])

    sol = MeshToVtk(res )
    print(sol)
    return "OK"


def CheckIntegrity_CreateMeshOfTriangles():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    print(res)
    return "OK"

def CheckIntegrity_CreateMeshFromConstantRectilinearMesh():
    from BasicTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh()
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

def CheckIntegrity_QuadToLin():
    myMesh = UnstructuredMesh()
    myMesh.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0.5,0,0],[0.5,0.5,0],[0,0.5,0],[0,0,0.5],[0.5,0,0.5],[0,0.5,0.5]] ,dtype=np.float)
    tag = myMesh.GetNodalTag("linPoints")
    tag.AddToTag(0)
    tag.AddToTag(1)
    tag.AddToTag(2)
    tag.AddToTag(3)
    import BasicTools.FE.ElementNames as ElementNames

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

def CheckIntegrity_CleanDoubleNodes():
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh, nodesToTestMask= np.array([True,False,True,False,True],dtype=np.bool) )
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover


    return "ok"

def CheckIntegrity_CleanLonelyNodes():
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanLonelyNodes(mesh)
    if mesh.GetNumberOfNodes() != 4:
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_GetVolume():
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    vol = GetVolume(mesh)
    if vol != (1./6.):
        raise Exception('Error en the calculation of the volumen')# pragma: no cover

    from BasicTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([3,3,3]);
    myMesh.SetSpacing([0.5, 0.5,0.5]);
    print(myMesh)
    vol  = GetVolume(CreateMeshFromConstantRectilinearMesh(myMesh))
    if abs(vol-1.) > 1e-8 :
        raise Exception('Error en the calculation of the volumen')# pragma: no cover

    return "ok"

def CheckIntegrity_MirrorMesh():
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

def CheckIntegrity_ExtractElementByTags():
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


def CheckIntegrity_CleanEmptyTags():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetNodalTag("Point0").AddToTag(0)
    res.GetNodalTag("EmptyTagN")

    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("EmptyTagE")
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("Element0").AddToTag(0)

    CleanEmptyTags(res)
    if len(res.nodesTags) != 1 :
        raise# pragma: no cover

    if len(res.GetNamesOfElemTags()) != 1 :
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_AddTagPerBody():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    resII = MirrorMesh(res,x=0,y=0,z=0)
    print( resII.nodes)
    print( resII.GetElementsOfType(ElementNames.Triangle_3))
    print(AddTagPerBody(resII))

    if len(resII.nodesTags) != 8:
        raise # pragma: no cover
    print(resII.nodesTags)

    return "ok"

def CheckIntegrity_ExtractElementsByMask():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("tri1").AddToTag(0)

    #tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([False,True]))
    #print(tri.connectivity)

    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0],dtype=np.int))
    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0,1],dtype=np.bool))

    print(tri.connectivity)
    return "ok"

def CheckIntegrity_GetDualGraph():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    dg, nused = GetDualGraph(res)

    return "ok"


def CheckIntegrity():


    CheckIntegrity_DeleteInternalFaces()
    CheckIntegrity_ExtractElementsByMask()
    CheckIntegrity_GetVolume()
    CheckIntegrity_CreateMeshOfTriangles()
    CheckIntegrity_CreateMeshFromConstantRectilinearMesh()
    CheckIntegrity_QuadToLin()
    CheckIntegrity_CleanDoubleNodes()
    CheckIntegrity_CleanLonelyNodes()
    CheckIntegrity_MirrorMesh()
    CheckIntegrity_ExtractElementByTags()
    CheckIntegrity_CleanEmptyTags()
    CheckIntegrity_AddTagPerBody()
    CheckIntegrity_MeshToVtk()
    CheckIntegrity_VtkToMesh()
    CheckIntegrity_GetDualGraph()
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
