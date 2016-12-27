# -*- coding: utf-8 -*-
import numpy as np

from OTTools.FE.UnstructuredMesh import UnstructuredMesh
import OTTools.FE.ElementNames as ElementNames



def CreateMeshOfTriangles(points,tris):
    
    res = UnstructuredMesh()
    
    res.nodes = np.array(points, dtype=np.double)
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int)

    elements = res.GetElementsOfType(ElementNames.Triangle_3)
    elements.connectivity = np.array(tris,dtype=np.int)
    elements.originalIds = np.arange(0,elements.connectivity.shape[0],dtype=np.int)
    elements.cpt = elements.connectivity.shape[0]
    return res
    
def CreateMeshFromConstantRectilinearMesh(CRM, ofTetras= False):
    res = UnstructuredMesh()
    res.nodes = CRM.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int);
    
    
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

        for elem in xrange(CRM.GetNumberOfElements()):
            
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
   
    res = UnstructuredMesh()
    res.nodes = inputmesh.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int);
    import copy 
    res.nodesTags = copy.deepcopy(inputmesh.nodesTags)
    
    for elementName in inputmesh.elements.keys():
        quadElement = inputmesh.elements[elementName]
        if elementName == ElementNames.Tetrahedron_10:
            
            lineelements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 8
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*8)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*8
                for i in xrange(quadElement.GetNumberOfElements()):
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
                for i in xrange(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2,3]];
                
                
        elif elementName == ElementNames.Triangle_6:
            
            lineelements = res.GetElementsOfType(ElementNames.Triangle_3)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 4
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*4)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*4
                for i in xrange(quadElement.GetNumberOfElements()):
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
                for i in xrange(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1,2]];
    
        elif elementName == ElementNames.Bar_3:
            
            lineelements = res.GetElementsOfType(ElementNames.Bar_2)
            initNbElem = lineelements.GetNumberOfElements();
            if divideQuadElements:
                nbOfNewElements = 2
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements()*2)
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()*2
                for i in xrange(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i*2+0,:] = quadConn[[0,2]];
                    lineelements.connectivity[initNbElem+i*2+1,:] = quadConn[[2,1]];
                    if lineariseMiddlePoints :
                        res.nodes[quadConn[2],:] = (res.nodes[quadConn[0],:] + res.nodes[quadConn[1],:] )/2
            else:
                nbOfNewElements = 1
                lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
                lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
                for i in xrange(quadElement.GetNumberOfElements()):
                    quadConn = quadElement.connectivity[i,:];
                    lineelements.connectivity[initNbElem+i,:] = quadConn[[0,1]];                
        elif elementName == ElementNames.Bar_2:
            lineelements = res.GetElementsOfType(ElementNames.Bar_2)
            initNbElem = lineelements.GetNumberOfElements();
            
            lineelements.Reserve(initNbElem+quadElement.GetNumberOfElements())
            nbOfNewElements = 1
            lineelements.connectivity[initNbElem:initNbElem+quadElement.GetNumberOfElements(),:] = quadElement.connectivity
            lineelements.cpt = initNbElem+quadElement.GetNumberOfElements()
        else:
            raise Exception('Error : not coded yet for this type of elements ' + str(elementName))# pragma: no cover
        #copy of tags 
        for originaltagName in quadElement.tags :
            destinationtag = lineelements.GetTag(originaltagName)
            originaltag = quadElement.tags[originaltagName]
            for i in xrange(originaltag.cpt):
                for t in xrange(nbOfNewElements):
                    destinationtag.AddToTag(initNbElem+originaltag.id[i]*nbOfNewElements+t)
                    
            destinationtag.tighten()
            
        res.ComputeGlobalOffset()            


    if divideQuadElements == False:
        CleanLonelyNodes(res)      
    
    return res
    
    
def CleanLonelyNodes(res):
    
    usedNodes = np.zeros(res.GetNumberOfNodes(),dtype=np.bool )
    NewIndex =  np.zeros(res.GetNumberOfNodes(),dtype=np.int )
    for elementName in res.elements.keys():
        elements = res.elements[elementName]
        usedNodes[elements.connectivity.flatten()] = True;

    cpt = 0 ;
    for n in xrange(res.GetNumberOfNodes()):
        NewIndex[n] = cpt
        cpt += usedNodes[n]
        
    #filter the nodes
    res.nodes = res.nodes[usedNodes ,:]
    res.originalIDNodes = res.originalIDNodes[usedNodes ]
    
    #node tags
    for tagName in res.nodesTags :
        tag = res.nodesTags[tagName]
        tag.tighten()
        tag.SetIds(NewIndex[np.extract(usedNodes[tag.id],tag.id )])
    
    #renumbering the connectivity matrix
    for elementName in res.elements.keys():
        elements = res.elements[elementName]
        elements.connectivity = NewIndex[elements.connectivity]
    
    
def CheckIntegrity():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    CleanLonelyNodes(res)
    
    ###########################
    from OTTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh

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
    
    
    myMesh = UnstructuredMesh()
    myMesh.nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0.5,0,0],[0.5,0.5,0],[0,0.5,0],[0,0,0.5],[0.5,0,0.5],[0,0.5,0.5]] ,dtype=np.float)
    tag = myMesh.GetNodalTag("linPoints")
    tag.AddToTag(0)
    tag.AddToTag(1)
    tag.AddToTag(2)
    tag.AddToTag(3)
    import OTTools.FE.ElementNames as ElementNames

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
    
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover