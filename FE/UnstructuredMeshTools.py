# -*- coding: utf-8 -*-
import numpy as np

from OTTools.FE.UnstructuredMesh import UnstructuredMesh
import OTTools.FE.ElementNames as ElementNames



def CreateMeshOfTriangles(points,tris):
    return CreateMeshOf(points,tris,elemName = ElementNames.Triangle_3 )

def CreateMeshOf(points,connectivity,elemName = None):

    if elemName is None:
        raise Exception("Need a element name ")# pragma: no cover

    res = UnstructuredMesh()

    res.nodes = np.array(points, dtype=np.double)
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int)

    elements = res.GetElementsOfType(elemName)
    elements.connectivity = np.array(connectivity,dtype=np.int)
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
        for originaltag in quadElement.tags :
            destinationtag = lineelements.GetTag(originaltag.name)
            for i in xrange(originaltag.cpt):
                for t in xrange(nbOfNewElements):
                    destinationtag.AddToTag(initNbElem+originaltag.id[i]*nbOfNewElements+t)

            destinationtag.tighten()

        res.ComputeGlobalOffset()


    if divideQuadElements == False:
        CleanLonelyNodes(res)

    return res

def CleanDoubleNodes(res, tol = None):
    if tol is None:
        res.ComputeBoundingBox()
        tol = np.linalg.norm(res.boundingMax - res.boundingMin)*1e-7

    nbnodes = res.GetNumberOfNodes()
    toKeep = np.zeros(nbnodes, dtype=np.bool )
    newindex = np.zeros(nbnodes, dtype=np.int )

    cpt =0
    for i in xrange(nbnodes):
        for j in xrange(i):
            dist =np.linalg.norm(res.nodes[i,:]-res.nodes[j,:])
            if dist < tol and toKeep[j]:
                newindex[i] = newindex[j]
                break
        else:
            newindex[i] = cpt
            cpt += 1
            toKeep[i] = True;

    res.nodes = res.nodes[toKeep,:]

    for elementName in res.elements.keys():
        elements = res.elements[elementName]
        elements.connectivity = newindex[elements.connectivity]

def CleanLonelyNodes(res):

    usedNodes = np.zeros(res.GetNumberOfNodes(),dtype=np.bool )
    for elementName in res.elements.keys():
        elements = res.elements[elementName]
        usedNodes[elements.connectivity.flatten()] = True;

    cpt = 0 ;
    NewIndex =  np.zeros(res.GetNumberOfNodes(),dtype=np.int )-1
    for n in xrange(res.GetNumberOfNodes()):
        if usedNodes[n]:
            NewIndex[n] = cpt
            cpt += 1

    #filter the nodes
    res.nodes = res.nodes[usedNodes ,:]
    #print(res.originalIDNodes.shape)
    #print(usedNodes.shape)
#    res.originalIDNodes = res.originalIDNodes[usedNodes ]

    #node tags
    for tag in res.nodesTags :
        tag.tighten()
        tag.SetIds(NewIndex[np.extract(usedNodes[tag.id],tag.id )])

    #renumbering the connectivity matrix
    for elementName in res.elements.keys():
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
    cpt = nbpoints

    if x is not None:
        vec = np.array([ [  -1,1,1], ],dtype=np.float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[x,0,0])*vec+ [x,0,0]
        cpt = cpt*2

    if y is not None:
        vec = np.array([ [  1,-1,1], ],dtype=np.float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[0,y,0])*vec+ [0,y,0]
        cpt = cpt*2

    if z is not None:
        vec = np.array([ [  1,1,-1], ],dtype=np.float)
        outmesh.nodes[cpt:(2*cpt),:] = (outmesh.nodes[0:cpt,:]-[0,0,z])*vec+ [0,0,z]
        cpt = cpt*2

    for name,vals in inmesh.elements.iteritems():
        nbelements = vals.GetNumberOfElements()
        outelements = outmesh.GetElementsOfType(name)
        outelements.Reserve(nbelements*(2**d))
        outelements.connectivity[0:nbelements,:] = vals.connectivity
        cpt = nbelements
        pcpt = nbpoints
        permutation = ElementNames.mirrorPermutation[name]
        if x is not None:
            outelements.connectivity[cpt:(2*cpt),:] = (outelements.connectivity[0:cpt,:]+pcpt)[:,permutation]
            pcpt = pcpt *2
            cpt = cpt*2

        if y is not None:
            outelements.connectivity[cpt:(2*cpt),:] = (outelements.connectivity[0:cpt,:]+pcpt)[:,permutation]
            pcpt = pcpt *2
            cpt = cpt*2

        if z is not None:
            outelements.connectivity[cpt:(2*cpt),:] = (outelements.connectivity[0:cpt,:]+pcpt)[:,permutation]
            pcpt = pcpt *2
            cpt = cpt*2

        outelements.cpt = cpt


    return outmesh

def ExtractElementByTags(inmesh,tagsToKeep, allNodes=False):

    outmesh = type(inmesh)()

    outmesh.nodes = np.copy(inmesh.nodes)
    outmesh.originalIDNodes = np.copy(inmesh.originalIDNodes)

    import copy
    for tag in inmesh.nodesTags:
        outmesh.nodesTags.AddTag(copy.deepcopy(tag) )


    nodalMask = np.zeros(inmesh.GetNumberOfNodes(),dtype = np.bool)
    for name,elems in inmesh.elements.iteritems():

       if np.any([x in elems.tags.keys() for x in tagsToKeep] ) == False:
           if np.any([x in inmesh.nodesTags.keys() for x in tagsToKeep]) == False:
               continue# pragma: no cover


       toKeep = np.zeros(elems.GetNumberOfElements(), dtype=np.bool)

       # check elements tags
       for tagToKeep in tagsToKeep:
           if elems.tags.has_key(tagToKeep):
               elems.tags[tagToKeep].tighten()
               toKeep[elems.tags[tagToKeep].id] = True

       # check for nodes tags
       for tagToKeep in tagsToKeep:
           if inmesh.nodesTags.has_key(tagToKeep):
             nodalMask.fill(False)
             tag = inmesh.GetNodalTag(tagToKeep)
             tag.tighten()
             nodalMask[tag.id] = True
             elemMask = np.sum(nodalMask[elems.connectivity],axis=1)
             if allNodes :
                 toKeep[elemMask == elems.GetNumberOfNodesPerElement()] = True
             else:
                 toKeep[elemMask > 0] = True


       newIndex = np.empty(elems.GetNumberOfElements(), dtype=np.int )
       cpt =0;
       for i in xrange(elems.GetNumberOfElements()):
           newIndex[i] = cpt
           cpt += 1 if toKeep[i] else 0



       outelem = outmesh.GetElementsOfType(name)
       nbTokeep = np.sum(toKeep)
       outelem.Allocate(nbTokeep)
       outelem.connectivity = elems.connectivity[toKeep,:]

       for tag in elems.tags  :
           tag.tighten()
           #print(tag.id)
           temp = np.extract(toKeep[tag.id],tag.id)
           #print(temp)
           newid = newIndex[temp]
           #print(newid)
           outelem.tags.CreateTag(tag.name).SetIds(newid)


    CleanLonelyNodes(outmesh)
    return outmesh

def GetVolume(inmesh) :

    def VolumeOfTetrahedrons(nodes,tets):
        lvol = 0.
        conn = elems.connectivity
        a = inmesh.nodes[conn[:,0],:]
        b = inmesh.nodes[conn[:,1],:]
        c = inmesh.nodes[conn[:,2],:]
        d = inmesh.nodes[conn[:,3],:]
        e = np.cross(b-d,c-d)
        f = (a-d)
        for n in xrange(elems.GetNumberOfElements()):
            lvol += np.abs( np.dot(f[n,:],e[n,:])  )
        #vol += (1./6.)* np.sum(np.abs( np.dot( a-d,np.cross(b-d,c-d) )  ) )
        return lvol*(1./6.)


    vol = 0;
    for name,elems in inmesh.elements.iteritems():
        if ElementNames.dimension[name] != 3:
            continue# pragma: no cover
        if name == ElementNames.Tetrahedron_4:
            vol += VolumeOfTetrahedrons(inmesh.nodes,elems.connectivity)
        else:
            raise Exception('code me please...')# pragma: no cover
    return vol


def CleanEmptyTags(inmesh):
    for name,elems in inmesh.elements.iteritems():
        elems.tags.RemoveEmptyTags()
    inmesh.nodesTags.RemoveEmptyTags()

# to generate one tag per body
# body is defines by all the nodes connected by the elements
def AddTagPerBody(inmesh):

    # generation of the dual graph
    maxNumConnections = 100

    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=int )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=int );

    for name,elems in inmesh.elements.iteritems():
        for i in xrange(elems.GetNumberOfElements()):
            coon = elems.connectivity[i,:]
            size = elems.GetNumberOfNodesPerElement()
            for j in xrange(size):
                myIndex = coon[j]
                for k in xrange(size):
                    if k == j:
                        continue
                    dualGraph[myIndex,usedPoints[myIndex]] =  coon[k]
                    usedPoints[myIndex] += 1
                    if usedPoints[myIndex] == maxNumConnections :
                        c = np.unique(dualGraph[myIndex,:])
                        dualGraph[myIndex,0:len(c)] = c
                        usedPoints[myIndex] = len(c)


    # Connectivity walk
    nbOfNodes = inmesh.GetNumberOfNodes()
    treated = np.zeros(inmesh.GetNumberOfNodes(),dtype=np.bool)
    nextpoint = np.zeros(inmesh.GetNumberOfNodes(),dtype=np.int_)
    nextpointcpt = 0

    #print(dualGraph)
    # loop over the bodies
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
            #in the case we have only False in the treated
            if treated[initialPoint]:
                break
            treated[initialPoint] = True
            nextpoint[nextpointcpt] = initialPoint
            nextpointcpt +=1



            tagName = "Body_"+str(bodyCpt)
            #print("New Tag")
            tag = inmesh.GetNodalTag(tagName)
            bodyCpt += 1
            pointsInThisBody = 1


        else:
            raise Exception ("Error in this function")# pragma: no cover


        while cpt < nextpointcpt:

            # cover all the nodes connected to this node
            workingIndex = nextpoint[cpt]
            indexes = dualGraph[workingIndex,0:usedPoints[workingIndex]]
            #print(indexes)
            for index in indexes:
                #print(index)

                if not treated[index] :
                    #sol[index] = fillField[index]
                    treated[index] = True
                    #if phi[index] <= 0.:
                    nextpoint[nextpointcpt] = index
                    nextpointcpt +=1
                    tag.AddToTag(index)
                    #print(tag.id)
                    #raise
                    pointsInThisBody +=1

            cpt += 1
        pointsPerBody.append(pointsInThisBody)
    return pointsPerBody


def MeshToVtk(mesh, vtkobject=None):

    vtknumbers = {}
    vtknumbers[ElementNames.Bar_2] = 3

    vtknumbers[ElementNames.Triangle_3] = 5
    vtknumbers[ElementNames.Quadrangle_4] = 9
    vtknumbers[ElementNames.Hexaedron_8] = 12
    vtknumbers[ElementNames.Bar_3] = 21
    vtknumbers[ElementNames.Triangle_6] = 22
    vtknumbers[ElementNames.Tetrahedron_4] = 10
    vtknumbers[ElementNames.Tetrahedron_10] = 24

    try:
        from paraview import vtk
    except :
        import vtk

    if vtkobject is not None:
        output = vtkobject
    else:
        output = vtk.vtkUnstructuredGrid()


    output.Allocate(mesh.GetNumberOfElements())
    ##copy points
    pts = vtk.vtkPoints()
    pts.Allocate(mesh.GetNumberOfNodes())
    if mesh.nodes.shape[1] == 3 :
        for p in xrange(mesh.GetNumberOfNodes()):
            point = mesh.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],point[2])
    else:
        #2DCase
        for p in xrange(mesh.GetNumberOfNodes()):
            point = mesh.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],0.0)

    output.SetPoints(pts)

    if hasattr(mesh,"nodeFields"):
        for name,data in mesh.nodeFields.iteritems():
            pd = vtk.vtkFloatArray()
            pd.SetName(name)
            pd.SetNumberOfComponents(data.shape[1])
            pd.SetNumberOfTuples(mesh.GetNumberOfNodes())
            cpt = 0
            for i in xrange(mesh.GetNumberOfNodes()):
                for j in xrange(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1

            output.GetPointData().AddArray(pd)

    for elementsname,elementContainer in mesh.elements.iteritems():
        pointIds = vtk.vtkIdList()
        npe = elementContainer.GetNumberOfNodesPerElement()
        pointIds.SetNumberOfIds(npe)
        vtknumber = vtknumbers[elementsname]
        for e in xrange(elementContainer.GetNumberOfElements()):
            for i in xrange(npe):
                pointIds.SetId(i,elementContainer.connectivity[e,i])
            output.InsertNextCell(vtknumber, pointIds)

    if hasattr(mesh,"elemFields"):
        for name,data in mesh.elemFields.iteritems():
            pd = vtk.vtkFloatArray()
            pd.SetName(name)
            pd.SetNumberOfComponents(data.shape[1])
            pd.SetNumberOfTuples(mesh.GetNumberOfElements())
            cpt = 0
            for i in xrange(mesh.GetNumberOfElements()):
                for j in xrange(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1

            output.GetCellData().AddArray(pd)
    return output

###############################################################################
def CheckIntegrity_MeshToVtk():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    sol = MeshToVtk(res)
    return "OK"


def CheckIntegrity_CreateMeshOfTriangles():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    return "OK"

def CheckIntegrity_CreateMeshFromConstantRectilinearMesh():
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

    return "OK"

def CheckIntegrity_QuadToLin():
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

def CheckIntegrity_CleanDoubleNodes():
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    CleanDoubleNodes(mesh)
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
    return "ok"

def CheckIntegrity_MirrorMesh():
    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
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


def CheckIntegrity_CleanEmptyTags():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetNodalTag("Point0").AddToTag(0)
    res.GetNodalTag("EmptyTagN")

    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("EmptyTagE")
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("Element0").AddToTag(0)

    CleanEmptyTags(res)
    if len(res.nodesTags) != 1 :
        raise# pragma: no cover

    if len(res.GetNamesOfCellTags()) != 1 :
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_AddTagPerBody():
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    resII = MirrorMesh(res,x=0,y=0,z=0)
    print(AddTagPerBody(resII))

    if len(resII.nodesTags) != 8:
        raise# pragma: no cover
    print(resII.nodesTags)

    return "ok"


def CheckIntegrity():
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
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover