# -*- coding: utf-8 -*-

import OTTools.FE.ElementNames as ElementNames


vtknumbers = {} 
vtknumbers[ElementNames.Bar_2] = 3

vtknumbers[ElementNames.Triangle_3] = 5
vtknumbers[ElementNames.Quadrangle_4] = 9
vtknumbers[ElementNames.Hexaedron_8] = 12
vtknumbers[ElementNames.Bar_3] = 21
vtknumbers[ElementNames.Triangle_6] = 22
vtknumbers[ElementNames.Tetrahedron_10] = 24


def ReadMesh(filename):
    extention = filename.split(".")[-1].lower()
    if extention ==  "asc":
        import OTTools.IO.AscReader as AscReader
        return AscReader.ReadAsc(filename)
    elif extention ==  "geof":
        import OTTools.IO.GeofReader as GeofReader
        return GeofReader.ReadGeof(filename)
    elif extention ==  "msh":
        import OTTools.IO.GmshReader as GmshReader
        return GmshReader.ReadGmsh(filename)
    elif extention ==  "inp":
        import OTTools.IO.InpReader as ImpReader
        return ImpReader.ReadInp(filename)
    else:
        raise Exception ("Unkown file extention : " + str(extention))


## to use this function add this lines to the 
## programmmable source in paraview
##
#
#from OTTools.IO.UniversalReader import ReadMeshAnPopulateVtkObject as ReadMeshAnPopulateVtkObject
#filename = "here you put your filename"
#
#ReadMeshAnPopulateVtkObject(filename,self.GetOutput())

def ReadMeshAndPopulateVtkObject(filename, vtkobject= None):
    res = ReadMesh(filename)
    try:
        from paraview import vtk
    except :
        import vtk
        
    if vtkobject is not None:
        output = vtkobject
    else:
        output = vtk.vtkUnstructuredGrid()
    
    
    output.Allocate(res.GetNumberOfElements())
    ##copy points        
    pts = vtk.vtkPoints()
    pts.Allocate(res.GetNumberOfNodes())
    if res.nodes.shape[1] == 3 :
        for p in xrange(res.GetNumberOfNodes()):
            point = res.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],point[2])
    else:
        #2DCase
        for p in xrange(res.GetNumberOfNodes()):
            point = res.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],0.0)
    
    output.SetPoints(pts)
    
    
    for elementsname,elementContainer in res.elements.iteritems():
        pointIds = vtk.vtkIdList()
        npe = elementContainer.GetNumberOfNodesPerElement()
        pointIds.SetNumberOfIds(npe)
        vtknumber = vtknumbers[elementsname]
        for e in xrange(elementContainer.GetNumberOfElements()):
            for i in xrange(npe):
                pointIds.SetId(i,elementContainer.connectivity[e,i])
            output.InsertNextCell(vtknumber, pointIds)



