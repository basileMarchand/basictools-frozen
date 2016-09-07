# -*- coding: utf-8 -*-
import OTTools.FE.ElementNames as EN
import numpy as np

def ReadStl(fileName=None,string=None):
    from cStringIO import StringIO
    import OTTools.FE.UnstructuredMesh as UM

    #from collections import OrderedDict as OD

    if fileName is not None:
        f = open(fileName, 'r')
        string = f.read()
        f.close()
    
    string = StringIO(string)
    
    resUM = UM.UnstructuredMesh()

    l = string.readline();
    name = l.strip('\n').lstrip().rstrip().split()[1]
        
    p = []
    for line in string:
        l = line.strip('\n').lstrip().rstrip()
        if l.find("facet")>-1 : 
            if l.find("normal")>-1 : 
             #normal = np.fromstring(l.split("normal")[1],sep=" ")
             continue
        if l.find("outer loop")>-1 :  
          for i in range(3):
            line = string.readline()
            l = line.strip('\n').lstrip().rstrip()
            if l.find("vertex")>-1 : 
              p.append(np.fromstring(l.split("vertex")[1],sep=" ") )
          if len(p) == 3:
            resUM.nodes = np.vstack((resUM.nodes,p[0][np.newaxis,:],p[1][np.newaxis,:],p[2][np.newaxis,:]))
            p = []
          else:
            print("error: outer loop with less than 3 vertex")
            raise

    elements = resUM.GetElementsOfType(EN.Triangle_3)    
    elements.connectivity = np.array(xrange(resUM.GetNumberOfNodes()),dtype=np.int)
    elements.connectivity.shape = (resUM.GetNumberOfNodes()/3,3)
    #elements.connectivity = elements.connectivity.T 
    elements.originalIds = np.arange(resUM.GetNumberOfNodes()/3,dtype=np.int )

    return resUM


        
def LoadSTLWithVTK(filenameSTL):
    import vtk
    readerSTL = vtk.vtkSTLReader()
    readerSTL.SetFileName(filenameSTL)
    # 'update' the reader i.e. read the .stl file
    readerSTL.Update()

    polydata = readerSTL.GetOutput()

    # If there are no points in 'vtkPolyData' something went wrong
    if polydata.GetNumberOfPoints() == 0:
        raise ValueError(
            "No point data could be loaded from '" + filenameSTL)
        return None
    
    return polydata
    

def CheckIntegrity():
    data = """   solid cube_corner
          facet normal 0.0 -1.0 0.0
            outer loop
              vertex 0.0 0.0 0.0
              vertex 1.0 0.0 0.0
              vertex 0.0 0.0 1.0
            endloop
          endfacet
          facet normal 0.0 0.0 -1.0
            outer loop
              vertex 0.0 0.0 0.0
              vertex 0.0 1.0 0.0
              vertex 1.0 0.0 0.0
            endloop
          endfacet
          facet normal -1.0 0.0 0.0
            outer loop
              vertex 0.0 0.0 0.0
              vertex 0.0 0.0 1.0
              vertex 0.0 1.0 0.0
            endloop
          endfacet
          facet normal 0.577 0.577 0.577
            outer loop
              vertex 1.0 0.0 0.0
              vertex 0.0 1.0 0.0
              vertex 0.0 0.0 1.0
            endloop
          endfacet
        endsolid"""
        
    
    res = ReadStl(string=data)
    print(res)
    if res.GetNumberOfNodes() != 12: raise Exception()
    if res.GetNumberOfElements() != 4: raise Exception()
        
    
    try:
        import vtk
        import OTTools.TestData as T2
        print('reading mesh using vtk')
        mesh = LoadSTLWithVTK(T2.GetTestDataPath()+"stlexample.stl")
        print(mesh)
    except:
        pass
        
    return 'ok'
  
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
    

    
    
