# -*- coding: utf-8 -*-

import numpy as np

class STL(object):
    def __init__(self):
       self.tris = []
       self.name = ""
       self.boundingMin = np.array([1e10, 1e10, 1e10],dtype= float)
       self.boundingMax = np.array([-1e10, -1e10, -1e10],dtype= float)
     

    class Tri():
        def __init__(self,p0=None,p1=None,p2=None,normal=None):
            if p0 is None:
                self.p0 = np.zeros((3,),dtype= float)
            else:
                self.p0 = np.array(p0,dtype= float)
                
            if p1 is None:
                self.p2 = np.zeros((3,),dtype= float)
            else:
                self.p1 = np.array(p1,dtype= float)
                
            if p2 is None:
                self.p2 = np.zeros((3,),dtype= float)
            else:
                self.p2 = np.array(p2,dtype= float)
                
            if normal is None:
                self.normal = np.zeros((3,),dtype= float)
            else:
                self.normal = np.array(normal,dtype= float)

        def Renormalise(self):
            self.normal = np.cross(self.p1-self.p0,self.p2-self.p0)
            self.normal = self.normal/np.linalg.norm(self.normal)
                
    def AddTri(self,p0,p1,p2,normal):
        self.tris.append(STL.Tri(p0,p1,p2,normal));
        self.boundingMin = np.minimum(np.minimum(np.minimum(p0,p1),p2),self.boundingMin);
        self.boundingMax = np.maximum(np.maximum(np.maximum(p0,p1),p2),self.boundingMax);
        
    def ComputeNormals(self):
        for t in self.tris:
            t.Renormalise();
        
def WriteStl(data,output= None):    
    import sys
    if output is None:
        output = sys.stdout
    
    output.write("solid {}\n".format(data.name))          
    for t in data.tris:
        output.write(" facet normal {}\n".format(" ".join(map(str,t.normal)) ))
        output.write("  outer loop\n")
        output.write("   vertex {}\n".format(" ".join(map(str,t.p0))))
        output.write("   vertex {}\n".format(" ".join(map(str,t.p1))))
        output.write("   vertex {}\n".format(" ".join(map(str,t.p2))))
        output.write("  endloop\n")
        output.write(" endfacet\n")
    output.write("endsolid\n")          
                 
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

    res = STL()
    l = string.readline();
    res.name = l.strip('\n').lstrip().rstrip().split()[1]
        
    p = []
    for line in string:
        l = line.strip('\n').lstrip().rstrip()
        if l.find("facet")>-1 : 
            if l.find("normal")>-1 : 
             normal = np.fromstring(l.split("normal")[1],sep=" ")
             continue
        if l.find("outer loop")>-1 :  
          for i in range(3):
            line = string.readline()
            l = line.strip('\n').lstrip().rstrip()
            if l.find("vertex")>-1 : 
              p.append(np.fromstring(l.split("vertex")[1],sep=" ") )
          if len(p) == 3:
            print(resUM.nodes.shape)
            print(p[0][np.newaxis,:].shape)
            print(p[1][np.newaxis,:].shape)
            print(p[2][np.newaxis,:].shape)
            
            resUM.nodes = np.vstack((resUM.nodes,p[0][np.newaxis,:],p[1][np.newaxis,:],p[2][np.newaxis,:]))
            res.AddTri(p0 =p[0],p1=p[1],p2=p[2],normal=normal ) 
            p = []
          else:
              print("error: outer loop with less than 3 vertex")

    elements = resUM.GetElementsOfType('tri')    
    elements.connectivity = np.array(xrange(resUM.GetNumberOfNodes()),dtype=np.int)
    elements.connectivity.shape = (3,resUM.GetNumberOfNodes()/3)
    elements.connectivity = elements.connectivity.T 
    elements.originalIds = np.arange(resUM.GetNumberOfNodes()/3,dtype=np.int )
    print(resUM)

    return res

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

    return 'ok'
        
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
    



def vtk_show(renderer, width=400, height=300):
    import vtk
    from IPython.display import Image
    """
    Takes vtkRenderer instance and returns an IPython Image with the rendering.
    """
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(width, height)
    renderWindow.Render()
     
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()
     
    writer = vtk.vtkPNGWriter()
    writer.SetWriteToMemory(1)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
    data = str(buffer(writer.GetResult()))
    
    return Image(data)
    
if __name__ == '__main__':
    #print(CheckIntegrity())# pragma: no cover   

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
    print("------------------------")
    WriteStl(res)
    print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---------")
    res.ComputeNormals()
    WriteStl(res)
    
    
    import vtk
    import OTTools.TestData as T2
    tempPath = T2.GetTestDataPath()
    mesh = LoadSTLWithVTK(tempPath+"../../stlexample.stl")


    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(0.25)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1.0, 1.0, 1.0)
    vtk_show(renderer)
    
    from pycaster import pycaster
    caster = pycaster.rayCaster.fromSTL(tempPath+"../../stlexample.stl", scale=1)