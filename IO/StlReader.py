# -*- coding: utf-8 -*-
import numpy as np

__author__ = "Felipe Bordeu"
import BasicTools.FE.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase


def ReadStl(fileName=None,string=None):
    obj = StlReader()
    obj.SetFileName(fileName)
    obj.SetStringToRead(string)
    res = obj.Read()
    return res

def LoadSTLWithVTK(filenameSTL):
    import vtk
    readerSTL = vtk.vtkSTLReader()
    readerSTL.SetFileName(filenameSTL)
    # 'update' the reader i.e. read the .stl file
    readerSTL.Update()

    polydata = readerSTL.GetOutput()

    # If there are no points in 'vtkPolyData' something went wrong
    if polydata.GetNumberOfPoints() == 0:
        raise ValueError("No point data could be loaded from '" + filenameSTL)# pragma: no cover

    return polydata

class StlReader(ReaderBase):
    def __init__(self,fileName = None):
        super(StlReader,self).__init__()

    def Read(self, fileName=None,string=None,out=None):

      if fileName is not None:
          self.SetFileName(fileName)

      if string is not None:
          self.SetStringToRead(string)


      # check if binary
      self.readFormat = "rb"
      self.StartReading()

      header = ""
      #read the first non space caracters
      while len(header) < 5:
          header += "".join([ x.decode("utf-8",'ignore') for x in self.filePointer.read(1).split() ])
      self.EndReading()

      if header == "solid":
          self.PrintVerbose("Ascii File")
          return self.ReadStlAscii()
      else:
          self.PrintVerbose("Binary File")
          return self.ReadStlBinary()


    def ReadStlBinary(self):
        # from https://en.wikipedia.org/wiki/STL_(file_format)#Binary_STL

        self.readFormat = "rb"
        self.StartReading()

        import BasicTools.FE.UnstructuredMesh as UM

        header = self.readData(80,np.int8)
        header = ''.join([chr(item) for item in header])

        print("HEADER : '" + header + "'")
        nbTriangles = self.readInt32()
        print("reading  : " + str(nbTriangles) + " triangles")

        resUM = UM.UnstructuredMesh()

        #resUM.nodes = np.empty((nbTriangles*3,3), dtype=float)

        dt = np.dtype([('normal', (np.float32,3)),
                       ('points', (np.float32,9)),
                       ('att', (np.int16)),
                       ])

        data = np.fromfile(self.filePointer,dtype=dt,count=nbTriangles,sep="")
        normals = np.array(data["normal"])
        resUM.nodes = np.array(data["points"])
        resUM.nodes.shape = (nbTriangles*3,3)

#        cpt = 0
#        for i in range(nbTriangles):
#            data = self.readFloats32(12)
#            normals[i,:] = data[0:3]
#            resUM.nodes[cpt,:] = data[3:6]
#            cpt +=1
#            resUM.nodes[cpt,:] = data[6:9]
#            cpt +=1
#            resUM.nodes[cpt,:] = data[9:12]
#            cpt +=1
#
#            #normals[i,:] = self.readFloats32(3)
#            #for j in range(3):
#            #    resUM.nodes[cpt,:] = self.readFloats32(3)
#            #    cpt += 1
#            # trash the  Attribute byte count
#            self.readData(1,np.int16)


        elements = resUM.GetElementsOfType(EN.Triangle_3)
        elements.connectivity = np.array(range(resUM.GetNumberOfNodes()),dtype=np.int)
        elements.connectivity.shape = (nbTriangles,3)
        elements.originalIds = np.arange(nbTriangles,dtype=np.int )
        elements.cpt = nbTriangles
        resUM.elemFields["normals"] = normals
        self.EndReading()

        return resUM

    def ReadStlAscii(self):

        self.readFormat = "r"
        self.StartReading()


        import BasicTools.FE.UnstructuredMesh as UM

        resUM = UM.UnstructuredMesh()

        name = self.ReadCleanLine().split()[1];

        p = []
        normals = np.empty((0,3), dtype=float)
        nodesbuffer = []
        while True:
            line = self.ReadCleanLine()
            if not line:
                break

            l = line.strip('\n').lstrip().rstrip()
            if l.find("facet")>-1 :
                if l.find("normal")>-1 :
                 normals = np.concatenate((normals, np.fromstring(l.split("normal")[1],sep=" ")[np.newaxis]),axis=0)
                 continue
            if l.find("outer loop")>-1 :
              for i in range(3):
                line = self.ReadCleanLine()
                l = line.strip('\n').lstrip().rstrip()
                if l.find("vertex")>-1 :
                  p.append(np.fromstring(l.split("vertex")[1],sep=" ") )
              if len(p) == 3:
                nodesbuffer.extend(p)
                #resUM.nodes = np.vstack((resUM.nodes,p[0][np.newaxis,:],p[1][np.newaxis,:],p[2][np.newaxis,:]))
                p = []
              else:# pragma: no cover
                print("error: outer loop with less than 3 vertex")
                raise

        self.EndReading()


        resUM.nodes = np.array(nodesbuffer)
        del nodesbuffer
        elements = resUM.GetElementsOfType(EN.Triangle_3)
        elements.connectivity = np.array(range(resUM.GetNumberOfNodes()),dtype=np.int)
        elements.connectivity.shape = (resUM.GetNumberOfNodes()//3,3)
        #elements.connectivity = elements.connectivity.T
        elements.originalIds = np.arange(resUM.GetNumberOfNodes()/3,dtype=np.int )
        elements.cpt = elements.connectivity.shape[0]
        resUM.elemFields["normals"] = normals
        return resUM






def CheckIntegrity():
    data = u"""   solid cube_corner
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
        endsolid
"""


    res = ReadStl(string=data)
    print(res)
    if res.GetNumberOfNodes() != 12: raise Exception()
    if res.GetNumberOfElements() != 4: raise Exception()

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    f =open(tempdir+"test_input_stl_data.stl","w")
    f.write(data)
    f.close()
    res = ReadStl(fileName=tempdir+"test_input_stl_data.stl")

    try:
        import vtk
        print('reading mesh using vtk')
        mesh = LoadSTLWithVTK(tempdir+"test_input_stl_data.stl")
        print(mesh)
    except:# pragma: no cover
        print("warning : error importing vtk, disabling some tests ")
        pass
    from BasicTools.TestData import GetTestDataPath


    print("Binary reading")
    res1 = ReadStl(fileName=GetTestDataPath()+"coneBinary.stl")
    print (res1)

    print("Ascii reading")
    res2 = ReadStl(fileName=GetTestDataPath()+"coneAscii.stl")
    print (res2)

    if res1.GetNumberOfNodes() != res2.GetNumberOfNodes(): raise Exception()

    # 1e-6 because the ascii file has only 6 decimals
    if np.max(abs(res1.nodes -  res2.nodes)) > 1e-6 : raise Exception()

    if res1.GetNumberOfElements() !=res2.GetNumberOfElements(): raise Exception()

    conn1 = res1.GetElementsOfType(EN.Triangle_3).connectivity
    conn2 = res2.GetElementsOfType(EN.Triangle_3).connectivity

    if not np.all(np.equal(conn1,conn2)) : raise Exception()

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover






