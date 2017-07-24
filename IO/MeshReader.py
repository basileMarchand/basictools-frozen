# -*- coding: utf-8 -*-
import struct
import numpy as np

__author__ = "Felipe Bordeu"
import BasicTools.FE.UnstructuredMesh  as UM
import BasicTools.FE.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase

# for ascii files
meditName = {}
meditName['Edges'] = EN.Bar_2
meditName['Triangles'] = EN.Triangle_3
meditName['Quadrilaterals'] = EN.Quadrangle_4
meditName['Tetrahedra'] = EN.Tetrahedron_4
meditName['Hexahedra'] = EN.Hexaedron_8

# for binary files
meditNumber = {}
meditNumber[5] = EN.Bar_2
meditNumber[6] = EN.Triangle_3
meditNumber[7] = EN.Quadrangle_4
meditNumber[8] = EN.Tetrahedron_4

def ReadMesh(fileName=None,string=None,ReadRefsAsField=False ):
    reader = MeshReader()
    reader.ReadRefsAsField(ReadRefsAsField)
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read()
    return reader.output


class MeshReader(ReaderBase):
    def __init__(self):
        super(MeshReader,self).__init__()
        self.refsAsAField = False
        self.dim = 3


    def ReadRefsAsField(self, val):
        self.refsAsAField = bool(val)

    def Read(self,out=None):
        if self.fileName is not None and self.fileName[-1] == "b":
            return self.ReadMeshBinary(out=out)
        else:
            return self.ReadMeshAscii(out=out)

    def ReadExtraField(self,fileName):
        if  fileName[-1] == "b":
            return self.ReadExtraFieldBinary(fileName)
        else:
            return self.ReadExtraFieldAscii(fileName)

    def ReadExtraFieldAscii(self,fileName):

        self.PrintVerbose("Reading ExtraField")
        f =  open(fileName, "r")

        res = {}


        dataType = float
        while(True):
            line = f.readline()
            if line == "" :
                break

            l = line.strip('\n').lstrip().rstrip()
            if len(l) == 0: continue

            if l.find("MeshVersionFormatted")>-1 :
                if len(l.split()) > 1:
                    formatFile = int(l.split()[1])
                else:
                    formatFile = int(f.readline())

                if formatFile == 2:
                    dataType = np.float64
                #print ('formatFile  : '+ str(dataType))
                continue

            if l.find("Dimension")>-1 :
                s = l.split()
                if len(s)>1:
                    s.pop(0)
                    l = " ".join(s)
                else:
                    l = f.readline()
                self.dim = int(l)
                #print ('Dimension : '+ str(dimension))
                continue

            def ReadFieldsASCII(f):
                datares = []
                line = f.readline()
                #print(line)
                l = line.strip('\n').lstrip().rstrip()
                #print(l)
                nbentries = int(l.split()[0])
                line = f.readline()
                l = line.strip('\n').lstrip().rstrip()
                nbfields = int(l.split()[0])
                #print ("l"),
                #print (l)
                fieldSizes = [ int(x) for x in l.split()[1:] ]
                for i in range(len(fieldSizes)):
                    if fieldSizes[i] == 2 :
                       fieldSizes[i] = self.dim
                    elif fieldSizes[i] == 3 :
                        raise


                #print ("fieldSizes"),
                #print (fieldSizes)
                ncoomp = np.sum(fieldSizes)
                self.PrintDebug("nbfields" + str(fieldSizes))
                #print("nbentries "+str(nbentries)+ " ")
                #print("ncoomp "+str(ncoomp)+ " ")
                #print("Reading "+str(ncoomp*nbentries)+ " Entries")
                data = np.fromfile(f,dtype=float,count=int(ncoomp*nbentries),sep=" ")
                data.shape = (nbentries,ncoomp)
                cpt = 0
                for i in range(nbfields):
                    datares.append(data[:,cpt:cpt+fieldSizes[i]])
                    cpt += fieldSizes[i]
                return datares

            for fieldName in  ["SolAtVertices", "SolAtTetrahedra" ]:
                found = False
                if l.find(fieldName)>-1 and len(l) == len(fieldName):
                    #print("reading SolAtVertices")
                    data = ReadFieldsASCII(f)
                    for i in range(len(data)) :
                        res[fieldName+str(i)]  =  data[i]
                    found = True
                    break

            if found:
                continue

            if l.find("End")>-1 :
                break


            print("Ignoring entrie " + str(l) )


        return res

    def _readExtraFieldBinary(self,f,dataSize):
      data = f.read(4)
      endOfInformation = struct.unpack("i", data)[0]

      data = f.read(4)
      nbNodes = struct.unpack("i", data)[0]

      data = f.read(4)
      nbfields = struct.unpack("i", data)[0]

      fieldSizes = np.empty(nbfields, dtype=np.int)
      #self.PrintDebug("nbfields" + str(nbfields))
      res = []
      for i in range(nbfields):
          data = f.read(4)
          fieldSizes[i] = struct.unpack("i", data)[0]

      for i in range(nbfields):
          dataformat = ('f' if dataSize==4 else "d"    )

          if fieldSizes[i] == 2 :
             fieldSizes[i] = self.dim
          elif fieldSizes[i] == 3 :
             raise

          data = f.read(dataSize*nbNodes*fieldSizes[i])

          field = np.empty((nbNodes,fieldSizes[i]),dtype=np.float)

          data = struct.unpack(dataformat*(nbNodes*fieldSizes[i]), data)
          #print(data)
          cpt =0
          for n in range(nbNodes):
              for j in range(fieldSizes[i]):
                  field[n,j] = data[cpt]
                  cpt += 1

          if fieldSizes[i] == 1 :
             field.shape = (len(field))
          res.append(field)
      return res



    def ReadExtraFieldBinary(self,fileName):
      self.PrintVerbose("Reading ExtraField")
      f =  open(fileName, "rb")


      dataType = np.float32
      dataSize = 4;

      data = f.read(4)
      # end of file
      if len(data) == 0:
          raise Exception("Empty file")
      key = struct.unpack("i", data)[0]

     #MeshVersionFormatted
      if key == 1:
          data = f.read(4)
          udata = struct.unpack("i", data)[0]
          if udata == 2 :
              dataType = np.float64
              dataSize = 8
              #print("data in double")
      else:
          raise Exception('Expected key value equal to 1 (binary), are you sure this file is binary??')

      #if self.refsAsAField:
      #   elemRefs = np.zeros(0 ,dtype=np.int)
      res = {}
      while True:
          data = f.read(4)

          # end of file
          if len(data) == 0:
             #print("EOF")
             break
          key = struct.unpack("i", data)[0]
          #print("key : " +str(key))
          #Dimension
          if key == 3:
              data = f.read(4)
              endOfInformation = struct.unpack("i", data)[0]
              #f.seek(endOfInformation-4)
              data = f.read(4)
              dimension = struct.unpack("i", data)[0]
              f.seek(endOfInformation)
              #print("dimension : " + str(dimension))
              continue

          #SolAtVertices
          if key == 62:
              data = self._readExtraFieldBinary(f,dataSize)
              for i in range(len(data)) :
                  res["SolAtVertices"+str(i)]  =  data[i]
              continue
          #SolAtTriangles
          if key == 64:
              data = self._readExtraFieldBinary(f,dataSize)
              for i in range(len(data)) :
                  res["SolAtTriangles"+str(i)]  =  data[i]
              continue
          #SolAtTetrahedra
          if key == 66:
              data = self._readExtraFieldBinary(f,dataSize)
              for i in range(len(data)) :
                  res["SolAtTetrahedra"+str(i)]  =  data[i]
              continue

          #End key
          if key == 54:
             break

          self.PrintVerbose("skiping key : "  + str(key) )
          data = f.read(4)
          endOfInformation = struct.unpack("i", data)[0]
          f.seek(endOfInformation)
      return res

    def ReadMeshBinary(self,out=None):
      self.readFormat = 'rb'
      self.StartReading()


      f = self.filePointer
      if out is None:
          res = UM.UnstructuredMesh()
      else:
          res = out

      self.output = res
      dataType = np.float32
      dataSize = 4;

      data = f.read(4)
      # end of file
      if len(data) == 0:
          raise Exception("Empty file")
      key = struct.unpack("i", data)[0]

     #MeshVersionFormatted
      if key == 1:
          data = f.read(4)
          udata = struct.unpack("i", data)[0]
          if udata == 2 :
              dataType = np.float64
              dataSize = 8
              #print("data in double")
      else:
          raise Exception('Expected key value equal to 1 (binary), are you sure this file is binary??')

      globalElementCounter = 0;
      #if self.refsAsAField:
      #   elemRefs = np.zeros(0 ,dtype=np.int)
      elemRefsDic = {}

      while True:
          data = f.read(4)
          # end of file
          if len(data) == 0:
             break
          key = struct.unpack("i", data)[0]

          #Dimension
          if key == 3:
              data = f.read(4)
              endOfInformation = struct.unpack("i", data)[0]
              #f.seek(endOfInformation-4)
              data = f.read(4)
              dimension = struct.unpack("i", data)[0]
              f.seek(endOfInformation)
              #print("dimension : " + str(dimension))
              continue

          #Vertices
          if key == 4:
              data = f.read(4)
              endOfInformation = struct.unpack("i", data)[0]
              #print(endOfInformation)
              # I dont undestant why I have to trash 4 bytes
              data = f.read(4)
              nbNodes = struct.unpack("i", data)[0]

              res.nodes = np.empty((nbNodes,dimension),dtype=dataType)
              res.originalIDNodes= np.empty((nbNodes,),dtype=np.int)

              self.PrintVerbose("Reading " + str(nbNodes) + " nodes ")
              cpt =0;
              dataformat = ('f' if dataSize==4 else "d"    )*dimension + "i"
              if self.refsAsAField:
                  nodalrefs = np.zeros(nbNodes ,dtype=np.int)
                  self.nodalFields['refs'] = nodalrefs


              while(True):
                  data = f.read(dataSize*dimension+4)
                  #print(('f' if dataSize==4 else "d"    )*dimension + "i")
                  udata = struct.unpack(dataformat,data)
                  res.nodes[cpt,:] = udata[0:dimension]
                  res.originalIDNodes[cpt] = int(cpt)

                  #print(udata)
                  if self.refsAsAField:
                      ref = udata[dimension]
                      nodalrefs[cpt] = ref
                      #print (str(cpt) +' : '+ str(ref) + ' | '),
                  else:
                      ref = str(udata[dimension])
                      res.GetNodalTag("NTag"+ref).AddToTag(cpt)

                  cpt +=1
                  #self.PrintProgress(cpt,maxx = nbNodes)

                  #if (cpt/10000 == float(cpt)/10000):
                    #self.PrintProgress()
                  #  print(100*float(cpt)/nbNodes)
                  if cpt == nbNodes:
                    break
              #print(f.tell())
              del cpt
              continue



          # all kind of elements
          if key in meditNumber:
              elements = res.GetElementsOfType(meditNumber[key])
              self.PrintVerbose("Reading elements of type " + elements.elementType )
              nbNodes = elements.GetNumberOfNodesPerElement()

              data = f.read(4)
              endOfInformation = struct.unpack("i", data)[0]


              # I dont undestant why I have to trask 4 bites

              data = f.read(4)
              nbElements = struct.unpack("i", data)[0]
              elements.Reserve(nbElements)
              #(endOfInformation- beginOfInformation)/((nbNodes+1)*5)
              self.PrintVerbose('Reading '+ str(nbElements) + " Elements ")
              elements.cpt =0;
              dataformat = "i"*(1+nbNodes)

              if self.refsAsAField:
                    elemRefs = np.zeros(nbElements ,dtype=np.int)
                    elemRefsDic[elements.elementType] = elemRefs

              while(elements.cpt < nbElements):
                data = f.read((nbNodes+1)*4)
                udata = struct.unpack(dataformat,data)
                conn = udata[:nbNodes]
                elements.connectivity[elements.cpt,:] = conn;
                elements.originalIds[elements.cpt] = globalElementCounter;


                #elements.AddNewElement(conn, cpt)



                if self.refsAsAField:
                    ref = int(udata[nbNodes])
                    #if globalElementCounter >= len(elemRefs) :
                    #    elemRefs = np.resize(elemRefs, globalElementCounter*2+1)
                    elemRefs[elements.cpt] = ref
                else:
                    ref = str(udata[nbNodes])
                    elements.GetTag("ETag"+ref).AddToTag(elements.cpt)

                globalElementCounter +=1
                elements.cpt +=1
                #if nbElements == cpt:
                #    break

                #self.PrintProgress(elements.cpt, maxx = nbElements )
                #if (elements.cpt/10000 == float(elements.cpt)/10000):
                #    print(elements.cpt)


              elements.connectivity -= 1;
              self.PrintVerbose('Reading '+ str(nbElements) + " Elements Done ")


              continue

          #ridges
          if key == 14:
              endOfInformation = struct.unpack("i", f.read(4))[0]
              nbentries = struct.unpack("i", f.read(4))[0]
              #ids = np.fromfile(f,dtype=int,count=int(nbentries),sep="")
              data = f.read(nbentries*4)
              ids = np.array(struct.unpack("i"*nbentries,data),dtype=np.int)-1
              #print(ids)
              #raise
              bars = res.GetElementsOfType(EN.Bar_2)
              bars.tags.CreateTag("Ridges").SetIds(ids)
              continue


          if key == 62:
             endOfInformation = struct.unpack("i", data)[0]
             f.seek(endOfInformation)
             continue


          #End key
          if key == 54:
             break

          self.PrintVerbose("skiping key : "  + str(key) )
          data = f.read(4)
          endOfInformation = struct.unpack("i", data)[0]
          f.seek(endOfInformation)

      res.GenerateManufacturedOriginalIDs()

      if self.refsAsAField:
         elemRefs = np.resize(elemRefs,globalElementCounter )
         cpt =0
         for name,val in res.elements.items():
             elemRefs[cpt:cpt+val.GetNumberOfElements()] = elemRefsDic[name]
             cpt += val.GetNumberOfElements()
         self.elementsFields['refs'] = elemRefs

      self.EndReading()
      return res





    def ReadMeshAscii(self, _fileName=None,_string=None, fieldFileName= None,out=None):

        if _fileName is not None:
            self.SetFileName(_fileName)

        if _string is not None:
            self.SetStringToRead(_string)

        self.readFormat = 'r'
        self.StartReading()

        if out is None:
            res = UM.UnstructuredMesh()
        else:
            res = out

        self.output = res

        dataType = float
        while(True):
            line = self.filePointer.readline()
            if line == "" :
                break
        #for line in self.filePointer:
            l = line.strip('\n').lstrip().rstrip()
            if len(l) == 0: continue

            if l.find("MeshVersionFormatted")>-1 :
                if len(l.lstrip().rstrip().split()) > 1:
                    formatFile = (l.lstrip().rstrip().split()[1])
                else:
                    formatFile = int(self.filePointer.readline())
                if formatFile == 2:
                    dataType = np.float64
                #print ('formatFile  : '+ str(dataType))
                continue

            if l.find("Dimension")>-1 :
                s = l.split()
                if len(s)>1:
                    s.pop(0)
                    l = " ".join(s)
                else:
                    l = self.filePointer.readline()
                dimension = int(l)
                self.PrintVerbose('Dimension : '+ str(dimension))
                continue

            if l.find("Vertices")>-1 and len(l) == 8:
                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()

                nbNodes = int(l.split()[0])
                self.PrintVerbose("Reading "+str(nbNodes)+ " Nodes")
                res.nodes = np.empty((nbNodes,dimension),dtype = dataType)
                #ref = UM.Tag('ref')
                res.originalIDNodes= np.empty((nbNodes,),dtype=np.int)
                cpt =0;
                while(True):
                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0: continue
                    s = l.split()

                    res.nodes[cpt,:] = list(map(float,s[0:dimension]))
                    res.originalIDNodes[int(cpt)] = int(cpt)

                    ref = s[dimension]

                    res.GetNodalTag("NTag"+ref).AddToTag(cpt)

                    cpt +=1
                    if cpt == nbNodes:
                        break
                continue


            if l in meditName:
                elements = res.GetElementsOfType(meditName[l])
                nbNodes = elements.GetNumberOfNodesPerElement()
                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()

                nbElements = int(l.split()[0])
                self.PrintVerbose("Reading "+str(nbElements)+ " Elements")
                cpt =0;
                while(True):
                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0:
                        continue

                    s = l.split()

                    elements.AddNewElement(s[0:nbNodes], cpt)
                    ref = s[nbNodes]
                    elements.GetTag("ETag"+ref).AddToTag(cpt)

                    cpt +=1
                    if nbElements == cpt:# pragma: no cover
                        break

                elements.connectivity -= 1;
                continue

            if l.find("End")>-1 :
                break

            self.PrintVerbose("ignoring line : " + l )

        self.EndReading()
        return res



def CheckIntegrity():


    __teststring = u"""
MeshVersionFormatted
2

Dimension 3
Vertices
4
0 0 0 1
1 0 0 1
0 1 0 1
0 0 1 0

Tetrahedra
1
1 2 3 4 52
End
"""


    __teststringField=u"""
MeshVersionFormatted 2

Dimension 3

SolAtVertices
4
1 1

1.
2.
3.
4.

End
"""

    res = ReadMesh(string=__teststring)
    print(res)

    from BasicTools.Helpers.Tests import TestTempDir
    newFileName = TestTempDir().GetTempPath()+"mshFile.mesh"
    open(newFileName,'w').write(__teststring)
    res = ReadMesh(fileName=newFileName)

    newFileName = TestTempDir().GetTempPath()+"mshFile.sol"
    open(newFileName,'w').write(__teststringField)
    resfield = ReadMesh(fileName=newFileName)


    from BasicTools.IO.MeshWriter import WriteMesh as WriteMesh

    newFileName = TestTempDir().GetTempPath()+"mshFile.meshb"
    WriteMesh(newFileName, res,binary=True)

    res = ReadMesh(newFileName,ReadRefsAsField=True)
    print(res)


    sol = MeshReader().ReadExtraField(TestTempDir().GetTempPath()+"mshFile.sol")
    print(sol)
    from BasicTools.IO.MeshWriter import MeshWriter as MeshWriter
    mw = MeshWriter()
    mw.SetBinary(True)
    mw.SetFileName(TestTempDir().GetTempPath()+"mshFile.solb")
    mw.OpenSolutionFile(res)
    mw.WriteSolutionsFields(res,SolsAtVertices=sol['SolAtVertices0'])
    mw.CloseSolutionFile()

    sol = MeshReader().ReadExtraField(TestTempDir().GetTempPath()+"mshFile.solb")

    return 'ok'

if __name__ == '__main__':# pragma: no cover
    from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
    BaseOutputObject.SetGlobalDebugMode(True)
    print(CheckIntegrity())# pragma: no cover
