# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
import struct
import numpy as np

import OTTools.FE.UnstructuredMesh  as UM
import OTTools.FE.ElementNames as EN
from OTTools.IO.ReaderBase import ReaderBase
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

def ReadMesh(fileName=None,string=None ):
    reader = MeshReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read()
    return reader.output


class MeshReader(ReaderBase):
    def __init__(self):
        super(MeshReader,self).__init__()
        self.refsAsAField = True


    def ReadRefsAsField(self, val):
        self.refsAsAField = val

    def Read(self):
        if self.fileName is not None and self.fileName[-1] == "b":
            self.ReadMeshBinary()
        else:
            self.ReadMeshAscii()

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
              data = f.read(4)
              endOfInformation = struct.unpack("i", data)[0]

              data = f.read(4)
              nbNodes = struct.unpack("i", data)[0]

              data = f.read(4)
              nbfields = struct.unpack("i", data)[0]

              fieldSizes = np.empty(nbfields, dtype=np.int)
              self.Print("nbfields")
              self.Print(str(nbfields))
              for i in xrange(nbfields):
                  data = f.read(4)
                  fieldSizes[i] = struct.unpack("i", data)[0]

              for i in xrange(nbfields):
                  dataformat = ('f' if dataSize==4 else "d"    )

                  data = f.read(dataSize*nbNodes*fieldSizes[i])

                  field = np.empty((nbNodes,fieldSizes),dtype=np.float)
                  data = struct.unpack(dataformat*(nbNodes*fieldSizes[i]), data)
                  #print(data)
                  cpt =0
                  for n in xrange(nbNodes):
                      for j in xrange(fieldSizes[i]):
                          field[n,j] = data[cpt]
                          cpt += 1

                  res["SolAtVertices"+str(i)]  =  field
              continue
          #End key
          if key == 54:
             break

          self.PrintVerbose("skiping key : "  + str(key) )
          data = f.read(4)
          endOfInformation = struct.unpack("i", data)[0]
          f.seek(endOfInformation)
      return res

    def ReadMeshBinary(self):
      self.readFormat = 'rb'
      self.StartReading()


      f = self.filePointer
      res = UM.UnstructuredMesh()
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
                      ref = int(udata[dimension])
                      nodalrefs[cpt] = ref
                      #print (str(cpt) +' : '+ str(ref) + ' | '),
                  else:
                      ref = str(udata[dimension])
                      res.GetNodalTag("NTag"+ref).AddToTag(cpt)

                  cpt +=1
                  self.PrintProgress(cpt,maxx = nbNodes)

                  #if (cpt/10000 == float(cpt)/10000):
                    #self.PrintProgress()
                  #  print(100*float(cpt)/nbNodes)
                  if cpt == nbNodes:
                    break
              #print(f.tell())
              del cpt
              continue



          # all kind of elements
          if meditNumber.has_key(key):
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

                self.PrintProgress(elements.cpt, maxx = nbElements )
                #if (elements.cpt/10000 == float(elements.cpt)/10000):
                #    print(elements.cpt)


              elements.connectivity -= 1;
              continue

          if key == 62:
             endOfInformation = struct.unpack("i", data)[0]
             f.seek(endOfInformation)
             break


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
         for name,val in res.elements.iteritems():
             elemRefs[cpt:cpt+val.GetNumberOfElements()] = elemRefsDic[name]
             cpt += val.GetNumberOfElements()
         self.elementsFields['refs'] = elemRefs

      self.EndReading()
      return res





    def ReadMeshAscii(self, _fileName=None,_string=None, fieldFileName= None):


        if _fileName is not None:
            self.SetFileName(_fileName)

        if _string is not None:
            self.SetStringToRead(_string)

        self.readFormat = 'r'
        self.StartReading()


        res = UM.UnstructuredMesh()
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
                #print ('Dimension : '+ str(dimension))
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

                    res.nodes[cpt,:] = map(float,s[0:dimension])
                    res.originalIDNodes[int(cpt)] = int(cpt)

                    ref = s[dimension]

                    res.GetNodalTag("NTag"+ref).AddToTag(cpt)

                    cpt +=1
                    if cpt == nbNodes:
                        break
                continue


            if meditName.has_key(l):
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


    __teststring = """
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


    __teststringField="""
MeshVersionFormatted
2

Dimension
3

SolAtVertices
4
1.
2.
3.
4.

End
"""

    res = ReadMesh(string=__teststring)
    print(res)

    from OTTools.Helpers.Tests import TestTempDir
    newFileName = TestTempDir().GetTempPath()+"mshFile.mesh"
    open(newFileName,'w').write(__teststring)
    res = ReadMesh(fileName=newFileName)

    newFileName = TestTempDir().GetTempPath()+"mshFile.sol"
    open(newFileName,'w').write(__teststringField)
    resfield = ReadMesh(fileName=newFileName)


    from OTTools.IO.MeshWriter import WriteMesh as WriteMesh

    newFileName = TestTempDir().GetTempPath()+"mshFile.meshb"
    WriteMesh(newFileName, res,binary=True)

    res = ReadMesh(newFileName)
    print(res)

    return 'ok'

if __name__ == '__main__':
    from OTTools.Helpers.BaseOutputObject import BaseOutputObject
    BaseOutputObject.SetGlobalDebugMode(True)
    print(CheckIntegrity())# pragma: no cover
