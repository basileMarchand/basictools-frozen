# -*- coding: utf-8 -*-
import struct
import numpy as np

import BasicTools.Containers.UnstructuredMesh  as UM
import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase

from BasicTools.IO.MeshTools import ASCIITypes
from BasicTools.IO.MeshTools import ASCIITags
from BasicTools.IO.MeshTools import BinaryTypes
from BasicTools.IO.MeshTools import BinaryKeywords as BKeys
from BasicTools.IO.MeshTools import BinaryTags
from BasicTools.IO.MeshTools import BinaryFields
import BasicTools.IO.MeshTools as MT

def ReadMesh(fileName=None,string=None,ReadRefsAsField=False ):
    reader = MeshReader()
    reader.SetReadRefsAsField(ReadRefsAsField)
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read()
    return reader.output

def ReadSol(fileName, out=None):
    reader = MeshSolutionReaderWrapper()
    reader.SetFileName(fileName)
    reader.Read(out=out)
    return reader.output

    import os.path

    dirname = os.path.dirname(fileName)
    basename,extention = os.path.splitext(os.path.basename(fileName))

    if extention[-1] == "b":
        f = os.path.join(dirname,basename+".meshb")
    else:
        f = os.path.join(dirname,basename+".mesh")


# we check if the file exist, if not we try the other type
    if not os.path.isfile(f):
        if extention[-1] == "b":
            f = os.path.join(dirname,basename+".mesh")
        else:
            f = os.path.join(dirname,basename+".meshb")

    if not os.path.isfile(f):
        raise Exception("unable to find a mesh file")
    reader = MeshReader()
    reader.SetFileName(fileName=f)
    reader.Read()
    mesh = reader.output
    fields = reader.ReadExtraFields(fileName);
    mesh.nodeFields = {k:v for k,v in fields.items() if k.find("SolAtVertices") != -1  }
    if 'SolAtTetrahedra0' in fields:

        if mesh.GetElementsOfType(EN.Tetrahedron_4).GetNumberOfElements() == mesh.GetNumberOfElements():
            mesh.elemFields = {k:v for k,v in fields.items() if k.find("SolAtTetrahedra") != -1  }
    return mesh


class MeshReader(ReaderBase):
    def __init__(self):
        super(MeshReader,self).__init__()
        self.refsAsAField = False
        self.dim = 3
        self.version = -1

    def SetReadRefsAsField(self, val):
        self.refsAsAField = bool(val)

    def Read(self,out=None):
        if self.binary :
            return self.ReadMeshBinary(out=out)
        else:
            return self.ReadMeshAscii(out=out)

    def SetFileName(self,fileName):
        super(MeshReader,self).SetFileName(fileName)

        if fileName is not None and fileName[-1] == "b":
            self.SetBinary(True)
        else:
            self.SetBinary(False)

    def ReadExtraFields(self,fileName):
        self.SetFileName(fileName)

        mesh = self.Read()

        res = {}
        res.update(mesh.nodeFields)
        res.update(mesh.elemFields)
        return res

#        if  fileName[-1] == "b":
#            return self.ReadExtraFieldsBinary(fileName)
#        else:
#            return self.ReadExtraFieldsAscii(fileName)
##ASCII PART #################################################################
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

            if len(l) > 7 and l.find("Vertices") == 0 :
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


            if l in ASCIITypes:
                elements = res.GetElementsOfType(ASCIITypes[l])
                nbNodes = elements.GetNumberOfNodesPerElement()
                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()

                nbElements = int(l.split()[0])
                if nbElements == 0:
                    continue
                self.PrintVerbose("Reading "+str(nbElements)+ " Elements")
                cpt =0;
                while(True):
                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0:
                        continue

                    s = list(map(int,l.split()))
                    elements.AddNewElement(s[0:nbNodes], cpt)
                    ref = s[nbNodes]
                    elements.GetTag("ETag"+str(ref)).AddToTag(cpt)

                    cpt +=1
                    if nbElements == cpt:# pragma: no cover
                        break

                elements.connectivity -= 1;
                continue

            if l in ASCIITags:
                elements = res.GetElementsOfType(ASCIITags[l][0])
                tag = elements.GetTag(ASCIITags[l][1])

                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()
                nbIds = int(l.split()[0])
                self.PrintVerbose("Reading tags Elements")
                cpt=0
                ids = []
                while(True):
                    if nbIds == cpt:# pragma: no cover
                        break

                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()

                    if len(l) == 0:
                        continue

                    newids = list(map(int,l.split()))
                    cpt += len(newids)

                    ids.extend(newids)

                tag.SetIds(np.array(ids,dtype=int)-1)
                continue

            if l in  ["SolAtVertices"]:
                fieldname = l
                data = self._ReadFieldsASCII(self,dimension)
                for i in range(len(data)) :
                    res.nodeFields[fieldname+str(i)]  =  data[i]
                continue

            if l in  ["SolAtTetrahedra" ]:
                fieldname = l
                data = self._ReadFieldsASCII(self,dimension)
                for i in range(len(data)) :
                    res.elemFields[fieldname+str(i)]  =  data[i]
                continue

            if l.find("End")>-1 :
                break

            self.PrintVerbose("ignoring line :->" + l +"<-")
        res.PrepareForOutput()
        self.EndReading()
        return res

    def _ReadFieldsASCII(self,myFile,dim):
        datares = []
        line = myFile.ReadCleanLine(withError = True)

        nbentries = int(line.split()[0])
        line = myFile.ReadCleanLine(withError = True)
        nbfields = int(line.split()[0])
        fieldSizes = [ int(x) for x in line.split()[1:] ]

        for i in range(len(fieldSizes)):
            if fieldSizes[i] == 2 :
               fieldSizes[i] = dim
            elif fieldSizes[i] == 3 :
                raise

        ncoomp = np.sum(fieldSizes)

        data = np.fromfile(file=myFile.filePointer,dtype=float,count=int(ncoomp*nbentries),sep=" \n")

        data.shape = (nbentries,ncoomp)
        cpt = 0
        for i in range(nbfields):
            datares.append(data[:,cpt:cpt+fieldSizes[i]])
            cpt += fieldSizes[i]
        return datares

##BINARY PART ################################################################################
    def _ReadBinaryInt(self):
        if self.version <= 3 :
            return self.readInt32()
        else:
            return self.readInt64()


    def _ReadBinaryHeader(self):
      key = self.readInt32()

      dataType = np.float32
      dataSize = 4;

      if key == BKeys["GmfVersionFormatted"]:
          self.version = self.readInt32()
          if self.version >= 2 :
              dataType = np.float64
              dataSize = 8
              #print("data in double")
      else:
          raise Exception('Expected key value equal to '+str(BKeys["GmfVersionFormatted"])+' (GmfVersionFormatted), got : '+str(key)+' are you sure this file is binary??')

      key = self.readInt32()
      if key == BKeys["GmfDimension"]:
          if self.version == 3:
              endOfInformation = self.readInt64()
          else:
              endOfInformation = self.readInt32()

          dimension = self.readInt32()
          self.filePointer.seek(endOfInformation)
      else:
          raise Exception('Expected key value equal to '+str(BKeys["GmfDimension"])+' (GmfDimension), got : '+str(key)+' are you sure this file is binary??')

      return dataType,dataSize,dimension

    def ReadMeshBinary(self,out=None):
      self.readFormat = 'rb'
      self.StartReading()

      f = self.filePointer

      if out is None:
          res = UM.UnstructuredMesh()
      else:
          res = out

      self.output = res

      dataType,dataSize,dimension = self._ReadBinaryHeader()

      globalElementCounter = 0;
      elemRefsDic = {}

      while True:
          key = self.readInt32()
          self.PrintVerbose("key" + str(key))
          if key == BKeys["GmfEnd"]:
              break

          if self.version == 3:
               endOfInformation = self.readInt64()
          else:
               endOfInformation = self.readInt32()
          #Vertices
          if key == BKeys["GmfVertices"]:
              nbNodes = self.readInt32()

              self.PrintVerbose("Reading " + str(nbNodes) + " nodes ")

              res.nodes = np.empty((nbNodes,dimension),dtype=dataType)
              res.originalIDNodes= np.empty((nbNodes,),dtype=np.int)

              if dataSize == 4:
                 dt =  np.dtype([('pos', np.float32,(dimension,) ), ('ref', np.int32, (1,))])
              else:
                 dt =  np.dtype([('pos', np.float64,(dimension,) ), ('ref', np.int32, (1,))])

              data = np.fromfile(self.filePointer,dtype=dt,count=nbNodes, sep="")

              res.nodes = data[:]["pos"]
              res.originalIDNodes = np.arange(nbNodes)

              refs = data[:]["ref"]
              if self.refsAsAField:
                  res.nodeFields['refs']  = refs
              else:
                  cpt =0
                  while(cpt < len(refs)):
                      ref = refs[cpt][0]
                      #print(ref)
                      res.GetNodalTag("NTag"+str(ref)).AddToTag(cpt)
                      cpt +=1
              continue

          if key == BKeys["GmfCorners"]:
              nbCorners = self.readInt32()
              data = np.fromfile(self.filePointer,dtype=np.int32,count=nbCorners, sep="")
              res.nodesTags.CreateTag("Corners").SetIds(data-1)
              continue

          # all kind of elements
          if key in BinaryTypes:

              elements = res.GetElementsOfType(BinaryTypes[key])
              self.PrintVerbose("Reading elements of type " + elements.elementType )
              nbNodes = elements.GetNumberOfNodesPerElement()
              nbElements = self.readInt32()

              elements.Reserve(nbElements)

              self.PrintVerbose('Reading '+ str(nbElements) + " Elements ")
              elements.cpt = 0;

              dt =  np.dtype([('conn', np.int32,(nbNodes,) ), ('ref', np.int32, (1,))])

              data = np.fromfile(self.filePointer,dtype=dt,count=nbElements, sep="")

              elements.connectivity = (data[:]["conn"]-1).astype(np.int_)

              elements.originalIds = np.arange(globalElementCounter,globalElementCounter+nbElements);
              elements.cpt = nbElements
              globalElementCounter +=nbElements


              refs = data[:]["ref"]
              if self.refsAsAField:
                  elemRefsDic[elements.elementType] = refs
              else:
                  cpt =0
                  for ref in refs:
                      #print(ref[0])
                      elements.GetTag("ETag"+str(ref[0])).AddToTag(cpt)
                      cpt +=1
              continue

              self.PrintVerbose('Reading '+ str(nbElements) + " Elements Done ")


              continue

          if key == BKeys["GmfRequiredVertices"]:
              tagname = MT.RequiredVertices
              self.PrintVerbose("Reading " + str(tagname) )
              nbentries = self.readInt32()
              ids = np.fromfile(self.filePointer,dtype=np.int32,count=nbentries, sep="")-1
              res.nodesTags.CreateTag(tagname).SetIds(ids)
              continue

          if key in BinaryTags:
              elemtype,tagname = BinaryTags[key]
              self.PrintVerbose("Reading " + str(tagname) )
              nbentries = self.readInt32()
              ids = np.fromfile(self.filePointer,dtype=np.int32,count=nbentries, sep="")-1
              elements = res.GetElementsOfType(elemtype)
              elements.tags.CreateTag(tagname).SetIds(ids)
              continue

          if key in BinaryFields:
              data = self._readExtraFieldBinary(self,dataSize,dimension)
              for i in range(len(data)) :
                  res.elemFields[BinaryFields[key]+str(i)]  =  data[i]
              continue

          if key ==  BKeys["GmfSolAtVertices"]:
             data = self._readExtraFieldBinary(self,dataSize,dimension)
             for i in range(len(data)) :
                  res.nodeFields["SolAtVertices"+str(i)]  =  data[i]
             continue


          if key not in [BKeys[x] for x in ["GmfNormals", "GmfNormalAtVertices", "GmfTangents", "GmfTangentAtVertices"] ]:
              self.PrintVerbose("skiping key : "  + str(key) )
          f.seek(endOfInformation)


      res.GenerateManufacturedOriginalIDs()

      if self.refsAsAField:
         elemRefs = np.empty(globalElementCounter,dtype=np.int )
         cpt =0
         for name,val in res.elements.items():
             elemRefs[cpt:cpt+val.GetNumberOfElements()] = elemRefsDic[name].ravel()
             cpt += val.GetNumberOfElements()
         res.elemFields['refs'] = elemRefs

      self.EndReading()
      return res



    def _readExtraFieldBinary(self,myFile,dataSize,dim):

      nbEntities = myFile._ReadBinaryInt()
      nbfields = myFile._ReadBinaryInt()

      fieldSizes = np.empty(nbfields, dtype=np.int)
      res = []
      for i in range(nbfields):
          fieldSizes[i] =  myFile._ReadBinaryInt()
          if fieldSizes[i] == 2 :
             fieldSizes[i] = dim
          elif fieldSizes[i] == 3 :
             raise

      ncoomp = np.sum(fieldSizes)
      if dataSize == 4:
          data = myFile.readFloats32(nbEntities*ncoomp)
      else:
          data = myFile.readFloats64(nbEntities*ncoomp)

      data.shape = (nbEntities,ncoomp)

      cpt =0
      for i in range(nbfields):
          res.append(data[:,cpt:cpt+fieldSizes[i]])
          cpt += fieldSizes[i]
      return res


#    def ReadExtraFieldsBinary(self,fileName):
##      self.PrintVerbose("Reading ExtraField")
#
#      myFile = ReaderBase()
#      myFile.SetFileName(fileName)
#      myFile.readFormat = "rb"
#      myFile.StartReading()
#
#      self._ReadBinaryHeader()
#
#      dataType = np.float32
#      dataSize = 4;
#
#      self.ReadHeader()
#      key = myFile.readInt32()
#
#     #MeshVersionFormatted
#      if key == 1:
#          udata = myFile.readInt32()
#          if udata == 2 :
#              dataType = np.float64
#              dataSize = 8
#              #print("data in double")
#      else:
#          raise Exception('Expected key value equal to 1 (binary), are you sure this file is binary??')
#
#      #if self.refsAsAField:
#      #   elemRefs = np.zeros(0 ,dtype=np.int)
#      res = {}
#      KeyNumberName = {}
#      KeyNumberName[62] =  "SolAtVertices"
#      KeyNumberName[64] =  "SolAtTriangles"
#      KeyNumberName[66] =  "SolAtTetrahedra"
#      while True:
#          key = myFile.readInt32()
#          print("Key " +str(key))
#          #End key
#          if key == 54:
#             break
#          endOfInformation = myFile.readInt32()
#
#          print("End of information " + str(endOfInformation) )
#          if key == 3:
#              #f.seek(endOfInformation-4)
#              dimension = myFile.readInt32()
#              myFile.seek(endOfInformation)
#              #print("dimension : " + str(dimension))
#              continue
#
#          #SolAtVertices
#          if key in KeyNumberName:
#              data = self._readExtraFieldBinary(myFile,dataSize,dimension)
#              for i in range(len(data)) :
#                  res[KeyNumberName[key]+str(i)]  =  data[i]
#              continue
#
#
#          print("skiping key : "  + str(key) )
#          myFile.seek(endOfInformation)
#      myFile.EndReading()
#      return res



#    def ReadExtraFieldsAscii(self,fileName):
#        #self.PrintVerbose("Reading ExtraField")
#
#        myFile = ReaderBase()
#        myFile.SetFileName(fileName)
#        print(myFile.readFormat)
#        myFile.StartReading()
#        res = {}
#
#        dataType = float
#        while(True):
#            line = myFile.ReadCleanLine()
#            print("-> "+line)
#            if line is None :
#                break
#
#            if line.find("MeshVersionFormatted")>-1 :
#                if len(line.split()) > 1:
#                    formatFile = int(line.split()[1])
#                else:
#                    line = myFile.ReadCleanLine(withError = True)
#                    formatFile = int(line)
#
#                if formatFile == 2:
#                    dataType = np.float64
#                #print ('formatFile  : '+ str(dataType))
#                continue
#
#            if line.find("Dimension")>-1 :
#                s = line.split()
#                if len(s)>1:
#                    s.pop(0)
#                    strdim = " ".join(s)
#                else:
#                    strdim = myFile.ReadCleanLine(withError = True)
#                dim = int(strdim)
#                print ('Dimension : '+ str(dim))
#                continue
#                print ('Dimension : '+ str(dim))
#
#
#
#            for fieldName in  ["SolAtVertices", "SolAtTetrahedra" ]:
#                found = False
#                if line.find(fieldName)>-1 and len(line) == len(fieldName):
#                    #print("reading SolAtVertices")
#                    data = self._ReadFieldsASCII(myFile,dim)
#                    for i in range(len(data)) :
#                        res[fieldName+str(i)]  =  data[i]
#                    found = True
#                    break
#
#            if found:
#                continue
#
#            if line.find("End")>-1 :
#                break
#
#            print("Ignoring entrie " + str(line) )
#            raise
#
#        myFile.EndReading()
#        return res
class MeshSolutionReaderWrapper():
    def __init__(self):
       super(MeshSolutionReaderWrapper,self).__init__()
       self.canHandleTemporal = False

    def SetFileName(self,fileName):
        import os.path
        self.fileName = fileName
        dirname = os.path.dirname(fileName)
        basename,extention = os.path.splitext(os.path.basename(fileName))

        if extention[-1] == "b":
            f = os.path.join(dirname,basename+".meshb")
        else:
            f = os.path.join(dirname,basename+".mesh")


        # we check if the file exist, if not we try the other type
        if not os.path.isfile(f):
            if extention[-1] == "b":
                f = os.path.join(dirname,basename+".mesh")
            else:
                f = os.path.join(dirname,basename+".meshb")

        if not os.path.isfile(f):
            raise Exception("unable to find a mesh file")
        self.reader = MeshReader()
        self.reader.SetFileName(fileName=f)

    def Read(self,out=None):
        if out :
            raise

        self.reader.Read()
        mesh = self.reader.output
        fields = self.reader.ReadExtraFields(self.fileName);
        mesh.nodeFields = {k:v for k,v in fields.items() if k.find("SolAtVertices") != -1  }
        if 'SolAtTetrahedra0' in fields:
            if mesh.GetElementsOfType(EN.Tetrahedron_4).GetNumberOfElements() == mesh.GetNumberOfElements():
                mesh.elemFields = {k:v for k,v in fields.items() if k.find("SolAtTetrahedra") != -1  }
        return mesh




from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".mesh",MeshReader)
RegisterReaderClass(".meshb",MeshReader)

RegisterReaderClass(".sol",MeshSolutionReaderWrapper)
RegisterReaderClass(".solb",MeshSolutionReaderWrapper)

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
    f = open(newFileName,'w')
    f.write(__teststringField)
    f.close()

    print('Reading : ' +str(newFileName))

    resfield = MeshReader().ReadExtraFields(fileName=newFileName)


    from BasicTools.IO.MeshWriter import WriteMesh as WriteMesh

    print("###########################################################")
    print(res)
    newFileName = TestTempDir().GetTempPath()+"mshFile.meshb"
    WriteMesh(newFileName, res,binary=True)

    res = ReadMesh(newFileName,ReadRefsAsField=True)
    print(res)


    sol = MeshReader().ReadExtraFields(TestTempDir().GetTempPath()+"mshFile.sol")
    print(sol)

    from BasicTools.IO.MeshWriter import MeshWriter as MeshWriter
    mw = MeshWriter()
    mw.SetBinary(True)
    mw.SetFileName(TestTempDir().GetTempPath()+"mshFile.solb")
    mw.OpenSolutionFile(res)

    mw.WriteSolutionsFields(res,PointFields=[sol['SolAtVertices0']])
    mw.filePointer.write(struct.pack('i', 54))
    mw.Close()

    sol = MeshReader().ReadExtraFields(TestTempDir().GetTempPath()+"mshFile.solb")

    return 'ok'

if __name__ == '__main__':# pragma: no cover
    from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
    BaseOutputObject.SetGlobalDebugMode(True)
    print(CheckIntegrity())# pragma: no cover
