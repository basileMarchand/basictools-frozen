# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import struct
import numpy as np

from BasicTools.Containers.MeshBase import Tag as Tag
from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.Containers.ElementNames as EN

from BasicTools.IO.MeshTools import BinaryKeywords,BinaryNumber,ASCIIName, ASCIITags
from BasicTools.IO.MeshTools import Corners,Ridges,RequiredEdges,RequiredTriangles,RequiredVertices
import BasicTools.IO.MeshTools as MT
from BasicTools.IO.MeshTools import BinaryKeywords as BKeys
from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType

def WriteMesh(filename,mesh,PointFields=None,solutionOnOwnFile= False, binary=True, nodalRefNumber= None,elemRefNumber=None):
    OW = MeshWriter()
    OW.SetBinary(binary);
    OW.Open(filename)
    OW.Write(mesh,PointFields = PointFields,solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber,elemRefNumber=elemRefNumber)
    OW.Close()

class MeshWriter(WriterBase):
    def __init__(self):
        super(MeshWriter,self).__init__()
        self.dataType = np.float32
        self.dataSize = 4;
        self.SetSinglePrecission()

    def __str__(self):
        res  = 'MeshWriter : \n'
        res += '   FileName : '+str(self.fileName) + '\n'
        res += '   Binary : ' + ('True' if self.isBinary() else 'False') + "\n"
        return res

    def SetFileName(self,fileName):
        self.fileName = fileName;

    def SetSinglePrecission(self,single=True):
        if single:
            self.dataType = np.float32
            self.dataSize = 4;
        else:
            self.dataType = np.float64
            self.dataSize = 8;

    def Write(self,meshObject,PointFields=None, solutionOnOwnFile= False, nodalRefNumber= None, elemRefNumber=None,PointFieldsNames=None,CellFieldsNames=None,CellFields=None):
        if self.isBinary():
            return self.WriteBINARY(meshObject,PointFields=PointFields, solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber,elemRefNumber=elemRefNumber )
        else:
            return self.WriteASCII(meshObject,PointFields=PointFields, solutionOnOwnFile= solutionOnOwnFile, nodalRefNumber =nodalRefNumber,elemRefNumber=elemRefNumber )

    def WriteBINARY(self,meshObject,PointFields=None, solutionOnOwnFile= False, nodalRefNumber = None,elemRefNumber=None):
        self.SetSinglePrecission(single=False)
        #key MeshVersionFormatted
        self.filePointer.write(struct.pack('i', 1))
        self.filePointer.write(struct.pack('i', self.dataSize//4))

        #key Dimension (3)
        self.filePointer.write(struct.pack('i', 3))

        flat= True
        if meshObject.nodes.shape[1] == 3:
            mmax = np.max(meshObject.nodes[:,2])
            mmin = np.min(meshObject.nodes[:,2])
            flat = (mmax == mmin)
        if meshObject.nodes.shape[1] == 3 and  not flat :
            dimension = 3
        else:
            dimension = 2

        self.filePointer.write(struct.pack('i', self.filePointer.tell()+8))# end of information
        self.filePointer.write(struct.pack('i', dimension)) #dimension

        #key Vertices (4)
        self.filePointer.write(struct.pack('i', 4))
        currentposition = self.filePointer.tell()
        numberofpoints = meshObject.GetNumberOfNodes()
        endOfInformation = currentposition + 2*4+ numberofpoints*(self.dataSize*dimension+4)
        self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
        self.filePointer.write(struct.pack('i', numberofpoints))# numberofpoints

        dataformat = ('f' if self.dataSize == 4 else 'd')
        posn = meshObject.GetPosOfNodes()
        for n in range(numberofpoints):
             posn[n,0:dimension].astype(self.dataType).tofile(self.filePointer, sep='')
             if nodalRefNumber is  None :
                 self.filePointer.write(struct.pack("i", 0))# refs
             else :
                 self.filePointer.write(struct.pack("i", nodalRefNumber[n]))# refs

        if MT.RequiredVertices  in meshObject.nodesTags :
            ids = meshObject.nodesTags[MT.RequiredVertices].GetIds()+1
            nbids = len(ids)
            if nbids:
               self.filePointer.write(struct.pack('i', 15 ))
               currentposition = self.filePointer.tell()
               endOfInformation = currentposition+ (2+len(ids))*4
               self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
               self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
               ids.astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')

        globalOffset =0
        for elementContainer in meshObject.elements:
            self.PrintDebug("Output of " + str(elementContainer ))
            data = meshObject.elements[elementContainer]
            nbelements = data.GetNumberOfElements()
            if nbelements == 0:
                self.Print("Empty elemnt container (skipping) + " + elementContainer)
                continue

            if elementContainer == EN.Point_1:
               print("MeshWriter warning: ignoring EN.Point_1 elements ")
               globalOffset += data.GetNumberOfElements()
               continue 

            elemtype = BinaryNumber[elementContainer]

            self.filePointer.write(struct.pack('i', elemtype ))
            nbNodes = data.GetNumberOfNodesPerElement()

            currentposition = self.filePointer.tell()
            endOfInformation = currentposition+ (2+data.GetNumberOfElements()*(nbNodes+1))*4
            self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
            self.filePointer.write(struct.pack('i', nbelements))# GetNumberOfElements

            tempcoon =  np.zeros( (data.connectivity.shape[0],data.connectivity.shape[1]+1),dtype=np.int32,order="c"  )
            tempcoon[:,0:data.connectivity.shape[1]] = data.connectivity.astype(np.int32)
            tempcoon[:,0:data.connectivity.shape[1]] += 1
            if elemRefNumber is not None:
                tempcoon[:,data.connectivity.shape[1]] = elemRefNumber[globalOffset: globalOffset+nbelements]
            tempcoon.tofile(self.filePointer, format=dataformat,sep='')


            globalOffset += data.GetNumberOfElements()

        if "Corners" in meshObject.nodesTags:
            tag = meshObject.nodesTags['Corners']
            ids = tag.GetIds()
            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', BKeys["GmfCorners"] ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')

        bars = meshObject.GetElementsOfType(EN.Bar_2)
        if "Ridges" in bars.tags  and len(bars.tags["Ridges"]):

            self.PrintDebug("output Ridges" )
            ids = np.empty(0,dtype=int)
            if "Ridges" in bars.tags :
                tag = bars.tags["Ridges"]

                ids = np.concatenate((ids,tag.GetIds()))

            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', 14 ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')

        bars = meshObject.GetElementsOfType(EN.Bar_2)
        if RequiredEdges in bars.tags  and len(bars.tags[RequiredEdges]):

            self.PrintDebug("output RequiredEdges" )
            ids = np.empty(0,dtype=int)
            if RequiredEdges in bars.tags :
                tag = bars.tags[RequiredEdges]

                ids = np.concatenate((ids,tag.GetIds()))
            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', BinaryKeywords["GmfRequiredEdges"] ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')

        tris = meshObject.GetElementsOfType(EN.Triangle_3)
        if "RequiredTriangles" in tris.tags  and len(tris.tags["RequiredTriangles"]):

            self.PrintDebug("output RequiredTriangles" )
            ids = np.empty(0,dtype=int)
            if "RequiredTriangles" in tris.tags :
                tag = tris.tags["RequiredTriangles"]

                ids = np.concatenate((ids,tag.GetIds()))
            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', 17 ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')

        if PointFields is not None and len(PointFields)>0 :
            if solutionOnOwnFile :
                self.Close();
                self.OpenSolutionFileBinary(support=meshObject)
            self.WriteSolutionsFieldsBinary(meshObject,PointFields)

    def OpenSolutionFileAscii(self,support):
        self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol" , 'w')
        self._isOpen = True

        self.filePointer.write("MeshVersionFormatted 2\n\n")
        dimension = 3

        if support.nodes.shape[1] == 3 :
            mmax = np.max(support.nodes[:,2])
            mmin = np.min(support.nodes[:,2])
            if  mmax ==  mmin and mmax == 0.:
               dimension = 2
        self.filePointer.write("Dimension {}\n\n".format(dimension))

    def WriteSolutionsFieldsAscii(self,meshObject,PointFields=None,SolsAtTriangles=None,SolsAtTetrahedra=None):
        if PointFields is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtVertices",PointFields)

        if SolsAtTriangles is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtTriangles",SolsAtTriangles)

        if SolsAtTetrahedra is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtTetrahedra",SolsAtTetrahedra)

    def _WriteSolutionsFieldsAsciiUsingKey(self,meshObject,keyword,Sols):
        if len(Sols) == 0:
            return
        nbentries = Sols[0].shape[0]

        #key SolAtVertices = 62
        self.filePointer.write("{}\n".format(keyword))

        self.filePointer.write("{}\n".format(nbentries))

        nbfields = len(Sols)

        self.filePointer.write("{} ".format(nbfields))

        for sol in Sols:
            if len(sol.shape)== 1:
                sol = sol[:,np.newaxis]
            size = sol.shape[1]

            if size == 1:
                self.filePointer.write("1 ")
            elif size == meshObject.GetDimensionality():
                self.filePointer.write("2 ")
            else:
                self.filePointer.write("3 ")

        self.filePointer.write("\n\n")

        for i in range(nbentries):
            for sol in Sols:
                if len(sol.shape)== 1:
                    sol[i].tofile(self.filePointer, sep=' ' )
                else:
                    sol[i,:].tofile(self.filePointer, sep=' ' )
                self.filePointer.write(" ")
            self.filePointer.write("\n")
        self.filePointer.write("\n")

    def CloseSolutionFileAscii(self):
        self.filePointer.write("End\n") #dimension

    def OpenSolutionFile(self,support):
        if self.isBinary():
            self.OpenSolutionFileBinary(support)
        else:
            self.OpenSolutionFileAscii(support)

    def WriteSolutionsFields(self,meshObject,PointFields=None,SolsAtTriangles=None,SolsAtTetrahedra=None):
        if self.isBinary():
            self.WriteSolutionsFieldsBinary(meshObject,PointFields=PointFields,SolsAtTriangles=SolsAtTriangles,SolsAtTetrahedra=SolsAtTetrahedra)
        else:
            self.WriteSolutionsFieldsAscii(meshObject,PointFields=PointFields,SolsAtTriangles=SolsAtTriangles,SolsAtTetrahedra=SolsAtTetrahedra)

    def Close(self):

        if self.isBinary():
            self.CloseSolutionFileBinary()
        else:
            self.CloseSolutionFileAscii()
        WriterBase.Close(self)

    def OpenSolutionFileBinary(self,support):

        self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".solb" , 'wb',0)
        self._isOpen = True

        #key MeshVersionFormatted
        self.filePointer.write(struct.pack('i', 1))
        self.filePointer.write(struct.pack('i', self.dataSize//4))
        #
        #key Dimension (3)
        self.filePointer.write(struct.pack('i', 3))
        dimension = 3 #support.GetDimensionality()

        if support.nodes.shape[1] == 3 :
            mmax = np.max(support.nodes[:,2])
            mmin = np.min(support.nodes[:,2])
            if  mmax ==  mmin:
               dimension = 2

        self.filePointer.write(struct.pack('i', self.filePointer.tell()+4*2))# end of information
        self.filePointer.write(struct.pack('i', dimension)) #dimension


    def WriteSolutionsFieldsBinary(self,meshObject,PointFields=None,SolsAtTriangles=None,SolsAtTetrahedra=None):
        if PointFields is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject,62,PointFields)

        if SolsAtTriangles is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject,64,SolsAtTriangles)

        if SolsAtTetrahedra is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject,66,SolsAtTetrahedra)


    def _WriteSolutionsFieldsBinaryUsingKey(self,meshObject,key,Sols):

        if len(Sols) == 0:
            return

        NumberOfEntries = Sols[0].shape[0]

        self.filePointer.write(struct.pack('i', key))
        nbfields = len(Sols)

        nbcoomp = 0
        for sol in Sols:
            if len(sol.shape) == 1:
                nbcoomp += 1
            else:
                nbcoomp += sol.shape[-1]

        endOfInformation = self.filePointer.tell()+4*(1+1+1+nbfields)+nbcoomp*NumberOfEntries*self.dataSize
        self.filePointer.write(struct.pack('i', endOfInformation))# end of information

        self.filePointer.write(struct.pack('i', NumberOfEntries))# numberofpoints
        self.filePointer.write(struct.pack('i',nbfields))# numberofpoints

        for sol in Sols:
            # we add a extra axis for scalar field stored in a vector (only one index)
            if len(sol.shape)== 1:
                sol = sol[:,np.newaxis]

            size = sol.shape[1]
            if size == 1:
                self.filePointer.write(struct.pack('i',1))#
            else:
                ## vectors
                self.filePointer.write(struct.pack('i',2))#

        for i in range(NumberOfEntries):
            for sol in Sols:
                if len(sol.shape)== 1:
                    sol = sol[:,np.newaxis]
                sol[i,:].astype(self.dataType).tofile(self.filePointer, sep='')

        if ( not (self.filePointer.tell() == endOfInformation) ) : raise Exception("Error in the writing code, please debug me!!!")

    def CloseSolutionFileBinary(self):
        self.filePointer.write(struct.pack('i', 54)) #dimension

    def WriteASCII(self,meshObject,PointFields=None, solutionOnOwnFile= False, nodalRefNumber = None,elemRefNumber=None):
        meshObject.ComputeGlobalOffset()
        self.filePointer.write("MeshVersionFormatted 2 \n")

        self.filePointer.write("Dimension " + str(meshObject.GetDimensionality()) + "\n\n")
        self.filePointer.write("Vertices\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} \n\n".format(numberofpoints) )
        posn = meshObject.GetPosOfNodes()

        for n in range(numberofpoints):
               posn[np.newaxis,n,:] .tofile(self.filePointer, sep=" ")
               if nodalRefNumber is None:
                   self.filePointer.write(" 0 \n")
               else:
                   self.filePointer.write(" {} \n".format(int(nodalRefNumber[n])) )
        self.filePointer.write("\n" )

        if meshObject.IsConstantRectilinear():
            import BasicTools.Containers.ElementNames
            elements = [ BasicTools.Containers.ElementNames.Hexaedron_8 ]
        else:
            elements = meshObject.elements

        globalOffset = 0
        for name, elementContainer in elements.items():
            if elementContainer.GetNumberOfElements() == 0:
                continue

            elemtype = ASCIIName.get(name,None)
            if elemtype is None:
                print("(MeshWriter) skiping this type of elements : " + name )
                continue

            if meshObject.IsConstantRectilinear():
                #hack to make the constantrectilinear api as the unstructured #TODO clean me
                class T(): pass
                data = T()
                data.connectivity = meshObject.GenerateFullConnectivity()
                nelem = meshObject.GetNumberOfElements()
            else:
                data = elementContainer
                nelem = data.GetNumberOfElements()
            if nelem == 0:
                continue
            self.filePointer.write("{} \n".format(elemtype) )
            self.filePointer.write("{}\n".format(nelem) )
            for i in range(nelem):
                    connectivity = (data.connectivity[i,:]+1).astype(int)
                    connectivity.tofile(self.filePointer, sep=" ")

                    if elemRefNumber is None :
                        self.filePointer.write(" 0\n")#
                    else:
                        self.filePointer.write(" " + str(elemRefNumber[globalOffset+i]) + "\n")#

            globalOffset += data.GetNumberOfElements()
            self.filePointer.write("\n")

        nTags = [RequiredVertices,Corners]

        for tagname in nTags:
            if tagname in meshObject.nodesTags:
               tag = meshObject.nodesTags[tagname]
               if len(tag):
                  self.filePointer.write(tagname+" \n");
                  self.filePointer.write("{} \n\n".format(len(tag)) )
                  (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
                  self.filePointer.write("\n" )

#        if "Corners" in meshObject.nodesTags:
#            tag = meshObject.nodesTags['Corners']
#            if len(tag):
#                self.filePointer.write("Corners\n");
#                self.filePointer.write("{} \n\n".format(len(tag)) )
#                (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
#                self.filePointer.write("\n" )

        for TagNameInFile,(ElementType,TagName) in ASCIITags.items():
            elements = meshObject.GetElementsOfType(ElementType)
            if TagName in elements.tags:
                tag = elements.tags[TagName ]
                if len(tag):
                    self.filePointer.write(str(TagNameInFile)+"\n");
                    self.filePointer.write("{} \n\n".format(len(tag)) )
                    (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
                    self.filePointer.write("\n" )

        self.filePointer.write("\n")

        if solutionOnOwnFile :
            self.Close();
            self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol" , 'w')
            self._isOpen = True
            self.filePointer.write("#Written by BasicTools package\n")
            self.filePointer.write("MeshVersionFormatted\n2 \n")
            self.filePointer.write("Dimension\n" + str(meshObject.GetDimensionality()) + "\n\n")

        if PointFields is not None and len(PointFields)>0 :
            if solutionOnOwnFile :
                self.Close();
                self.OpenSolutionFileAscii(support=meshObject)
            self.WriteSolutionsFieldsAscii(meshObject,PointFields)

from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".mesh",MeshWriter)

def CreateMeshWriterBinary(ops):
    obj = MeshWriter()
    obj.SetBinary()
    return obj
RegisterWriterClass(".meshb",MeshWriter,CreateMeshWriterBinary)

def CheckIntegrity():
    import BasicTools.Containers.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    with mymesh.WithModification():
        mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1]],dtype=PBasicFloatType)
        mymesh.originalIDNodes = np.array([1, 3, 4, 5,6],dtype=PBasicIndexType)

        mymesh.nodesTags.CreateTag("coucou").AddToTag(0)
        mymesh.nodesTags.CreateTag(Corners).AddToTag(0)
        mymesh.nodesTags.CreateTag(MT.RequiredVertices).SetIds([0])

        tets = mymesh.GetElementsOfType(EN.Tetrahedron_4)
        tets.AddNewElement([0,1,2,4],0)

        tris = mymesh.GetElementsOfType(EN.Triangle_3)
        tris.AddNewElement([0,1,2],0)
        tris.AddNewElement([2,1,3],3)
        tris.originalIds = np.array([3, 5],dtype=PBasicIndexType)

        tris.tags.CreateTag("RequiredTriangles").AddToTag(0)

        tris = mymesh.GetElementsOfType(EN.Bar_3)

        bars = mymesh.GetElementsOfType(EN.Bar_2)
        bars.AddNewElement([0,1],0)
        bars.tags.CreateTag("Ridges").AddToTag(0)
        bars.tags.CreateTag(MT.RequiredEdges).AddToTag(0)

    nodalRefNumber = np.arange(mymesh.GetNumberOfNodes())
    elemRefNumber = np.arange(mymesh.GetNumberOfElements())

    WriteMesh(tempdir+"CheckIntegrity_with_refs.mesh",mymesh,PointFields=[nodalRefNumber*10], elemRefNumber=elemRefNumber,nodalRefNumber=nodalRefNumber,solutionOnOwnFile=True)
    WriteMesh(tempdir+"CheckIntegrity_with_refs.mesh",mymesh,binary=False, PointFields=[nodalRefNumber*10], elemRefNumber=elemRefNumber,nodalRefNumber=nodalRefNumber,solutionOnOwnFile=True)

    OW = MeshWriter()
    OW.SetBinary(False)
    OW.SetGlobalDebugMode(True)
    OW.Open(tempdir+"Test_MmgWriter_ASCII.mesh")
    print(OW)
    OW.Write(mymesh)
    OW.Close()

    OWB = MeshWriter()
    OWB.SetBinary(True)
    OWB.Open(tempdir+"Test_MmgWriter_BIN.meshb")
    print(OWB)
    OWB.Write(mymesh)
    OWB.Close()

    print(mymesh)
    sol = np.arange(mymesh.GetNumberOfNodes(),dtype=PBasicFloatType)
    WriteMesh(tempdir+"Test_MmgWriter_II_Binary.mesh", mymesh ,PointFields=[sol ] )
    WriteMesh(tempdir+"Test_MmgWriter_II_Ascii.mesh", mymesh ,PointFields=[sol ], binary=False )

    res = CreateMeshWriterBinary({})
    print(res)

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
