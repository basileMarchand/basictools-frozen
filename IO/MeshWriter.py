# -*- coding: utf-8 -*-

import struct
import numpy as np
__author__ = "Felipe Bordeu"

from BasicTools.FE.MeshBase import Tag as Tag
from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.FE.ElementNames as EN


from MeshTools import BinaryNumber
from MeshTools import ASCIIName
from MeshTools import ASCIITags
from MeshTools import Corners


def WriteMesh(filename,mesh,SolsAtVertices=None,solutionOnOwnFile= False, binary=True, nodalRefNumber= None,elemRefNumber=None):
    OW = MeshWriter()
    OW.SetBinary(binary);
    OW.Open(filename)
    OW.Write(mesh,SolsAtVertices = SolsAtVertices,solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber,elemRefNumber=elemRefNumber)
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

    def Write(self,meshObject,SolsAtVertices=None, solutionOnOwnFile= False, nodalRefNumber= None, elemRefNumber=None):
        if self.isBinary():
            return self.WriteBINARY(meshObject,SolsAtVertices=SolsAtVertices, solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber,elemRefNumber=elemRefNumber )
        else:
            return self.WriteASCII(meshObject,SolsAtVertices=SolsAtVertices, solutionOnOwnFile= solutionOnOwnFile, nodalRefNumber =nodalRefNumber,elemRefNumber=elemRefNumber )



    def WriteBINARY(self,meshObject,SolsAtVertices=None, solutionOnOwnFile= False, nodalRefNumber = None,elemRefNumber=None):


        #key MeshVersionFormatted
        self.filePointer.write(struct.pack('i', 1))
        self.filePointer.write(struct.pack('i', self.dataSize//4))

        #key Dimension (3)
        self.filePointer.write(struct.pack('i', 3))
        dimension = meshObject.GetDimensionality()
        self.filePointer.write(struct.pack('i', self.filePointer.tell()+8))# end of information
        self.filePointer.write(struct.pack('i', dimension)) #dimension


        #key Vertices (4)
        self.filePointer.write(struct.pack('i', 4))
        currentposition = self.filePointer.tell()
        numberofpoints = meshObject.GetNumberOfNodes()
        #  endOfInformation + numberofpoints + (3 float + ref )
        endOfInformation = currentposition + 2*4+ numberofpoints*(self.dataSize*dimension+4)
        self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
        self.filePointer.write(struct.pack('i', numberofpoints))# numberofpoints

        dataformat = ('f' if self.dataSize == 4 else 'd')
        posn = meshObject.GetPosOfNodes()
        for n in range(numberofpoints):
             #self.filePointer.write(struct.pack(dataformat,  posn[n,:]))# end of information
             posn[n,:].astype(self.dataType).tofile(self.filePointer, sep='')

             if nodalRefNumber is  None :
                 self.filePointer.write(struct.pack("i", 0))# refs
             else :
                 self.filePointer.write(struct.pack("i", nodalRefNumber[n]))# refs



        self.PrintDebug("position at end " + str(self.filePointer.tell()))
        self.PrintDebug("calculate position at end " + str(endOfInformation))

        globalOffset =0
        for elementContainer in meshObject.elements:
            self.PrintDebug("Output of " + str(elementContainer ))
            elemtype = BinaryNumber[elementContainer]
            data = meshObject.elements[elementContainer]
            nbelements = data.GetNumberOfElements()
            if nbelements == 0:
                self.Print("Empty elemnt container (skipping) + " + elementContainer)
                continue



            self.filePointer.write(struct.pack('i', elemtype ))
            nbNodes = data.GetNumberOfNodesPerElement()

            currentposition = self.filePointer.tell()
            endOfInformation = currentposition+ (2+data.GetNumberOfElements()*(nbNodes+1))*4
            self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
            self.filePointer.write(struct.pack('i', nbelements))# GetNumberOfElements

            dataformat = "i"*nbNodes
            for i in range(nbelements):
                (data.connectivity[i,:]+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')
                if elemRefNumber is None :
                    self.filePointer.write(struct.pack("i", 0 ))#
                else:
                    self.filePointer.write(struct.pack("i", elemRefNumber[globalOffset+i] ))#

            globalOffset += data.GetNumberOfElements()

            self.PrintDebug("position at end " + str(self.filePointer.tell()))
            self.PrintDebug("calculate position at end " + str(endOfInformation))



        bars = meshObject.GetElementsOfType(EN.Bar_2)
        if "Ridges" in bars.tags  and len(bars.tags["Ridges"]):

            self.PrintDebug("output Ridges" )
            ids = np.empty(0,dtype=int)
            if "Ridges" in bars.tags :
                tag = bars.tags["Ridges"]

                ids = np.concatenate((ids,tag.GetIds()))

            #if "Meca" in bars.tags :
            #    tag = bars.tags["Meca"]
            #    ids = np.concatenate((ids,tag.GetIds()))

            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', 14 ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')
                self.PrintDebug("position at end " + str(self.filePointer.tell()))
                self.PrintDebug("calculate position at end " + str(endOfInformation))

        tris = meshObject.GetElementsOfType(EN.Triangle_3)
        if "RequiredTriangles" in tris.tags  and len(tris.tags["RequiredTriangles"]):

            self.PrintDebug("output RequiredTriangles" )
            ids = np.empty(0,dtype=int)
            if "RequiredTriangles" in tris.tags :
                tag = tris.tags["RequiredTriangles"]

                ids = np.concatenate((ids,tag.GetIds()))

            #if "Meca" in bars.tags :
            #    tag = bars.tags["Meca"]
            #    ids = np.concatenate((ids,tag.GetIds()))

            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', 17 ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')
                self.PrintDebug("position at end " + str(self.filePointer.tell()))
                self.PrintDebug("calculate position at end " + str(endOfInformation))


        if "Meca" in bars.tags :

            self.PrintDebug("output Mecas" )
            ids = np.empty(0,dtype=int)

            if "Meca" in bars.tags :
                tag = bars.tags["Meca"]
                ids = np.concatenate((ids,tag.GetIds()))

            nbids = len(ids)
            if nbids:
                self.filePointer.write(struct.pack('i', 16 ))
                currentposition = self.filePointer.tell()
                endOfInformation = currentposition+ (2+nbids)*4
                self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
                self.filePointer.write(struct.pack('i', nbids))# GetNumberOfElements
                (ids+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')
                self.PrintDebug("position at end " + str(self.filePointer.tell()))
                self.PrintDebug("calculate position at end " + str(endOfInformation))




        if SolsAtVertices is not None:
            if solutionOnOwnFile :
                #self.filePointer.write(struct.pack('i', 54)) #dimension
                self.Close();
                self.OpenSolutionFileBinary()
            self.WriteSolutionsFieldsBinary(meshObject,SolsAtVertices)

        #key End
        #self.filePointer.write(struct.pack('i', 54)) #dimension
        #self.Close()

    def OpenSolutionFileAscii(self,support):

        self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol" , 'wb',0)
        self._isOpen = True

        #key MeshVersionFormatted
        self.filePointer.write("MeshVersionFormatted 2\n\n")
        dimension = support.GetDimensionality()
        self.filePointer.write("Dimension {}\n\n".format(dimension))

    def WriteSolutionsFieldsAscii(self,meshObject,SolsAtVertices=None,SolsAtTriangles=None,SolsAtTetrahedra=None):
        if SolsAtVertices is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtVertices",SolsAtVertices)

        if SolsAtTriangles is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtTriangles",SolsAtTriangles)

        if SolsAtTetrahedra is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtTetrahedra",SolsAtTetrahedra)

    def _WriteSolutionsFieldsAsciiUsingKey(self,meshObject,keyword,Sols):
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
        #self._isOpen = False
        #self.filePointer.close()



    def OpenSolutionFile(self,support):
        if self.isBinary():
            self.OpenSolutionFileBinary(support)
        else:
            self.OpenSolutionFileAscii(support)

    def WriteSolutionsFields(self,meshObject,SolsAtVertices=None,SolsAtTriangles=None,SolsAtTetrahedra=None):
        if self.isBinary():
            self.WriteSolutionsFieldsBinary(meshObject,SolsAtVertices=SolsAtVertices,SolsAtTriangles=SolsAtTriangles,SolsAtTetrahedra=SolsAtTetrahedra)
        else:
            self.WriteSolutionsFieldsAscii(meshObject,SolsAtVertices=SolsAtVertices,SolsAtTriangles=SolsAtTriangles,SolsAtTetrahedra=SolsAtTetrahedra)

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
        dimension = support.GetDimensionality()
        self.filePointer.write(struct.pack('i', self.filePointer.tell()+4*2))# end of information
        self.filePointer.write(struct.pack('i', dimension)) #dimension


    def WriteSolutionsFieldsBinary(self,meshObject,SolsAtVertices=None,SolsAtTriangles=None,SolsAtTetrahedra=None):
        if SolsAtVertices is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject,62,SolsAtVertices)

        if SolsAtTriangles is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject,64,SolsAtTriangles)

        if SolsAtTetrahedra is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject,66,SolsAtTetrahedra)


    def _WriteSolutionsFieldsBinaryUsingKey(self,meshObject,key,Sols):

        NumberOfEntries = Sols[0].shape[0]

        self.filePointer.write(struct.pack('i', key))
        nbfields = len(Sols)
        #          1       +    1             1 +  nfields + numberofpoints*nfields*dataSize
        # endOfInformation + numberofpoints + 1 + nfields +  +
        self.PrintDebug(Sols[0].shape)
        self.PrintDebug(nbfields)
        self.PrintDebug(NumberOfEntries)
        self.PrintDebug([x.shape for x in Sols])
        self.PrintDebug(self.dataSize)
        self.PrintDebug(self.dataType)
        self.PrintDebug(self.filePointer.tell())

        nbcoomp = 0
        for sol in Sols:
            if len(sol.shape) == 1:
                nbcoomp += 1
            else:
                nbcoomp += sol.shape[-1]

        self.PrintDebug(nbcoomp)
        self.PrintDebug("**************************************")

        endOfInformation = self.filePointer.tell()+4*(1+1+1+nbfields)+nbcoomp*NumberOfEntries*self.dataSize
        self.filePointer.write(struct.pack('i', endOfInformation))# end of information


        self.filePointer.write(struct.pack('i', NumberOfEntries))# numberofpoints
        #print(numberofpoints)
        self.filePointer.write(struct.pack('i',nbfields))# numberofpoints
        #print(nbfields)

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

        self.PrintDebug("position at end " + str(self.filePointer.tell()))
        self.PrintDebug("calculated position at end " + str(endOfInformation))
        if ( not (self.filePointer.tell() == endOfInformation) ) : raise Exception("Error in the writing code, please debug me!!!")

    def CloseSolutionFileBinary(self):
        self.filePointer.write(struct.pack('i', 54)) #dimension
        #self.filePointer.write(struct.pack('i', 00)) #dimension
        #self._isOpen = False
        #self.filePointer.close()

    def WriteASCII(self,meshObject,SolsAtVertices=None, solutionOnOwnFile= False, nodalRefNumber = None,elemRefNumber=None):

        #self.filePointer.write("# This file has been writen by the python routine MmgWriter of the BasicTools package\n")
        #self.filePointer.write("# For any question about this routine, please contact SAFRAN TECH Pole M&S - CAM Team\n")
        meshObject.ComputeGlobalOffset()
        self.filePointer.write("MeshVersionFormatted 2 \n")
        #
        self.filePointer.write("Dimension " + str(meshObject.GetDimensionality()) + "\n\n")
        self.filePointer.write("Vertices\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} \n\n".format(numberofpoints) )
        #
        posn = meshObject.GetPosOfNodes()

        #print("------------------")
        #print(self.fileName)
        #print("------------------")
        #print(posn[0,:])
        for n in range(numberofpoints):
               #np.savetxt(self.filePointer, posn[np.newaxis,n,:] ,newline="" )
               posn[np.newaxis,n,:] .tofile(self.filePointer, sep=" ")
               if nodalRefNumber is None:
                   self.filePointer.write(" 0 \n")
               else:
                   self.filePointer.write(" {} \n".format(int(nodalRefNumber[n])) )

        self.filePointer.write("\n" )



#        celtags = meshObject.GetNamesOfCellTags()

        if meshObject.IsConstantRectilinear():
            import BasicTools.FE.ElementNames
            elements = [ BasicTools.FE.ElementNames.Hexaedron_8 ]
        else:
            elements = meshObject.elements
#        tagcounter = 2
#        #for tagname in celtags:
        globalOffset = 0
        for name, elementContainer in elements.items():
            elemtype = ASCIIName[name]

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

                    #np.savetxt(self.filePointer, connectivity,fmt="%d", newline=" ", delimiter=" " )
                    connectivity.tofile(self.filePointer, sep=" ")

                    if elemRefNumber is None :
                        self.filePointer.write(" 0\n")#
                    else:
                        self.filePointer.write(" " + str(elemRefNumber[globalOffset+i]) + "\n")#
                        #self.filePointer.write(struct.pack("i", elemRefNumber[globalOffset+i] ))#


            globalOffset += data.GetNumberOfElements()


            self.filePointer.write("\n")
#            #except KeyError:
#            #  continue
#          #tagcounter += 1

        if "Corners" in meshObject.nodesTags:
            tag = meshObject.nodesTags['Corners']
            if len(tag):
                self.filePointer.write("Corners\n");
                self.filePointer.write("{} \n\n".format(len(tag)) )
                #self.filePointer.write("{} \n".format(tag.GetIds()+1 ) )
                (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
                self.filePointer.write("\n" )

        for TagNameInFile,(ElementType,TagName) in ASCIITags.items():
            elements = meshObject.GetElementsOfType(ElementType)
            if TagName in elements.tags:
                tag = elements.tags[TagName ]
                if len(tag):
                    self.filePointer.write(str(TagNameInFile)+"\n");
                    self.filePointer.write("{} \n\n".format(len(tag)) )
                    #self.filePointer.write("{} \n".format(tag.GetIds()+1 ) )
                    (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
                    self.filePointer.write("\n" )

        self.filePointer.write("\n")

        if solutionOnOwnFile :
            #self.filePointer.write("End\n")
            self.Close();
            self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol" , 'w',0)
            self._isOpen = True
            self.filePointer.write("# This file has been writen by the python routine MmgWriter of the BasicTools package\n")
            self.filePointer.write("# For any question about this routine, please contact SAFRAN TECH Pole M&S Team OT\n")

            self.filePointer.write("MeshVersionFormatted\n2 \n")
            #
            self.filePointer.write("Dimension\n" + str(meshObject.GetDimensionality()) + "\n\n")


        if SolsAtVertices is not None:
             self._WriteSolutionsFieldsAsciiUsingKey(meshObject,"SolAtVertices",SolsAtVertices)

#            self.filePointer.write("SolAtVertices\n")
#            self.filePointer.write("{} \n".format(numberofpoints) )
#            self.filePointer.write("{} ".format(len(SolsAtVertices)) )
#            for sol in SolsAtVertices:
#                # we add a extra axis for scalar field stored in a vector (only one index)
#                if len(sol.shape)== 1:
#                    sol = sol[:,np.newaxis]
#                size = sol.shape[1]
#                if size == 1:
#                    self.filePointer.write("1 ")
#                elif size == meshObject.GetDimensionality():
#                    self.filePointer.write("2 ")
#                else:
#                    self.filePointer.write("3 ")
#            self.filePointer.write("\n ")
#
#            for i in xrange(numberofpoints):
#                for sol in SolsAtVertices:
#                    if len(sol.shape)== 1:
#                        sol = sol[:,np.newaxis]
#                    np.savetxt(self.filePointer, sol[i,:], newline=" ", delimiter=" " )
#                self.filePointer.write("\n ")


        #self.filePointer.write("End\n")

def CheckIntegrity():
    import BasicTools.FE.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir
    TestTempDir().SetTempPath("/tmp/BasicTools_Test_Directory_QUaZuy_safe_to_delete")

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5,6],dtype=np.int)

    mymesh.nodesTags.CreateTag("coucou").AddToTag(0)
    mymesh.nodesTags.CreateTag(Corners).AddToTag(0)

    tets = mymesh.GetElementsOfType(EN.Tetrahedron_4)
    tets.AddNewElement([0,1,2,4],0)

    tris = mymesh.GetElementsOfType(EN.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)

    tris.tags.CreateTag("RequiredTriangles").AddToTag(0)

    bars = mymesh.GetElementsOfType(EN.Bar_2)
    bars.AddNewElement([0,1],0)
    bars.tags.CreateTag("Ridges").AddToTag(0)


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
    sol = np.arange(mymesh.GetNumberOfNodes(),dtype=np.float)
   # sol.shape = (mymesh.GetNumberOfNodes(),1)
    WriteMesh(tempdir+"Test_MmgWriter_II_Binary.mesh", mymesh ,SolsAtVertices=[sol ] )
    WriteMesh(tempdir+"Test_MmgWriter_II_Ascii.mesh", mymesh ,SolsAtVertices=[sol ], binary=False )
    #print(tempdir)
    #print(open(tempdir+"Test_GmshWriter_II.geof").read())

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
