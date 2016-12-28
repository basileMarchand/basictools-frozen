# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import Tag as Tag
from OTTools.IO.WriterBase import WriterBase as WriterBase
import numpy as np
import OTTools.FE.ElementNames as EN
import struct

mmgName = {}
mmgName[EN.Bar_2]         = 'Edges'
mmgName[EN.Triangle_3]    = 'Triangles'
mmgName[EN.Quadrangle_4]  = 'Quadrilaterals'
mmgName[EN.Tetrahedron_4] = 'Tetrahedra'
mmgName[EN.Hexaedron_8]   = 'Hexahedra'
#mmgName[EN.Wedge_6]       = ''
#mmgName[EN.Pyramid_5]     = ''
#mmgName[EN.Point_1]   = ''


# for binary files
meshNumber = {}
meshNumber[EN.Bar_2] = 5
meshNumber[EN.Triangle_3] = 6
meshNumber[EN.Quadrangle_4] = 7
meshNumber[EN.Tetrahedron_4] = 8


def WriteMesh(filename,mesh,SolAtVertices=None,solutionOnOwnFile= False, binary=True, nodalRefNumber= None,elemRefNumber=None):
    OW = MeshWriter()
    OW.SetBinary(binary);
    OW.Open(filename)
    OW.Write(mesh,SolAtVertices = SolAtVertices,solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber,elemRefNumber=elemRefNumber)
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

    def Write(self,meshObject,SolAtVertices=None, solutionOnOwnFile= False, nodalRefNumber= None, elemRefNumber=None):
        if self.isBinary():
            return self.WriteBINARY(meshObject,SolAtVertices=SolAtVertices, solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber,elemRefNumber=elemRefNumber )
        else:
            return self.WriteASCII(meshObject,SolAtVertices=SolAtVertices, solutionOnOwnFile= solutionOnOwnFile, nodalRefNumber =nodalRefNumber,elemRefNumber=elemRefNumber )

    def WriteBINARY(self,meshObject,SolAtVertices=None, solutionOnOwnFile= False, nodalRefNumber = None,elemRefNumber=None):


        #key MeshVersionFormatted
        self.filePointer.write(struct.pack('i', 1))
        self.filePointer.write(struct.pack('i', self.dataSize/4))

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
        for n in xrange(numberofpoints):
             #self.filePointer.write(struct.pack(dataformat,  posn[n,:]))# end of information
             posn[n,:].astype(self.dataType).tofile(self.filePointer, sep='')

             if nodalRefNumber is  None :
                 self.filePointer.write(struct.pack("i", 0))# refs
             else :
                 self.filePointer.write(struct.pack("i", nodalRefNumber[n]))# refs



        print("position at end " + str(self.filePointer.tell()))
        print("calculate position at end " + str(endOfInformation))

        globalOffset =0
        for elementContainer in meshObject.elements.keys():
            print("Output of " + str(elementContainer ))
            elemtype = meshNumber[elementContainer]
            data = meshObject.elements[elementContainer]
            nbelements = data.GetNumberOfElements()
            if nbelements == 0:
                print("Empty elemnt container (is this corect???)")
                raise


            self.filePointer.write(struct.pack('i', elemtype ))
            nbNodes = data.GetNumberOfNodesPerElement()

            currentposition = self.filePointer.tell()
            endOfInformation = currentposition+ (2+data.GetNumberOfElements()*(nbNodes+1))*4
            self.filePointer.write(struct.pack('i', endOfInformation ))# end of information
            self.filePointer.write(struct.pack('i', nbelements))# GetNumberOfElements

            dataformat = "i"*nbNodes
            for i in xrange(nbelements):
                (data.connectivity[i,:]+1).astype(np.int32).tofile(self.filePointer, format=dataformat,sep='')
                if elemRefNumber is None :
                    self.filePointer.write(struct.pack("i", 0 ))#
                else:
                    self.filePointer.write(struct.pack("i", elemRefNumber[globalOffset+i] ))#

            globalOffset += data.GetNumberOfElements()

            print("position at end " + str(self.filePointer.tell()))
            print("calculate position at end " + str(endOfInformation))



        if SolAtVertices is not None:
            if solutionOnOwnFile :
                self.filePointer.write(struct.pack('i', 54)) #dimension
                self.Close();
                self._OpenSolutionFileBinary()
            self._WriteFieldBinary(SolAtVertices)

        #key End
        self.filePointer.write(struct.pack('i', 54)) #dimension
        self.Close()


    def OpenSolutionFileBinary(self,support):



        self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".solb" , 'wb',0)
        self.__isOpen = True

        #key MeshVersionFormatted
        self.filePointer.write(struct.pack('i', 1))
        self.filePointer.write(struct.pack('i', self.dataSize/4))
        #
        #key Dimension (3)
        self.filePointer.write(struct.pack('i', 3))
        dimension = support.GetDimensionality()
        self.filePointer.write(struct.pack('i', self.filePointer.tell()+4*2))# end of information
        self.filePointer.write(struct.pack('i', dimension)) #dimension

    def WriteSolutionsFieldsBinary(self,meshObject,SolAtVertices):
        numberofpoints = meshObject.GetNumberOfNodes()

        #key SolAtVertices
        self.filePointer.write(struct.pack('i', 62))
        nbfields = len(SolAtVertices)
        #          1       +    1             1 +  nfields + numberofpoints*nfields*dataSize
        # endOfInformation + numberofpoints + 1 + nfields +  +
        endOfInformation = self.filePointer.tell()+4*(1+1+1+nbfields)+nbfields*numberofpoints*self.dataSize
        self.filePointer.write(struct.pack('i', endOfInformation))# end of information


        self.filePointer.write(struct.pack('i', numberofpoints))# numberofpoints
        print(numberofpoints)
        self.filePointer.write(struct.pack('i',nbfields))# numberofpoints
        print(nbfields)

        for sol in SolAtVertices:
            self.filePointer.write(struct.pack('i',1))#

            # we add a extra axis for scalar field stored in a vector (only one index)
            if len(sol.shape)== 1:
                sol = sol[:,np.newaxis]
            size = sol.shape[1]
            if size > 1:
                raise Exception("for now only scalars field")

        for i in xrange(numberofpoints):
            for sol in SolAtVertices:
                if len(sol.shape)== 1:
                    sol = sol[:,np.newaxis]
                sol[i,:].astype(self.dataType).tofile(self.filePointer, sep='')

        self.PrintDebug("position at end " + str(self.filePointer.tell()))
        self.PrintDebug("calculate position at end " + str(endOfInformation))

    def CloseSolutionFileBinary(self):
        self.filePointer.write(struct.pack('i', 54)) #dimension
        self.Close()

    def WriteASCII(self,meshObject,SolAtVertices=None, solutionOnOwnFile= False, nodalRefNumber = None,elemRefNumber=None):
        if  nodalRefNumber is not None or elemRefNumber is not None:
            raise Exception('Not implemented yet')

        self.filePointer.write("# This file has been writen by the python routine MmgWriter of the OTTools package\n")
        self.filePointer.write("# For any question about this routine, please contact SAFRAN TECH Pole M&S Team OT\n")

        self.filePointer.write("MeshVersionFormatted \n2 \n")
        #
        self.filePointer.write("Dimension\n" + str(meshObject.GetDimensionality()) + "\n\n")
        self.filePointer.write("Vertices\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} \n".format(numberofpoints) )
        #
        posn = meshObject.GetPosOfNodes()

        print("------------------")
        print(self.fileName)
        print(posn[0,:])
        print("------------------")
        for n in xrange(numberofpoints):
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] ,newline="" )

               self.filePointer.write(" {} \n".format(n+1) )

        self.filePointer.write("\n" )



#        celtags = meshObject.GetNamesOfCellTags()

        if meshObject.IsConstantRectilinear():
            import OTTools.FE.ElementNames
            elements = [ OTTools.FE.ElementNames.Hexaedron_8 ]
        else:
            elements = meshObject.elements.keys()
#        tagcounter = 2
#        #for tagname in celtags:
        elementcpt = 1
        for elementContainer in elements:
            elemtype = mmgName[elementContainer]
            self.filePointer.write("{} \n".format(elemtype) )

            if meshObject.IsConstantRectilinear():
                #hack to make the constantrectilinear api as the unstructured #TODO clean me
                class T(): pass
                data = T()
                data.connectivity = meshObject.GenerateFullConnectivity()
                nelem = meshObject.GetNumberOfElements()
            else:
                data = meshObject.elements[elementContainer]
                nelem = data.GetNumberOfElements()

            self.filePointer.write("{} \n".format(nelem) )
            for i in xrange(nelem):
                    connectivity = (data.connectivity[i,:]+1).astype(int)

                    np.savetxt(self.filePointer, connectivity,fmt="%d", newline=" ", delimiter=" " )
                    self.filePointer.write(" {} \n".format(elementcpt) )
                    elementcpt += 1

            self.filePointer.write("\n")
#            #except KeyError:
#            #  continue
#          #tagcounter += 1

        if solutionOnOwnFile :
            self.filePointer.write("End\n")
            self.Close();
            self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol" , 'w',0)
            self.__isOpen = True
            self.filePointer.write("# This file has been writen by the python routine MmgWriter of the OTTools package\n")
            self.filePointer.write("# For any question about this routine, please contact SAFRAN TECH Pole M&S Team OT\n")

            self.filePointer.write("MeshVersionFormatted\n2 \n")
            #
            self.filePointer.write("Dimension\n" + str(meshObject.GetDimensionality()) + "\n\n")


        if SolAtVertices is not None:
            self.filePointer.write("SolAtVertices\n")
            self.filePointer.write("{} \n".format(numberofpoints) )
            self.filePointer.write("{} ".format(len(SolAtVertices)) )
            for sol in SolAtVertices:
                # we add a extra axis for scalar field stored in a vector (only one index)
                if len(sol.shape)== 1:
                    sol = sol[:,np.newaxis]
                size = sol.shape[1]
                if size == 1:
                    self.filePointer.write("1 ")
                elif size == meshObject.GetDimensionality():
                    self.filePointer.write("2 ")
                else:
                    self.filePointer.write("3 ")
            self.filePointer.write("\n ")

            for i in xrange(numberofpoints):
                for sol in SolAtVertices:
                    if len(sol.shape)== 1:
                        sol = sol[:,np.newaxis]
                    np.savetxt(self.filePointer, sol[i,:], newline=" ", delimiter=" " )
                self.filePointer.write("\n ")


        self.filePointer.write("End\n")

def CheckIntegrity():
    import OTTools.FE.UnstructuredMesh as UM

    from OTTools.Helpers.Tests import TestTempDir

    TestTempDir.SetTempPath("/home/fbordeu-weld/tmp/mmg/")
    tempdir = TestTempDir.GetTempPath()


    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=np.int)

    mymesh.nodesTags.CreateTag("coucou").AddToTag(0)

    tris = mymesh.GetElementsOfType('tri3')
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)


    OW = MeshWriter()
    print(OW)
    OW.Open(tempdir+"Test_MmgWriter.mesh")
    OW.Write(mymesh)
    OW.Close()

    OW.SetBinary(True)
    OW.Open(tempdir+"Test_MmgWriter.mesh")
    OW.Write(mymesh)
    OW.Close()



    print(mymesh)
    WriteMesh(tempdir+"Test_MmgWriter_II.mesh", mymesh  )
    #print(tempdir)
    #print(open(tempdir+"Test_GmshWriter_II.geof").read())

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
