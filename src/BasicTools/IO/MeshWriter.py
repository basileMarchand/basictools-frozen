# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Mesh file writer
"""

from BasicTools.IO.IOFactory import RegisterWriterClass
import struct
import numpy as np

from BasicTools.Containers.MeshBase import Tag as Tag
from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.Containers.ElementNames as EN

from BasicTools.IO.MeshTools import BinaryKeywords, BinaryNumber, ASCIIName, ASCIITags, FieldTypes
from BasicTools.IO.MeshTools import Corners, Ridges, RequiredEdges, RequiredTriangles, RequiredVertices, GetTypesForVersion, BinaryName
import BasicTools.IO.MeshTools as MT
from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType



def WriteMesh(fileName, mesh, PointFields=None, solutionOnOwnFile=False, binary=True, nodalRefNumber=None, elemRefNumber=None):
    """Export Mesh to disk in the mesh format file.
        A file is created using the path and name of fileName

    Parameters
    ----------
    fileName : str
        name with path to the file to be created (relative or absolute)
    mesh : UnstructuredMesh
        the mesh to be written
    PointFields : list[np.ndarray], optional
        list of fields defined at the vertices of the mesh, to be included
        in the output file, by default None
    solutionOnOwnFile : bool, optional
        if True, the solution is written in a separate file, by default False
    binary : bool, optional
        if True, a binary file is produced, by default True
    nodalRefNumber : np.ndarray, optional
        a provided numbering for the vertices, by default None
    elemRefNumber : np.ndarray, optional
        a provided numbering for the elements, by default None
    """

    OW = MeshWriter()
    OW.SetBinary(binary)
    OW.Open(fileName)
    OW.Write(mesh, PointFields=PointFields, solutionOnOwnFile=solutionOnOwnFile,
            nodalRefNumber=nodalRefNumber, elemRefNumber=elemRefNumber)
    OW.Close()


class MeshWriter(WriterBase):
    """Class to write Unstructured mesh on disk in the mesh format file
    """
    def __init__(self):
        super(MeshWriter, self).__init__()

        self.SetVersion(2)

    def __str__(self):
        res = 'MeshWriter : \n'
        res += '   FileName : '+str(self.fileName) + '\n'
        res += '   Binary : ' + ('True' if self.isBinary() else 'False') + "\n"
        return res

    def SetFileName(self, fileName):
        """Sets the name of the output file

        Parameters
        ----------
        fileName : str
            name of the output file
        """
        self.fileName = fileName

    def SetSinglePrecision(self, single=True):
        """Sets the output precision to single precision

        Parameters
        ----------
        single : bool, optional
            if True, output is save in singe precision, by default True
        """
        if single :
            self.SetVersion(1)

    def SetVersion(self, version: int ):
        """Sets the output file format version

        Parameters
        ----------
        version : int
            format version
        """
        self.version = version
        posData, intData, floatData = GetTypesForVersion(version)

        self.posSize,   self.posFormat,   self.posType = posData
        self.intSize,   self.intFormat,   self.intType = intData
        self.floatSize, self.floatFormat, self.floatType = floatData


    def Write(self, meshObject, PointFields=None, solutionOnOwnFile=False, nodalRefNumber=None, elemRefNumber=None, PointFieldsNames=None, CellFieldsNames=None, CellFields=None):
        """Write mesh to file in mesh format

        Parameters
        ----------
        meshObject : UnstructuredMesh
            the mesh to be written
        PointFields : list[np.ndarray], optional
            list of fields defined at the vertices of the mesh, to be included
            in the output file, by default None
        solutionOnOwnFile : bool, optional
            if True, the solution is written in a separate file, by default False
        nodalRefNumber : np.ndarray, optional
            a provided numbering for the vertices, by default None
        elemRefNumber : np.ndarray, optional
            a provided numbering for the elements, by default None
        PointFieldsNames : None
            Not Used, by default None
        CellFieldsNames : None
            Not Used, by default None
        CellFields : None
            Not Used, by default None
        """
        if self.isBinary():
            return self.WriteBINARY(meshObject, PointFields=PointFields, solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber, elemRefNumber=elemRefNumber)
        else:
            return self.WriteASCII(meshObject, PointFields=PointFields, solutionOnOwnFile=solutionOnOwnFile, nodalRefNumber=nodalRefNumber, elemRefNumber=elemRefNumber)

    def GetDimensionFromMesh(self, meshObject):
        """Returns dimension of the mesh

        Parameters
        ----------
        meshObject : UnstructuredMesh
            mesh for which the dimension is returned

        Returns
        -------
        int
            dimension of the mesh
        """
        flat = True
        if meshObject.nodes.shape[1] == 3:
            mmax = np.max(meshObject.nodes[:, 2])
            mmin = np.min(meshObject.nodes[:, 2])
            flat = mmax == mmin == 0

        if meshObject.nodes.shape[1] == 3 and not flat:
            dimension = 3
        else:
            dimension = 2
        return dimension

    def WriteKeyWord(self, keyword):
        self.filePointer.write(struct.pack('i', BinaryKeywords[keyword]))

    def WriteInt(self, value):
        self.filePointer.write(struct.pack(self.intFormat, value))

    def WritePos(self, value):
        self.filePointer.write(struct.pack(self.posFormat, value))

    def CheckPos(self,value):
        if value != self.filePointer.tell():
            print(f"current endOfInformation: {value}")
            print(f"tell() {self.filePointer.tell()}")
            raise Exception("Error in the writing code, please debug me!!!")

    def WriteIntArray(self, intArray):
        nbids = len(intArray)
        currentposition = self.filePointer.tell()
        endOfInformation = currentposition +self.posSize+ self.intSize+ (nbids)*self.intSize
        self.WritePos(endOfInformation)
        self.WriteInt(nbids)
        (intArray).astype(self.intType).tofile(self.filePointer, format=self.intFormat, sep='')
        self.CheckPos(endOfInformation)

    def WriteBINARY(self, meshObject, PointFields=None, solutionOnOwnFile=False, nodalRefNumber=None, elemRefNumber=None):
        """Write mesh to file in binary mesh format

        Parameters
        ----------
        meshObject : UnstructuredMesh
            the mesh to be written
        PointFields : list[np.ndarray], optional
            list of fields defined at the vertices of the mesh, to be included
            in the output file, by default None
        solutionOnOwnFile : bool, optional
            if True, the solution is written in a separate file, by default False
        nodalRefNumber : np.ndarray, optional
            a provided numbering for the vertices, by default None
        elemRefNumber : np.ndarray, optional
            a provided numbering for the elements, by default None
        """
        # key MeshVersionFormatted  (always int32)
        self.WriteKeyWord("GmfVersionFormatted")
        self.filePointer.write(struct.pack('i', self.version))

        # key Dimension (3)  (always int32)
        dimension = self.GetDimensionFromMesh(meshObject)
        self.WriteKeyWord("GmfDimension")
        endOfInformation = self.filePointer.tell()+self.posSize+4
        self.WritePos(endOfInformation)
        self.filePointer.write(struct.pack('i', dimension))

        # key Vertices (4)
        #self.filePointer.write(struct.pack('i', BKeys["GmfVertices"]))
        self.WriteKeyWord("GmfVertices")

        currentposition = self.filePointer.tell()
        numberofpoints = meshObject.GetNumberOfNodes()
        endOfInformation = currentposition + self.posSize + self.intSize + numberofpoints*(self.floatSize*dimension+self.intSize)
        #self.filePointer.write(struct.pack('i', endOfInformation))  # end of information
        self.WritePos(endOfInformation)

        #self.filePointer.write(struct.pack('i', numberofpoints))  # numberofpoints
        self.WriteInt(numberofpoints)

        posn = meshObject.GetPosOfNodes()[:, 0:dimension].astype(self.floatType)

        if nodalRefNumber is None:
            nrn = np.ones((numberofpoints, 1), dtype=self.intType).squeeze()
        else:
            nrn = nodalRefNumber.astype(dtype=self.intType).squeeze()

        names = ["c"+str(i) for i in range(dimension)]+["id"]
        dtype_out = np.dtype({"names": names, "formats": [self.floatType]*dimension + [self.intType]})
        data = np.empty((posn.shape[0]), dtype=dtype_out)
        for i in range(dimension):
            data[names[i]] = posn[:, i]
        data["id"] = nrn
        self.filePointer.write(data.tobytes())

        self.CheckPos(endOfInformation)

        if MT.RequiredVertices in meshObject.nodesTags:
            ids = meshObject.nodesTags[MT.RequiredVertices].GetIds()+1
            nbids = len(ids)
            if nbids:
                #self.filePointer.write(struct.pack('i', BinaryKeywords["GmfRequiredVertices"]))
                self.WriteKeyWord("GmfRequiredVertices")

                currentposition = self.filePointer.tell()
                endOfInformation = currentposition + self.posSize+ self.intSize+len(ids)*self.intSize
                #self.filePointer.write(struct.pack('i', endOfInformation))  # end of information
                self.WritePos(endOfInformation)

                #self.filePointer.write(struct.pack('i', nbids))  # GetNumberOfElements
                self.WriteInt(nbids)
                ids.astype(self.intType).tofile(self.filePointer, format=self.intFormat, sep='')
                self.CheckPos(endOfInformation)

        globalOffset = 0
        for elementContainer in meshObject.elements:
            self.PrintDebug("Output of " + str(elementContainer))
            data = meshObject.elements[elementContainer]
            nbelements = data.GetNumberOfElements()
            if nbelements == 0:
                self.PrintVerbose("Empty element container (skipping) + " + elementContainer)
                continue

            if elementContainer == EN.Point_1:
                print("MeshWriter warning: ignoring EN.Point_1 elements ")
                globalOffset += data.GetNumberOfElements()
                continue

            #self.filePointer.write(struct.pack('i', elemtype))
            self.WriteKeyWord(BinaryName[elementContainer])

            nbNodes = data.GetNumberOfNodesPerElement()

            currentposition = self.filePointer.tell()

            endOfInformation = currentposition + self.posSize+self.intSize+(data.GetNumberOfElements()*(nbNodes+1))*self.intSize
            self.WritePos(endOfInformation)
            self.WriteInt(nbelements)

            tempcoon = np.zeros((data.GetNumberOfElements(), nbNodes+1), dtype=self.intType, order="c")
            tempcoon[:, 0:nbNodes] = data.connectivity.astype(self.intType)
            tempcoon[:, 0:nbNodes] += 1

            if elemRefNumber is not None:
                tempcoon[:, nbNodes] = elemRefNumber[globalOffset: globalOffset+nbelements]
            tempcoon.tofile(self.filePointer, format=self.intFormat, sep='')

            globalOffset += data.GetNumberOfElements()
            self.CheckPos(endOfInformation)

        if "Corners" in meshObject.nodesTags:
            tag = meshObject.nodesTags['Corners']
            ids = tag.GetIds()
            nbids = len(ids)
            if nbids:
                self.WriteKeyWord("GmfCorners")
                self.WriteIntArray(ids+1)

                #currentposition = self.filePointer.tell()
                #endOfInformation = currentposition + self.posSize+ self.intSize+(nbids)*self.intSize
                #self.WritePos(endOfInformation)
                #self.WriteInt(nbids)
                #(ids+1).astype(np.int32).tofile(self.filePointer, format=self.floatFormat, sep='')
                #self.CheckPos(endOfInformation)


        bars = meshObject.GetElementsOfType(EN.Bar_2)
        if "Ridges" in bars.tags and len(bars.tags["Ridges"]):

            self.PrintDebug("output Ridges")

            tag = bars.tags["Ridges"]
            ids =  tag.GetIds().astype(self.intType).ravel()
            nbids = len(ids)
            if nbids:
                self.WriteKeyWord("GmfRidges")
                self.WriteIntArray(ids+1)
                #currentposition = self.filePointer.tell()
                #endOfInformation = currentposition +self.posSize+ self.intSize+ (nbids)*self.intSize
                #self.WritePos(endOfInformation)
                #self.WriteInt(nbids)
                #(ids+1).astype(np.int32).tofile(self.filePointer, format=self.intFormat, sep='')
                #self.CheckPos(endOfInformation)

        bars = meshObject.GetElementsOfType(EN.Bar_2)
        if RequiredEdges in bars.tags and len(bars.tags[RequiredEdges]):

            self.PrintDebug("output RequiredEdges")
            ids = np.empty(0, dtype=int)
            if RequiredEdges in bars.tags:
                tag = bars.tags[RequiredEdges]

                ids = np.concatenate((ids, tag.GetIds()))
            nbids = len(ids)
            if nbids:
                self.WriteKeyWord("GmfRequiredEdges")
                self.WriteIntArray(ids+1)
                #currentposition = self.filePointer.tell()
                #endOfInformation = currentposition + self.posSize+ self.intSize+ (nbids)*self.intSize
                #self.WritePos(endOfInformation)
                #self.WriteInt(nbids)
                #(ids+1).astype(self.intType).tofile(self.filePointer, format=self.intFormat, sep='')
                #self.CheckPos(endOfInformation)

        tris = meshObject.GetElementsOfType(EN.Triangle_3)
        if "RequiredTriangles" in tris.tags and len(tris.tags["RequiredTriangles"]):

            self.PrintDebug("output RequiredTriangles")
            ids = np.empty(0, dtype=int)
            if "RequiredTriangles" in tris.tags:
                tag = tris.tags["RequiredTriangles"]

                ids = np.concatenate((ids, tag.GetIds()))
            nbids = len(ids)
            if nbids:
                self.WriteKeyWord("GmfRequiredTriangles")
                self.WriteIntArray(ids+1)

        if PointFields is not None and len(PointFields) > 0:
            if solutionOnOwnFile:
                self.Close()
                self.OpenSolutionFileBinary(support=meshObject)
            self.WriteSolutionsFieldsBinary(meshObject, PointFields)

    def OpenSolutionFileAscii(self, support):
        self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol", 'w')
        self._isOpen = True

        self.filePointer.write("MeshVersionFormatted 2\n\n")

        dimension = self.GetDimensionFromMesh(support)
        self.filePointer.write("Dimension {}\n\n".format(dimension))

    def WriteSolutionsFieldsAscii(self, meshObject, PointFields=None, SolsAtTriangles=None, SolsAtTetrahedra=None):
        if PointFields is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject, "SolAtVertices", PointFields)

        if SolsAtTriangles is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject, "SolAtTriangles", SolsAtTriangles)

        if SolsAtTetrahedra is not None:
            self._WriteSolutionsFieldsAsciiUsingKey(meshObject, "SolAtTetrahedra", SolsAtTetrahedra)

    def _WriteSolutionsFieldsAsciiUsingKey(self, meshObject, keyword, Sols):
        if len(Sols) == 0:
            return
        nbentries = Sols[0].shape[0]

        # key SolAtVertices = 62
        self.filePointer.write("{}\n".format(keyword))

        self.filePointer.write("{}\n".format(nbentries))

        nbfields = len(Sols)

        self.filePointer.write("{} ".format(nbfields))

        for sol in Sols:
            if len(sol.shape) == 1:
                sol = sol[:, np.newaxis]
            size = sol.shape[1]

            if size == 1:
                self.filePointer.write("1 ")
            elif size == meshObject.GetDimensionality():
                self.filePointer.write("2 ")
            else:
                self.filePointer.write("3 ")

        self.filePointer.write("\n\n")
        composedData = np.column_stack(Sols)
        np.savetxt(self.filePointer, composedData, fmt="%g")

        self.filePointer.write("\n")
        self.PrintDebug("tata")

    def CloseSolutionFileAscii(self):
        self.filePointer.write("End\n")  # dimension

    def OpenSolutionFile(self, support):
        if self.isBinary():
            self.OpenSolutionFileBinary(support)
        else:
            self.OpenSolutionFileAscii(support)

    def WriteSolutionsFields(self, meshObject, PointFields=None, SolsAtTriangles=None, SolsAtTetrahedra=None):
        if self.isBinary():
            self.WriteSolutionsFieldsBinary(meshObject, PointFields=PointFields,
                                            SolsAtTriangles=SolsAtTriangles, SolsAtTetrahedra=SolsAtTetrahedra)
        else:
            self.WriteSolutionsFieldsAscii(meshObject, PointFields=PointFields,
                                        SolsAtTriangles=SolsAtTriangles, SolsAtTetrahedra=SolsAtTetrahedra)

    def Close(self):
        if self.isBinary():
            self.CloseSolutionFileBinary()
        else:
            self.CloseSolutionFileAscii()
        WriterBase.Close(self)

    def OpenSolutionFileBinary(self, support):

        self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".solb", 'wb', 0)
        self._isOpen = True

        # key MeshVersionFormatted
        self.WriteKeyWord("GmfVersionFormatted")
        self.filePointer.write(struct.pack('i', self.version))
        #
        # key Dimension (3)
        self.WriteKeyWord("GmfDimension")
        self.filePointer.write(struct.pack('i', self.filePointer.tell()+4*2))  # end of information

        dimension = self.GetDimensionFromMesh(support)
        self.filePointer.write(struct.pack('i', dimension))  # dimension

    def WriteSolutionsFieldsBinary(self, meshObject, PointFields=None, SolsAtTriangles=None, SolsAtTetrahedra=None):
        if PointFields is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject, "GmfSolAtVertices", PointFields)

        if SolsAtTriangles is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject, "GmfSolAtTriangles", SolsAtTriangles)

        if SolsAtTetrahedra is not None:
            self._WriteSolutionsFieldsBinaryUsingKey(meshObject, "GmfSolAtTetrahedra", SolsAtTetrahedra)

    def _WriteSolutionsFieldsBinaryUsingKey(self, meshObject, key, Sols):

        if len(Sols) == 0:
            return

        NumberOfEntries = Sols[0].shape[0]

        self.WriteKeyWord(key)
        nbfields = len(Sols)

        nbcoomp = 0
        for sol in Sols:
            if len(sol.shape) == 1:
                nbcoomp += 1
            else :
                nbcoomp += sol.shape[-1]

        endOfInformation = self.filePointer.tell()+self.posSize + self.intSize*(1+1+nbfields)+nbcoomp*NumberOfEntries*self.floatSize
        self.WritePos(endOfInformation)

        self.WriteInt(NumberOfEntries)# numberofpoints
        self.WriteInt(nbfields)


        from BasicTools.IO.MeshTools import FieldTypes
        for sol in Sols:
            # we add a extra axis for scalar field stored in a vector (only one index)
            if len(sol.shape) == 1:
                sol = sol[:, np.newaxis]

            size = sol.shape[1]
            dim = meshObject.GetPointsDimensionality()
            if dim == 3 and np.all(meshObject.nodes[:, 2] == 0):
                dim = 2

            if size == 1:
                self.WriteInt( FieldTypes["GmfSca"])
            elif size == dim:
                self.WriteInt( FieldTypes["GmfVec"])
            elif size == dim*(dim+1)//2:
                self.WriteInt( FieldTypes["GmfSymMat"])
            elif size == dim*dim:
                self.WriteInt( FieldTypes["GmfMat"])
            else:
                raise Exception("Solution fields must be scalars, vectors , sym tensor or tensors.")

        np.column_stack(tuple(x.astype(self.floatType) for x in Sols)).tofile(self.filePointer, sep='')

        self.CheckPos(endOfInformation)


    def CloseSolutionFileBinary(self):
        self.WriteKeyWord("GmfEnd")

    def WriteASCII(self, meshObject, PointFields=None, solutionOnOwnFile=False, nodalRefNumber=None, elemRefNumber=None):
        """Write mesh to file in binary mesh format

        Parameters
        ----------
        meshObject : UnstructuredMesh
            the mesh to be written
        PointFields : list[np.ndarray], optional
            list of fields defined at the vertices of the mesh, to be included
            in the output file, by default None
        solutionOnOwnFile : bool, optional
            if True, the solution is written in a separate file, by default False
        nodalRefNumber : np.ndarray, optional
            a provided numbering for the vertices, by default None
        elemRefNumber : np.ndarray, optional
            a provided numbering for the elements, by default None
        """
        meshObject.ComputeGlobalOffset()
        self.filePointer.write("MeshVersionFormatted 2 \n")

        self.filePointer.write("Dimension " + str(meshObject.GetDimensionality()) + "\n\n")
        self.filePointer.write("Vertices\n")
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{}\n".format(numberofpoints))
        posn = meshObject.GetPosOfNodes()

        if nodalRefNumber is None:
            composedData = np.column_stack((posn, np.zeros((numberofpoints, 1), dtype=int)))
        else:
            composedData = np.column_stack((posn, nodalRefNumber))
        np.savetxt(self.filePointer, composedData, fmt="%g "*posn.shape[1] + "%i")

        self.filePointer.write("\n")

        if meshObject.IsConstantRectilinear():
            import BasicTools.Containers.ElementNames
            elements = [BasicTools.Containers.ElementNames.Hexaedron_8]
        else:
            elements = meshObject.elements

        globalOffset = 0
        for name, elementContainer in elements.items():
            if elementContainer.GetNumberOfElements() == 0:
                continue

            elemtype = ASCIIName.get(name, None)
            if elemtype is None:
                print("(MeshWriter) skiping this type of elements : " + name)
                continue

            if meshObject.IsConstantRectilinear():
                # hack to make the constantrectilinear api as the unstructured #TODO clean me
                class T():
                    pass
                data = T()
                data.connectivity = meshObject.GenerateFullConnectivity()
                nelem = meshObject.GetNumberOfElements()
            else:
                data = elementContainer
                nelem = data.GetNumberOfElements()
            if nelem == 0:
                continue
            self.filePointer.write("{}\n".format(elemtype))
            self.filePointer.write("{}\n".format(nelem))
            connectivity = (data.connectivity+1).astype(int)
            if elemRefNumber is None:
                composit_data = np.column_stack((connectivity, np.zeros(nelem, dtype=PBasicIndexType)))
            else:
                composit_data = np.column_stack(
                    (connectivity, elemRefNumber[globalOffset:globalOffset+data.GetNumberOfElements()]))
            composit_data.tofile(self.filePointer, sep=" ")

            globalOffset += data.GetNumberOfElements()
            self.filePointer.write("\n")

        nTags = [RequiredVertices, Corners]

        for tagname in nTags:
            if tagname in meshObject.nodesTags:
                tag = meshObject.nodesTags[tagname]
                if len(tag):
                    self.filePointer.write(tagname+" \n")
                    self.filePointer.write("{} \n\n".format(len(tag)))
                    (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
                    self.filePointer.write("\n")

        for TagNameInFile, (ElementType, TagName) in ASCIITags.items():
            elements = meshObject.GetElementsOfType(ElementType)
            if TagName in elements.tags:
                tag = elements.tags[TagName]
                if len(tag):
                    self.filePointer.write(str(TagNameInFile)+"\n")
                    self.filePointer.write("{} \n\n".format(len(tag)))
                    (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
                    self.filePointer.write("\n")

        self.filePointer.write("\n")

        if solutionOnOwnFile:
            self.Close()
            self.filePointer = open(".".join(self.fileName.split(".")[0:-1])+".sol", 'w')
            self._isOpen = True
            self.filePointer.write("#Written by BasicTools package\n")
            self.filePointer.write("MeshVersionFormatted\n2 \n")
            self.filePointer.write("Dimension\n" + str(meshObject.GetDimensionality()) + "\n\n")

        if PointFields is not None and len(PointFields) > 0:
            if solutionOnOwnFile:
                self.Close()
                self.OpenSolutionFileAscii(support=meshObject)
            self.WriteSolutionsFieldsAscii(meshObject, PointFields)

RegisterWriterClass(".mesh", MeshWriter)


def CreateMeshWriterBinary(ops):
    obj = MeshWriter()
    obj.SetBinary()
    return obj


RegisterWriterClass(".meshb", MeshWriter, CreateMeshWriterBinary)


def CheckIntegrity(GUI=False):
    import BasicTools.Containers.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    with mymesh.WithModification():
        mymesh.nodes = np.array([[0.00000000001, 0, 0], [1, 0, 0], [0, 1, 0],
                                [1, 1, 0], [0, 0, 1]], dtype=PBasicFloatType)
        mymesh.originalIDNodes = np.array([1, 3, 4, 5, 6], dtype=PBasicIndexType)

        #mymesh.nodesTags.CreateTag("coucou").AddToTag(0)
        mymesh.nodesTags.CreateTag(Corners).AddToTag(0)
        mymesh.nodesTags.CreateTag(MT.RequiredVertices).SetIds([0])

        tets = mymesh.GetElementsOfType(EN.Tetrahedron_4)
        tets.AddNewElement([0, 1, 2, 4], 0)

        tris = mymesh.GetElementsOfType(EN.Triangle_3)
        tris.AddNewElement([0, 1, 2], 0)
        tris.AddNewElement([2, 1, 3], 3)
        tris.originalIds = np.array([3, 5], dtype=PBasicIndexType)

        tris.tags.CreateTag("RequiredTriangles").AddToTag(0)

        #tris = mymesh.GetElementsOfType(EN.Bar_3)

        bars = mymesh.GetElementsOfType(EN.Bar_2)
        bars.AddNewElement([0, 1], 0)
        bars.tags.CreateTag("Ridges").AddToTag(0)
        bars.tags.CreateTag(MT.RequiredEdges).AddToTag(0)

    nodalRefNumber = np.arange(mymesh.GetNumberOfNodes())
    elemRefNumber = np.arange(mymesh.GetNumberOfElements())

    WriteMesh(tempdir+"CheckIntegrity_with_refs.mesh", mymesh,
            PointFields=[nodalRefNumber*10], elemRefNumber=elemRefNumber, nodalRefNumber=nodalRefNumber, solutionOnOwnFile=True)
    WriteMesh(tempdir+"CheckIntegrity_with_refs.mesh", mymesh, binary=False,
            PointFields=[nodalRefNumber*10], elemRefNumber=elemRefNumber, nodalRefNumber=nodalRefNumber, solutionOnOwnFile=True)

    OW = MeshWriter()
    OW.SetBinary(False)
    if GUI:
        OW.SetGlobalDebugMode(True)
    OW.Open(tempdir+"Test_MmgWriter_ASCII.mesh")
    print(OW)
    OW.Write(mymesh)
    OW.Close()
    OW.SetGlobalDebugMode()
    for v in [1,2,3,4]:
        print(f"-- {v} ------------------------------------------------------------------")
        OWB = MeshWriter()
        OWB.SetBinary(True)
        OWB.SetVersion(v)
        fname = tempdir+f"Test_MmgWriter_BIN_V{v}.meshb"
        OWB.Open(fname)
        OWB.Write(mymesh)
        OWB.Close()
        print(mymesh)
        from BasicTools.IO.MeshReader import ReadMesh
        mesh2 = ReadMesh(fname)
        print(f"-- {v} ------------------------------------------------------------------")

        print(mesh2)
        from BasicTools.Containers.MeshTools import IsClose
        if not IsClose(mymesh, mesh2):
            raise
    print("--------------------------------------------------------------------")

    print(mymesh)
    sol = np.arange(mymesh.GetNumberOfNodes(), dtype=PBasicFloatType)
    WriteMesh(tempdir+"Test_MmgWriter_II_Binary.mesh", mymesh, PointFields=[sol])
    WriteMesh(tempdir+"Test_MmgWriter_II_Ascii.mesh", mymesh, PointFields=[sol], binary=False)

    res = CreateMeshWriterBinary({})
    print(res)

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))  # pragma: no cover
