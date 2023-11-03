# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Mesh file reader
"""

import struct
import numpy as np

from BasicTools.Containers.UnstructuredMesh import  UnstructuredMesh
import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.IO.IOFactory import RegisterReaderClass

from BasicTools.IO.MeshTools import ASCIITypes, ASCIITags, RequiredVertices, FieldTypes
from BasicTools.IO.MeshTools import BinaryTypes, BinaryTags, BinaryFields
from BasicTools.IO.MeshTools import BinaryKeywords as BKeys, GetTypesForVersion
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType

def ReadMesh(fileName:str=None, string:str=None, ReadRefsAsField:bool=False)-> UnstructuredMesh:
    """Function API for reading a Mesh file

    Parameters
    ----------
    fileName : str, optional
        file to read, by default None
    string : str, optional
        string to convert, by default None
    ReadRefsAsField : bool, optional
        if true the two fiels will be preset (at the nodes and element) with the nbame refs
        if false the values of the refs will be converted tags, using the prefix NTag or ETag plus the number of the ref,
        by default False

    Returns
    -------
    UnstructuredMesh
        The mesh with the field present in the file
    """
    reader = MeshReader()
    reader.SetReadRefsAsField(ReadRefsAsField)
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    return reader.Read()


def ReadSol(fileName, out=None):
    """Function API for reading a solution file associated to a Mesh file

    Parameters
    ----------
    fileName : str, optional
        file to read, by default None
    out : UnstructuredMesh, optional
        output unstructured mesh object containing reading result, by default None

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    reader = MeshSolutionReaderWrapper()
    reader.SetFileName(fileName)
    return reader.Read(out=out)

class MeshReader(ReaderBase):
    """Mesh Reader class
    """
    def __init__(self):
        super(MeshReader, self).__init__()
        self.refsAsAField = False
        self.dim = 3
        self.version = -1

    def SetReadRefsAsField(self, val):
        self.refsAsAField = bool(val)

    def Read(self, out=None):
        """Function that performs the reading of a Mesh file

        Parameters
        ----------
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        if self.binary:
            return self.ReadMeshBinary(out=out)
        else:
            return self.ReadMeshAscii(out=out)

    def SetFileName(self, fileName):
        """Sets the name of the file to read

        Parameters
        ----------
        fileName : str
            name of the file to read
        """
        super(MeshReader, self).SetFileName(fileName)

        if fileName is not None and fileName[-1] == "b":
            self.SetBinary(True)
        else:
            self.SetBinary(False)

    def ReadExtraFields(self, fileName):
        self.SetFileName(fileName)

        mesh = self.Read()

        res = {}
        res.update(mesh.nodeFields)
        res.update(mesh.elemFields)
        return res

##ASCII PART #################################################################
    def ReadMeshAscii(self, _fileName=None, _string=None, fieldFileName=None, out=None):
        """Function that performs the reading of an ascii Mesh file

        Parameters
        ----------
        _fileName : str, optional
            name of the file to be read, by default None
        _string : str, optional
            data to be read as a string instead of a file, by default None
        fieldFileName : str, optional
            Not Used, by default None
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        if _fileName is not None:
            self.SetFileName(_fileName)

        if _string is not None:
            self.SetStringToRead(_string)

        self.readFormat = 'r'
        self.StartReading()

        if out is None:
            res = UnstructuredMesh()
        else:
            res = out


        dataType = float
        refs_per_elementType = {}
        while (True):
            line = self.filePointer.readline()
            if line == "":
                break
            l = line.strip('\n').lstrip().rstrip()
            if len(l) == 0:
                continue

            if l.find("MeshVersionFormatted") > -1:
                if len(l.lstrip().rstrip().split()) > 1:
                    formatFile = (l.lstrip().rstrip().split()[1])
                else:
                    formatFile = int(self.filePointer.readline())
                if formatFile == 2:
                    dataType = np.float64
                continue

            if l.find("Dimension") > -1:
                s = l.split()
                if len(s) > 1:
                    s.pop(0)
                    l = " ".join(s)
                else:
                    l = self.filePointer.readline()
                dimension = int(l)
                self.PrintVerbose('Dimension : ' + str(dimension))
                continue

            if len(l) > 7 and l.find("Vertices") == 0:
                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()

                nbNodes = int(l.split()[0])
                self.PrintVerbose("Reading "+str(nbNodes) + " Nodes")
                res.nodes = np.empty((nbNodes, dimension), dtype=dataType)
                res.originalIDNodes = np.empty((nbNodes,), dtype=PBasicIndexType)
                cpt = 0
                refs = np.empty(nbNodes, dtype=PBasicIndexType)
                while (True):
                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0:
                        continue
                    s = l.split()

                    res.nodes[cpt, :] = list(map(float, s[0:dimension]))
                    #res.originalIDNodes[cpt] = cpt

                    #ref = s[dimension]
                    refs[cpt] = int(s[dimension])

                    #if not self.refsAsAField:
                    #    res.GetNodalTag("NTag"+ref).AddToTag(cpt)

                    cpt += 1
                    if cpt == nbNodes:
                        break

                res.originalIDNodes = np.arange(nbNodes, dtype=PBasicIndexType)

                if self.refsAsAField:
                    res.nodeFields['refs'] = refs
                else :
                    nb_refs = np.unique(refs)
                    for iref in nb_refs:
                        res.GetNodalTag("NTag"+str(iref)).AddToTag(np.where(refs==iref)[0])

                continue

            if l in ASCIITypes:
                elements = res.GetElementsOfType(ASCIITypes[l])
                nbNodes = elements.GetNumberOfNodesPerElement()
                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()

                nbElements = int(l.split()[0])
                if nbElements == 0:
                    continue
                self.PrintVerbose("Reading "+str(nbElements) + " Elements")
                cpt = 0
                refs = np.empty(nbElements, dtype=PBasicIndexType)
                while (True):
                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0:
                        continue

                    s = list(map(int, l.split()))
                    if (len(s)/(nbNodes+1)) != (len(s)//(nbNodes+1)):
                        print(len(s))
                        raise Exception("Incorrect Number of int in the file")
                    offset = 0
                    for elemNB in range(len(s)//(nbNodes+1)):

                        elements.AddNewElement(s[offset:nbNodes+offset], cpt)
                        ref = s[nbNodes+offset]
                        #elements.GetTag("ETag"+str(ref)).AddToTag(cpt)
                        offset += nbNodes + 1

                        refs[cpt] = int(ref)
                        cpt += 1
                    if nbElements == cpt:  # pragma: no cover
                        break
                refs_per_elementType[elements.elementType] = refs

                elements.connectivity -= 1
                continue

            if l in ASCIITags:
                elements = res.GetElementsOfType(ASCIITags[l][0])
                tag = elements.GetTag(ASCIITags[l][1])

                line = self.filePointer.readline()
                l = line.strip('\n').lstrip().rstrip()
                nbIds = int(l.split()[0])
                self.PrintVerbose("Reading tags Elements")
                cpt = 0
                ids = []
                while (True):
                    if nbIds == cpt:  # pragma: no cover
                        break

                    line = self.filePointer.readline()
                    l = line.strip('\n').lstrip().rstrip()

                    if len(l) == 0:
                        continue

                    newids = list(map(int, l.split()))
                    cpt += len(newids)

                    ids.extend(newids)

                tag.SetIds(np.array(ids, dtype=int)-1)
                continue

            if l in ["SolAtVertices"]:
                fieldname = l
                data = self._ReadFieldsASCII(self, dimension)
                for i in range(len(data)):
                    res.nodeFields[fieldname+str(i)] = data[i]
                continue

            if l in ["SolAtTetrahedra"]:
                fieldname = l
                data = self._ReadFieldsASCII(self, dimension)
                for i in range(len(data)):
                    res.elemFields[fieldname+str(i)] = data[i]
                continue

            if l in ["Tangents", "NormalAtVertices", "TangentAtVertices", "Normals"]:
                self.PrintVerbose("ignoring line :->" + l + "<-")
                line = self.filePointer.readline()
                nls = int(line.strip('\n').lstrip().rstrip())
                for i in range(nls):
                    line = self.filePointer.readline()
                continue

            if l.find("End") > -1:
                break

            self.PrintVerbose("ignoring line :->" + l + "<-")
        res.PrepareForOutput()
        if self.refsAsAField:
            refs = np.empty(res.GetNumberOfElements(), dtype=PBasicIndexType)
            cpt = 0
            for elementType, data in res.elements.items():
                refs[cpt:cpt+data.GetNumberOfElements()] = refs_per_elementType[elementType]
                cpt += data.GetNumberOfElements()
            res.elemFields['refs'] = refs
        else:
            for elementType, data in res.elements.items():
                refs = refs_per_elementType[elementType]
                nb_refs = np.unique(refs)
                for iref in nb_refs:
                    data.tags.CreateTag("ETag"+str(iref)).AddToTag(np.where(refs==iref)[0])

        self.EndReading()
        return res

    def _ReadFieldsASCII(self, myFile, dim):
        datares = []
        line = myFile.ReadCleanLine(withError=True)

        nbentries = int(line.split()[0])
        line = myFile.ReadCleanLine(withError=True)
        nbfields = int(line.split()[0])
        fieldSizes = [int(x) for x in line.split()[1:]]

        for i in range(len(fieldSizes)):
            if fieldSizes[i] == FieldTypes["GmfSca"]:
                fieldSizes[i] = 1
            elif fieldSizes[i] == FieldTypes["GmfVec"]:
                fieldSizes[i] = dim
            elif fieldSizes[i] == FieldTypes["GmfSymMat"]:
                fieldSizes[i] = dim*(dim+1)//2
            elif fieldSizes[i] == FieldTypes["GmfMat"]:
                fieldSizes[i] = dim*dim
            else:
                raise Exception(f"Field {i} not conform to format")

        ncoomp = np.sum(fieldSizes)

        data = np.fromfile(file=myFile.filePointer, dtype=float, count=int(ncoomp*nbentries), sep=" \n")

        data.shape = (nbentries, ncoomp)
        cpt = 0
        for i in range(nbfields):
            datares.append(data[:, cpt:cpt+fieldSizes[i]])
            cpt += fieldSizes[i]
        return datares

##BINARY PART ################################################################################
    def ReadBinaryInt(self):
        if self.intSize == 4:
            return self.readInt32()
        else:
            return self.readInt64()

    def ReadEndOfInformation(self):
        if self.posSize == 4:
            return self.readInt32()
        else:
            return self.readInt64()

    def ReadKeyWord(self):
        return self.readInt32()

    def ReadIntArray(self):
        nbEntries = self.ReadBinaryInt()
        ids = np.fromfile(self.filePointer, dtype=self.intType, count=nbEntries, sep="")
        return ids

    def ReadBinaryHeader(self):
        """Function that performs the reading of the header of a binary Mesh file

        Returns
        -------
        int
            dimension of the read mesh
        """
        key = self.ReadKeyWord()

        if key == BKeys["GmfVersionFormatted"]:
            self.version = self.readInt32()
            posData, intData, floatData = GetTypesForVersion(self.version)

            self.posSize, self.posFormat, self.posType = posData
            self.intSize, self.intFormat, self.intType = intData
            self.floatSize, self.floatFormat, self.floatType = floatData
        else:
            raise Exception('Expected key value equal to '+str(BKeys["GmfVersionFormatted"]) +
                            ' (GmfVersionFormatted), got : '+str(key)+' are you sure this file is binary??')


        key = self.ReadKeyWord()
        if key == BKeys["GmfDimension"]:
            endOfInformation = self.ReadEndOfInformation()
            dimension = self.readInt32()
            self.filePointer.seek(endOfInformation)
        else:
            raise Exception('Expected key value equal to ' +
                            str(BKeys["GmfDimension"])+' (GmfDimension), got : '+str(key)+' are you sure this file is binary??')

        return dimension

    def ReadMeshBinary(self, out=None):
        """Function that performs the reading of a binary Mesh file

        Parameters
        ----------
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        self.readFormat = 'rb'
        self.StartReading()

        f = self.filePointer

        if out is None:
            res = UnstructuredMesh()
        else:
            res = out

        dimension = self.ReadBinaryHeader()

        globalElementCounter = 0
        elemRefsDic = {}

        endOfInformation = self.filePointer.tell()
        while True:
            if endOfInformation != self.filePointer.tell():
                print(f"endOfInformation : {endOfInformation}")
                print(f"tell {self.filePointer.tell()}")
                print(f"key {key}")
                raise Exception("Error in endOfInformation")

            try:
                key = self.ReadKeyWord()
            except EOFError:
                break
            self.PrintVerbose("key" + str(key))
            if key == BKeys["GmfEnd"]:
                break

            endOfInformation = self.ReadEndOfInformation()
            # Vertices
            if key == BKeys["GmfVertices"]:
                nbNodes = self.ReadBinaryInt()

                self.PrintVerbose("Reading " + str(nbNodes) + " nodes ")

                res.nodes = np.empty((nbNodes, dimension), dtype=PBasicFloatType)
                res.originalIDNodes = np.empty((nbNodes,), dtype=PBasicIndexType)
                dt = np.dtype([('pos', self.floatType, (dimension,)), ('ref', self.intType, (1,))])

                data = np.fromfile(self.filePointer, dtype=dt, count=nbNodes, sep="")

                res.nodes[:, :] = data[:]["pos"].astype(PBasicFloatType)
                res.originalIDNodes[:] = np.arange(nbNodes, dtype=PBasicIndexType)

                refs = data[:]["ref"]
                if self.refsAsAField:
                    res.nodeFields['refs'] = refs
                else:
                    cpt = 0
                    while (cpt < len(refs)):
                        ref = refs[cpt][0]
                        res.GetNodalTag("NTag"+str(ref)).AddToTag(cpt)
                        cpt += 1
                continue

            if key == BKeys["GmfCorners"]:
                data =self.ReadIntArray()
                #nbCorners = self.ReadBinaryInt()
                #data = np.fromfile(self.filePointer, dtype=self.intType, count=nbCorners, sep="")
                res.nodesTags.CreateTag("Corners").SetIds(data-1)
                continue

            # all kind of elements
            if key in BinaryTypes:

                elements = res.GetElementsOfType(BinaryTypes[key])
                self.PrintVerbose("Reading elements of type " + elements.elementType)
                nbNodes = elements.GetNumberOfNodesPerElement()
                nbElements = self.ReadBinaryInt()
                self.PrintVerbose('Reading ' + str(nbElements) + " Elements ")
                if nbElements > 0:
                    elements.Reserve(nbElements)
                    elements.cpt = 0

                    dt = np.dtype([('conn', self.intType, (nbNodes,)), ('ref', self.intType, (1,))])

                    data = np.fromfile(self.filePointer, dtype=dt, count=nbElements, sep="")

                    elements.connectivity = (data[:]["conn"]-1).astype(PBasicIndexType)

                    elements.originalIds = np.arange(globalElementCounter, globalElementCounter+nbElements)
                    elements.cpt = nbElements
                    globalElementCounter += nbElements

                    refs = data[:]["ref"]
                    if self.refsAsAField:
                        elemRefsDic[elements.elementType] = refs
                    else:
                        urefs = np.unique(refs)
                        if len(urefs) > 1 or urefs[0] != 0:
                            for ref in urefs:
                                elements.GetTag("ETag"+str(ref)).SetIds(np.where(refs == ref)[0])
                continue

            if key == BKeys["GmfRequiredVertices"]:
                tagname = RequiredVertices
                self.PrintVerbose("Reading " + str(tagname))
                res.nodesTags.CreateTag(tagname).SetIds( self.ReadIntArray()-1)

                continue

            if key in BinaryTags:
                elemtype, tagname = BinaryTags[key]
                self.PrintVerbose("Reading " + str(tagname))
                elements = res.GetElementsOfType(elemtype)
                elements.tags.CreateTag(tagname).SetIds(self.ReadIntArray()-1)
                continue

            if key in BinaryFields:
                data = self._readExtraFieldBinary(dimension)
                for i in range(len(data)):
                    res.elemFields[BinaryFields[key]+str(i)] = data[i]
                continue

            if key == BKeys["GmfSolAtVertices"]:
                data = self._readExtraFieldBinary(dimension)
                for i in range(len(data)):
                    res.nodeFields["SolAtVertices"+str(i)] = data[i]
                continue

            if key not in [BKeys[x] for x in ["GmfNormals", "GmfNormalAtVertices", "GmfTangents", "GmfTangentAtVertices"]]:
                self.PrintVerbose("skiping key : " + str(key))
            f.seek(endOfInformation)

        res.GenerateManufacturedOriginalIDs()
        res.PrepareForOutput()

        if self.refsAsAField:
            elemRefs = np.empty(globalElementCounter, dtype=PBasicIndexType)
            cpt = 0
            for name, val in res.elements.items():
                elemRefs[cpt:cpt+val.GetNumberOfElements()] = elemRefsDic[name].ravel()
                cpt += val.GetNumberOfElements()
            res.elemFields['refs'] = elemRefs

        self.EndReading()
        res.ConvertDataForNativeTreatment()
        return res

    def _readExtraFieldBinary(self, dim):
        nbEntities = self.ReadBinaryInt()
        nbfields = self.ReadBinaryInt()

        fieldSizes = np.empty(nbfields, dtype=PBasicIndexType)
        res = []
        for i in range(nbfields):
            fieldType = self.ReadBinaryInt()
            if fieldType == FieldTypes["GmfSca"]:
                fieldSizes[i] = 1
            elif fieldType == FieldTypes["GmfVec"]:
                fieldSizes[i] = dim
            elif fieldType == FieldTypes["GmfSymMat"]:
                fieldSizes[i] = dim*(dim+1)//2
            elif fieldType == FieldTypes["GmfMat"]:
                fieldSizes[i] = dim*dim
            else:
                raise Exception(f"Field {i} not conform to format")

        ncoomp = np.sum(fieldSizes)

        data = np.fromfile(self.filePointer, dtype=self.floatType, count=nbEntities*ncoomp, sep="")
        data.shape = (nbEntities, ncoomp)

        cpt = 0
        for i in range(nbfields):
            res.append(data[:, cpt:cpt+fieldSizes[i]])
            cpt += fieldSizes[i]
        return res


class MeshSolutionReaderWrapper():
    """Class to handle a solution file associated to a Mesh file
    """
    def __init__(self):
        super(MeshSolutionReaderWrapper, self).__init__()
        import locale
        self.canHandleTemporal = False
        self.encoding = locale.getpreferredencoding(False)

    def SetFileName(self, fileName):
        """Sets the name of the file to read

        Parameters
        ----------
        fileName : str
            name of the file to read
        """
        import os.path
        self.fileName = fileName
        if fileName is None or fileName == "None":
            return
        dirname = os.path.dirname(fileName)
        basename, extention = os.path.splitext(os.path.basename(fileName))

        if extention[-1] == "b":
            f = os.path.join(dirname, basename+".meshb")
        else:
            f = os.path.join(dirname, basename+".mesh")

        # we check if the file exist, if not we try the other type
        if not os.path.isfile(f):
            if extention[-1] == "b":
                f = os.path.join(dirname, basename+".mesh")
            else:
                f = os.path.join(dirname, basename+".meshb")

        if not os.path.isfile(f):
            raise Exception("unable to find a mesh file")
        self.reader = MeshReader()
        self.reader.encoding = self.encoding
        self.reader.SetFileName(fileName=f)

    def Read(self, out=None):
        """Function that performs the reading of a solution file associated to a Mesh file

        Parameters
        ----------
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        if out:
            raise

        mesh = self.reader.Read()
        fields = self.reader.ReadExtraFields(self.fileName)
        mesh.nodeFields = {k: v for k, v in fields.items() if k.find("SolAtVertices") != -1}
        if 'SolAtTetrahedra0' in fields:
            if mesh.GetElementsOfType(EN.Tetrahedron_4).GetNumberOfElements() == mesh.GetNumberOfElements():
                mesh.elemFields = {k: v for k, v in fields.items() if k.find("SolAtTetrahedra") != -1}
        return mesh


RegisterReaderClass(".mesh", MeshReader)
RegisterReaderClass(".meshb", MeshReader)

RegisterReaderClass(".sol", MeshSolutionReaderWrapper)
RegisterReaderClass(".solb", MeshSolutionReaderWrapper)


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

    __teststringField = u"""
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

    __teststringFieldMatSym = u"""
MeshVersionFormatted 2

Dimension 3

SolAtVertices
4
1 3

1. 1. 1. 1. 1. 1.
2. 2. 2. 2. 2. 2.
3. 3. 3. 3. 3. 3.
4. 4. 4. 4. 4. 4.

End
"""
    res = ReadMesh(string=__teststring)
    print(res)

    from BasicTools.Helpers.Tests import TestTempDir

    newFileName = TestTempDir().GetTempPath()+"mshFile.mesh"
    open(newFileName, 'w').write(__teststring)
    res = ReadMesh(fileName=newFileName)

    newFileName = TestTempDir().GetTempPath()+"mshFile.sol"
    f = open(newFileName, 'w')
    f.write(__teststringField)
    f.close()

    print('Reading : ' + str(newFileName))

    resfield = MeshReader().ReadExtraFields(fileName=newFileName)
    resfield = ReadSol(fileName=newFileName)

    newFileName = TestTempDir().GetTempPath()+"mshFileMatSym.mesh"
    open(newFileName, 'w').write(__teststring)

    newFileName = TestTempDir().GetTempPath()+"mshFileMatSym.sol"
    f = open(newFileName, 'w')
    f.write(__teststringFieldMatSym)
    f.close()

    print('Reading : ' + str(newFileName))

    resfieldmatsym = MeshReader().ReadExtraFields(fileName=newFileName)
    resfieldmatsym = ReadSol(fileName=newFileName)

    from BasicTools.IO.MeshWriter import WriteMesh as WriteMesh

    print("###########################################################")
    print(res)
    newFileName = TestTempDir().GetTempPath()+"mshFile.meshb"
    WriteMesh(newFileName, res, binary=True)

    res = ReadMesh(newFileName, ReadRefsAsField=True)
    print(res)

    sol = MeshReader().ReadExtraFields(TestTempDir().GetTempPath()+"mshFile.sol")
    print(sol)

    newFileName = TestTempDir().GetTempPath()+"mshFileMatSym.meshb"
    WriteMesh(newFileName, res, [resfieldmatsym.nodeFields["SolAtVertices0"]], solutionOnOwnFile=True, binary=True)
    newFileName = TestTempDir().GetTempPath()+"mshFileMatSym.solb"
    resfieldmatsym = ReadSol(fileName=newFileName)

    from BasicTools.IO.MeshWriter import MeshWriter as MeshWriter
    mw = MeshWriter()
    mw.SetBinary(True)
    mw.SetFileName(TestTempDir().GetTempPath()+"mshFile.solb")
    mw.OpenSolutionFile(res)

    mw.WriteSolutionsFields(res, PointFields=[sol['SolAtVertices0']])
    mw.filePointer.write(struct.pack('i', 54))
    mw.Close()

    sol = MeshReader().ReadExtraFields(TestTempDir().GetTempPath()+"mshFile.solb")

    return 'ok'


if __name__ == '__main__':  # pragma: no cover
    from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
    BaseOutputObject.SetGlobalDebugMode(True)
    print(CheckIntegrity())  # pragma: no cover
