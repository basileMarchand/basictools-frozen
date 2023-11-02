# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
"""Stl file reader
"""

from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.IO.IOFactory import RegisterReaderClass
import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.NumpyDefs import PBasicIndexType
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanDoubleNodes


def ReadStl(fileName:str=None, string:str=None) -> UnstructuredMesh:
    """Read Stl file into a UnstructuredMesh

    Parameters
    ----------
    fileName : str, optional
        The name of the file to read, by default None
    string : str, optional
        the string to read in the case of reading from memory, by default None

    Returns
    -------
    UnstructuredMesh
        the stl surface
    """
    obj = StlReader()
    obj.SetFileName(fileName)
    obj.SetStringToRead(string)
    res = obj.Read()
    return res

class StlReader(ReaderBase):
    """Stl read class
    """
    def __init__(self, fileName=None):
        super().__init__(fileName=fileName)
        self.runCleanDoubleNodes = True

    def Read(self, fileName=None, string=None, out=None) -> UnstructuredMesh:
        """ Read a file or a string as a stl surface, ASCII and binary format are
        supported.

        Parameters
        ----------
        fileName : str, optional
            The name of the file to read, by default None
        string : str, optional
            the string to read in the case of reading from memory, by default None
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            the read stl surface
        """

        if fileName is not None:
            self.SetFileName(fileName)

        if string is not None:
            self.SetStringToRead(string)

        # check if binary
        self.readFormat = "rb"
        self.StartReading()

        header = ""
        # read the first non space characters to detect if binary or not
        while len(header) < 5:
            data = self.filePointer.read(1)
            if data[0] < 128:
                if chr(data[0]) == " ":
                    continue
                header += chr(data[0])

        if header == "solid":
            self.PrintVerbose("Ascii File")
            res = self.ReadStlAscii()
        else:
            self.PrintVerbose("Binary File")
            res = self.ReadStlBinary()

        if self.runCleanDoubleNodes:
            res.ComputeBoundingBox()
            tol = np.linalg.norm(res.boundingMax -res.boundingMin)*1e-15
            CleanDoubleNodes(res, tol=tol)
        return res

    def ReadStlBinary(self):
        # from https://en.wikipedia.org/wiki/STL_(file_format)#Binary_STL

        self.readFormat = "rb"
        self.StartReading()

        import BasicTools.Containers.UnstructuredMesh as UM

        header = self.readData(80, np.int8)
        hasMagicsColor = False
        try:
            header6 = ''.join([(chr(item) if item < 128 else " ") for item in header[0:6]])
            if header6 == "COLOR=":
                hasMagicsColor = True
            header = ''.join([(chr(item) if item < 128 else " ") for item in header])
            print("HEADER : '" + header + "'")

        except:
            pass
        nbTriangles = self.readInt32()
        print("reading  : " + str(nbTriangles) + " triangles")

        resUM = UM.UnstructuredMesh()

        #resUM.nodes = np.empty((nbTriangles*3,3), dtype=float)

        dt = np.dtype([('normal', (np.float32, 3)),
                        ('points', (np.float32, 9)),
                        ('att', (np.uint16)),
                        ])

        data = self.readData(nbTriangles, dt)
        normals = np.array(data["normal"])
        resUM.nodes = np.array(data["points"])
        resUM.nodes.shape = (nbTriangles*3, 3)

        elements = resUM.GetElementsOfType(EN.Triangle_3)
        elements.connectivity = np.array(range(resUM.GetNumberOfNodes()), dtype=PBasicIndexType)
        elements.connectivity.shape = (nbTriangles, 3)
        elements.originalIds = np.arange(nbTriangles, dtype=PBasicIndexType)
        elements.cpt = nbTriangles

        resUM.elemFields["normals"] = normals
        if hasMagicsColor:
            Color = np.array(data["att"])
            r = Color & 31
            g = (Color >> 5) & 31
            b = (Color >> 10) & 31
            active = (Color >> 15)
            resUM.elemFields["colors"] = np.vstack((r, g, b)).T
            resUM.elemFields["activecolors"] = active.astype(dtype=np.int8)

        self.EndReading()
        resUM.GenerateManufacturedOriginalIDs()
        resUM.PrepareForOutput()

        return resUM

    def ReadStlAscii(self):
        """ Read an ASCII stl file

        Returns
        -------
        UnstructuredMesh
            the read stl surface
        """
        self.readFormat = "r"
        self.StartReading()

        import BasicTools.Containers.UnstructuredMesh as UM

        resUM = UM.UnstructuredMesh()
        name = self.ReadCleanLine().split()[1]

        p = []
        normals = np.empty((0, 3), dtype=float)
        nodesbuffer = []
        while True:
            line = self.ReadCleanLine()
            if not line:
                break

            l = line.strip('\n').lstrip().rstrip()
            if l.find("facet") > -1:
                if l.find("normal") > -1:
                    normals = np.concatenate((normals, np.fromstring(
                        l.split("normal")[1], sep=" ")[np.newaxis]), axis=0)
                    continue
            if l.find("outer loop") > -1:
                for i in range(3):
                    line = self.ReadCleanLine()
                    l = line.strip('\n').lstrip().rstrip()
                    if l.find("vertex") > -1:
                        p.append(np.fromstring(l.split("vertex")[1], sep=" "))
                if len(p) == 3:
                    nodesbuffer.extend(p)
                    #resUM.nodes = np.vstack((resUM.nodes,p[0][np.newaxis,:],p[1][np.newaxis,:],p[2][np.newaxis,:]))
                    p = []
                else:  # pragma: no cover
                    print("error: outer loop with less than 3 vertex")
                    raise

        self.EndReading()

        resUM.nodes = np.array(nodesbuffer)
        del nodesbuffer
        elements = resUM.GetElementsOfType(EN.Triangle_3)
        elements.connectivity = np.array(range(resUM.GetNumberOfNodes()), dtype=PBasicIndexType)
        elements.connectivity.shape = (resUM.GetNumberOfNodes()//3, 3)
        elements.originalIds = np.arange(resUM.GetNumberOfNodes()/3, dtype=PBasicIndexType)
        elements.cpt = elements.connectivity.shape[0]
        resUM.elemFields["normals"] = normals
        resUM.GenerateManufacturedOriginalIDs()
        resUM.PrepareForOutput()

        return resUM


RegisterReaderClass(".stl", StlReader)


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
        endsolid
"""

    res = ReadStl(string=data)
    print(res)
    if res.GetNumberOfNodes() != 4:
        raise Exception()
    if res.GetNumberOfElements() != 4:
        raise Exception()

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    f = open(tempdir+"test_input_stl_data.stl", "w")
    f.write(data)
    f.close()
    res = ReadStl(fileName=tempdir+"test_input_stl_data.stl")

    from BasicTools.TestData import GetTestDataPath

    print("Binary reading")
    res1 = ReadStl(fileName=GetTestDataPath()+"coneBinary.stl")
    print(res1)

    print("Ascii reading")
    res2 = ReadStl(fileName=GetTestDataPath()+"coneAscii.stl")
    print(res2)

    if res1.GetNumberOfNodes() != res2.GetNumberOfNodes():
        raise Exception()

    # 1e-6 because the ascii file has only 6 decimals
    if np.max(abs(res1.nodes - res2.nodes)) > 1e-6:
        raise Exception()

    if res1.GetNumberOfElements() != res2.GetNumberOfElements():
        raise Exception()

    conn1 = res1.GetElementsOfType(EN.Triangle_3).connectivity
    conn2 = res2.GetElementsOfType(EN.Triangle_3).connectivity

    if not np.all(np.equal(conn1, conn2)):
        raise Exception()

    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())  # pragma: no cover
