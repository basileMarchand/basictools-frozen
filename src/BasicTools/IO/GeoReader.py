# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Geo file reader (Zset mesh file)
"""

import numpy as np
import struct

from BasicTools.IO.ReaderBase import ReaderBase
import BasicTools.Containers.UnstructuredMesh as UM
from BasicTools.IO.ZsetTools import GeofNumber,PermutationZSetToBasicTools, nbIntegrationsPoints
from BasicTools.NumpyDefs import PBasicIndexType

def ReadGeo(fileName=None,out=None,readElset=True,readFaset=True):
    """Function API for reading a geo mesh file

    Parameters
    ----------
    fileName : str, optional
        name of the file to be read, by default None
    out : UnstructuredMesh, optional
        output unstructured mesh object containing reading result, by default None
    readElset : bool, optional
        if False, ignores the elset informations, by default True
    readFaset : bool, optional
        if False, ignores the faset informations, by default True

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    reader = GeoReader()
    reader.SetFileName(fileName)
    return reader.Read(fileName=fileName, out=out,readElset=readElset,readFaset=readFaset)

def ReadMetaData(fileName=None):
    """Function API for reading the metadata of a geo mesh file

    Parameters
    ----------
    fileName : str, optional
        name of the file to be read, by default None

    Returns
    -------
    dict
        global information on the mesh to read
    """
    reader = GeoReader()
    reader.SetFileName(fileName)
    return reader.ReadMetaData()

class GeoReader(ReaderBase):
    """Geo Reader class
    """
    def __init__(self):
        super(GeoReader,self).__init__()
        self.readFormat = 'rb'

    def _getTag(self):
        nb = struct.unpack(">i", self.rawread(4))[0]
        rawtag = self.rawread(nb)
        return rawtag[0:-1].decode("utf-8")

    def readMagic(self):
            magic = self.rawread(13)
            if not (b'Z7BINARYGEOF\x00' == magic ) :
                raise(Exception("Bad file"))

    def readInt32(self):
        return struct.unpack(">i", self.rawread(4))[0]

    def readInts32(self,cpt):
        return np.array(struct.unpack(">"+"i"*cpt, self.rawread(4*cpt)))

    def ReadMetaData(self):
        """Function that performs the reading of the metadata of a geof mesh file

        Returns
        -------
        dict
            global information on the mesh to read
        """
        return self.Read(onlyMeta = True)

    def Read(self, fileName=None,out=None,readElset=True,readFaset=True, onlyMeta = False):
        """Function that performs the reading of a geo mesh file

        Parameters
        ----------
        fileName : str, optional
            name of the file to be read, by default None
        out : UnstructuredMesh, optional
            Not Used, by default None
        readElset : bool, optional
            Not Used, by default None
        readFaset : bool, optional
            Not Used, by default None
        onlyMeta : bool, optional
            if True, only read metadata, by default False

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        if fileName is not None:
            self.SetFileName(fileName)

        self.StartReading()

        res = UM.UnstructuredMesh()
        metadata = {}

        FENames = {}
        oidToElementContainer = {}
        oidToLocalElementNumber = {}

        cpt = 0
        while True:
            self.readMagic()
            tag = self._getTag()
            if tag == u"node":
                nbNodes = self.readInt32()
                metadata["nbNodes"] = int(nbNodes)
                dims = self.readInt32()
                metadata["dimensionality"] = int(dims)

                if onlyMeta :
                    self.filePointer.seek((dims*8+4)*nbNodes,1 )
                else:
                    res.nodes = np.zeros((nbNodes,3))
                    res.originalIDNodes = np.empty((nbNodes,),dtype=PBasicIndexType)

                    for i in range(nbNodes):
                        data = self.rawread((dims*8+4))
                        data = struct.unpack((">i"+"d"*dims), data)
                        res.originalIDNodes[i] = data[0]
                        res.nodes[i,0:dims] = data[1:]

            elif tag == u"element":
                n_elem = self.readInt32()
                metadata['nbElements'] = int(n_elem)
                IPPerElement = np.empty(metadata['nbElements'],dtype=PBasicIndexType)

                n_grp = self.readInt32()
                for i in range(n_grp):
                    n_in_grp = self.readInt32()
                    ltype =  self._getTag()
                    nametype = GeofNumber[ltype]
                    if not onlyMeta:
                        elements = res.GetElementsOfType(nametype)
                        elementsInContainerCpt = elements.GetNumberOfElements()
                        elements.Allocate(n_in_grp+elementsInContainerCpt)

                    n_node = self.readInt32()

                    if nametype not in FENames:
                        FENames[nametype] = []

                    FENames_etype = FENames[nametype]

                    nintegpoints =  nbIntegrationsPoints[ltype]
                    if onlyMeta :
                        self.filePointer.seek((4)*(n_node+2)*n_in_grp,1 )
                        for j in range(n_in_grp):
                            IPPerElement[cpt] = nintegpoints
                            cpt += 1
                    else:
                        perm = None
                        if ltype in PermutationZSetToBasicTools:
                            perm = PermutationZSetToBasicTools[ltype]


                        for j in range(n_in_grp):
                            idd = self.readInt32()
                            rank = self.readInt32()
                            conn = self.readInts32(n_node)
                            IPPerElement[cpt] = nintegpoints
                            if onlyMeta:
                                cpt += 1
                                continue

                            if perm:
                                conn =  [conn[x] for x in perm ]
                            elements.connectivity[j+elementsInContainerCpt,:] = conn
                            elements.originalIds[j+elementsInContainerCpt] = rank
                            FENames[nametype].append(ltype)
                            oidToElementContainer[rank] = elements
                            oidToLocalElementNumber[rank] = j+elementsInContainerCpt

                            cpt += 1

            elif tag == u"nset":
                n_nset = self.readInt32()
                for i in range(n_nset):
                    name = self._getTag()
                    size = self.readInt32()

                    if  onlyMeta :
                        self.filePointer.seek(size*4,1)
                    else:
                        ids = np.fromfile(self.filePointer,count=size, dtype=np.int32).byteswap()
                        res.nodesTags.CreateTag(name).SetIds(ids)

            elif tag == u"elset":
                n_nset = self.readInt32()
                for i in range(n_nset):
                    name = self._getTag()
                    size = self.readInt32()

                    if onlyMeta :
                        self.filePointer.seek(size*4,1)
                    else:
                        ids = np.fromfile(self.filePointer,count=size, dtype=np.int32).byteswap()
                        for oid in ids:
                            oidToElementContainer[oid].tags.CreateTag(name,False).AddToTag(oidToLocalElementNumber[oid])

            elif tag == u"bset":
                n_bset = self.readInt32()
                for i in range(n_bset):
                    bsetName = self._getTag()
                    eltype = self._getTag()
                    n_el = self.readInt32()

                    if eltype == "faset":
                        if onlyMeta :
                            for j in range(n_el):
                                how = self.readInt32()
                                self.filePointer.seek(how*4+5,1)
                        else:
                            for j in range(n_el):
                                how = self.readInt32()
                                conn = np.fromfile(self.filePointer,count=how, dtype=np.int32).byteswap()
                                rawname = self.rawread(5)
                                bsetelemtype = rawname[0:2].decode("utf-8")

                                nametype = GeofNumber[bsetelemtype]

                                perm = None
                                if bsetelemtype in PermutationZSetToBasicTools:
                                    conn =  [conn[x] for x in PermutationZSetToBasicTools[bsetelemtype] ]

                                elements = res.GetElementsOfType(nametype)
                                localId = elements.AddNewElement(conn,-1)
                                elements.GetTag(bsetName).AddToTag(localId-1)
                    elif eltype == "liset":
                        if onlyMeta :
                            for j in range(n_el):
                                how = self.readInt32()
                                self.filePointer.seek(how*4+5,1)
                        else:
                            for j in range(n_el):
                                how = self.readInt32()
                                conn = np.fromfile(self.filePointer,count=how, dtype=np.int32).byteswap()
                                rawname = self.rawread(5)
                                bsetelemtype = rawname[0:-1].decode("utf-8")
                                nametype = GeofNumber[bsetelemtype]

                                perm = None
                                if bsetelemtype in PermutationZSetToBasicTools:
                                    conn =  [conn[x] for x in PermutationZSetToBasicTools[bsetelemtype] ]

                                elements = res.GetElementsOfType(nametype)
                                localId = elements.AddNewElement(conn,-1)
                                elements.GetTag(bsetName).AddToTag(localId-1)


                    else:
                        raise(Exception("Error I dont know how to treat bset of type "   + str(eltype) ))

            elif tag == "ipset":
                d = self.readInt32()
                for i in range(d):
                    name = self.readTag()
                    l = self.readInt32()
                    for j in range(l):
                        nb = self.readInt32()
                        nb = self.readInt32()
                        self.filePointer.seek(4*nb)
            elif tag == "return":
                if onlyMeta:
                    metadata['nbIntegrationPoints'] = np.sum(IPPerElement)
                    metadata['IPPerElement'] = IPPerElement
                    self.EndReading()
                    res.PrepareForOutput()
                    return metadata
                else:
                    self.EndReading()
                    res.PrepareForOutput()

                    fenames = []
                    for elname in res.elements:
                        if elname not in FENames:
                            fenames.extend(["NA"]*res.elements[elname].GetNumberOfElements())
                        else:
                            fenames.extend(FENames[elname])

                    res.elemFields["FE Names"] = np.array(fenames,dtype=np.str_)
                    return res
            else:
                print(res)
                raise(Exception("Tag '"+ str(tag )+ "' not treated"))

        self.EndReading()
        res.PrepareForOutput()


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".geo",GeoReader)


def CheckIntegrity():

    from BasicTools.TestData import GetTestDataPath
    fileName  = GetTestDataPath()+"cube.geo"

    reader = GeoReader()
    reader.SetFileName(fileName)
    print(reader.ReadMetaData())
    print(ReadGeo(fileName=fileName))

    print(ReadMetaData(fileName=fileName))
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
