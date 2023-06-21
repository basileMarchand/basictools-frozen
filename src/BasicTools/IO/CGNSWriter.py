# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""CGNS file writer
"""
import os
import numpy as np

from BasicTools.Bridges.CGNSBridge import MeshToCGNS
from BasicTools.IO.WriterBase import WriterBase as WriterBase


class CGNSWriter(WriterBase):
    """Class to writes a CGNS file on disk
    """
    def __init__(self):
        super(CGNSWriter,self).__init__()
        self.canHandleTemporal = True
        self.canHandleAppend = False

    def __str__(self):
        res  = 'CGNSWriter'
        return res

    def Write(self, mesh, fileName, outpuPyTree = None, baseNumberOrName = 0, zoneNumberOrName = 0):
        """Function to writes a CGNS File on disk

        Parameters
        ----------
        mesh : UnstructuredMesh
            support of the data to be written
        fileName : str
            filename of the file to be read
        outpuPyTree : list
            existing pyTree in which the data structure in mesh will be appended
        baseNumberOrName : int or str, optional
            name of the base to use, by default 0 (first)
        zoneNumberOrName : int or str, optional
            name of the zone to be read, by default 0 (first)
        """

        newPyTree = MeshToCGNS(mesh, outpuPyTree, baseNumberOrName, zoneNumberOrName)

        from h5py import File, h5t
        import  h5py
        h5file = File(fileName, "w",libver=('v108', 'v108'), track_order=True)
        h5file.create_dataset(" format", data=np.array([78, 65, 84, 73, 86, 69, 0], dtype=np.int8) )
        h5file.create_dataset(" hdf5version", data=np.array([72, 68, 70, 53, 32, 86, 101, 114, 115, 105, 111, 110, 32, 49, 46, 49, 50, 46, 49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],dtype=np.int8) )

        buff_S33 = np.empty(1, '|S33')

        def ToStrPadded(string, n):
            return np.bytes_(bytes(string, "ascii").ljust(n, b"\x00"))

        def WriteToH5File(file, node):
            tid = h5t.C_S1.copy()
            tid.set_size(33)
            H5T_C_S1_33 = h5py.Datatype(tid)

            tid = h5t.C_S1.copy()
            tid.set_size(3)
            H5T_C_S1_3 = h5py.Datatype(tid)
            if node[3] == "CGNSTree_t":
                file.attrs.create("name", "HDF5 MotherNode" , dtype=H5T_C_S1_33 )
                file.attrs.create("label", "Root Node of HDF5 File" , dtype=H5T_C_S1_33 )
            else:
                file.attrs.create("name", node[0] , dtype=H5T_C_S1_33 )
                file.attrs.create("label",node[3] , dtype=H5T_C_S1_33 )
                file.attrs["flags"] = np.array([0], dtype=np.int32)

            if node[1] is None:
                file.attrs.create("type","MT" , dtype=H5T_C_S1_3 )

            else:
                numpyTypesToHdF5Types = {np.dtype("float32"):"R4", np.dtype("float64"):"R8", np.dtype("int32"):"I4", np.dtype("int64"):"I8", np.dtype("|S1"):"C1"}
                data = node[1]
                if node[1].dtype in numpyTypesToHdF5Types:
                    file.attrs.create("type",ToStrPadded(numpyTypesToHdF5Types[node[1].dtype], 3) , dtype=H5T_C_S1_3 )
                else:
                    raise Exception("For the moment cgns does not support string field")

                if node[1].dtype == np.dtype("|S1") :
                    idx1 = np.nonzero(data != b"")
                    idx2 = np.nonzero(data == b"")
                    tmp = data.copy()
                    data = np.zeros(data.shape, dtype=np.int8)
                    data[idx1] = np.vectorize(ord)(tmp[idx1])
                    data[idx2] = 0
                    data = data.astype(np.int8)

                data = data.transpose()
                file.create_dataset(" data", data=data)

            for child in node[2]:
                sgrp = file.create_group(child[0],track_order=True)
                WriteToH5File(sgrp, child)


        WriteToH5File(h5file,newPyTree)

def CheckIntegrity():

    try:
        import h5py
    except:
        return "skip! h5py module not available"

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    reader.atIntegrationPoints = False

    import BasicTools.IO.GeofReader as GR
    myMesh = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")

    reader.atIntegrationPoints = False
    for dn in myMesh.nodeFields:
        myMesh.nodeFields[dn] = reader.ReadField(fieldname=dn, timeIndex=1)

    myMesh.nodeFields["Nodes arange"] = np.arange(myMesh.GetNumberOfNodes(),dtype=float)
    myMesh.elemFields["Element arange"] = np.arange(myMesh.GetNumberOfElements(),dtype=float)
    from BasicTools.Containers.Tags import Tags

    ##################################
    # EXEMPLE SYNTAXE DU WRITER
    import BasicTools.IO.CGNSWriter as CW
    CgW = CW.CGNSWriter()
    CgW.Write(mesh = myMesh, fileName = tempdir+os.sep+"toto.cgns", baseNumberOrName = 0, zoneNumberOrName = 0)
    ##################################

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
