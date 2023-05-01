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

try:
    from CGNS.MAP import save
    cgnsLoaded = True
except:
    cgnsLoaded = False

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
        save(fileName, newPyTree)


def CheckIntegrity():

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

    myMesh.elemFields["arange"] = np.arange(myMesh.GetNumberOfElements())


    ##################################
    # EXEMPLE SYNTAXE DU WRITER
    import BasicTools.IO.CGNSWriter as CW
    CgW = CW.CGNSWriter()
    CgW.Write(mesh = myMesh, fileName = tempdir+os.sep+"toto.cgns", baseNumberOrName = 0, zoneNumberOrName = 0)
    ##################################

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
