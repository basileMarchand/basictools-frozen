# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


"""Vtu and Vtk file reader
"""
import os
import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Bridges.vtkBridge import VtkToMesh
from BasicTools.IO.IOFactory import RegisterReaderClass

class VtkReader(ReaderBase):
    """Vtk Reader class
    """
    def __init__(self,fileName = None) -> None:
        super().__init__(fileName=fileName)

    def Read(self, fileName:str=None) -> UnstructuredMesh:
        """Function that performs the reading mesh using vtk,
        the reader is selected based on the file_extension.
        Current files supported are vtk, stl

        Parameters
        ----------
        fileName : str, optional
            name of the file to be read, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, vtkPolyData
        from vtkmodules.vtkIOLegacy import vtkGenericDataObjectReader
        from vtkmodules.vtkIOGeometry import vtkSTLReader

        VTK_ReaderByExtention ={ "stl": vtkSTLReader}


        if fileName is not None:
            self.SetFileName(fileName)

        filename, file_extension = os.path.splitext(self.fileName)
        reader = VTK_ReaderByExtention.get(file_extension, vtkGenericDataObjectReader)()
        reader.SetFileName(self.fileName)
        reader.Update()
        res = reader.GetUnstructuredGridOutput()
        if res is None:
            res = reader.GetPolyDataOutput()

        self.output = VtkToMesh(res)
        self.output.ConvertDataForNativeTreatment()
        return self.output

try:
    import vtkmodules.vtkIOLegacy
    RegisterReaderClass(".vtk",VtkReader)
    RegisterReaderClass(".vti",VtkReader)
    RegisterReaderClass(".vtp",VtkReader)
    RegisterReaderClass(".vtr",VtkReader)
    RegisterReaderClass(".vts",VtkReader)
    RegisterReaderClass(".vtu",VtkReader)
except:
    pass

def CheckIntegrity():
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
