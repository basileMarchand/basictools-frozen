# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


"""Vtu and Vtk file reader
"""

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
        """Function that performs the reading of a vtk file

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
        from  vtkmodules.vtkIOLegacy import vtkGenericDataObjectReader

        if fileName is not None:
            self.SetFileName(fileName)

        reader = vtkGenericDataObjectReader()
        reader.SetFileName(self.fileName)
        reader.Update()
        res = reader.GetUnstructuredGridOutput()
        if res is None:
            res = reader.GetPolyDataOutput()
        self.output = VtkToMesh(res)
        self.output.ConvertDataForNativeTreatment()
        return self.output

RegisterReaderClass(".vtk",VtkReader)

def LoadVtuWithVTK(filename):
    """Function API for reading a vtu file using vtk

    Parameters
    ----------
    filename : str
        name of the file to be read

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, vtkPolyData
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

    readerVTU = vtkXMLUnstructuredGridReader()
    readerVTU.SetFileName(filename)
    readerVTU.Update()

    data = readerVTU.GetOutput()

    if data.GetNumberOfPoints() == 0:
        raise ValueError("No point data could be loaded from '" + filename)# pragma: no cover

    return data

class VtuReader(ReaderBase):
    """Vtu Reader class
    """
    def __init__(self,fileName = None):
        super(VtuReader,self).__init__(fileName =fileName)

    def Read(self, fileName=None,string=None,out=None):
        """Function that performs the reading of a vtu file

        Parameters
        ----------
        fileName : str, optional
            name of the file to be read, by default None
        string : str, optional
            data to be read as a string instead of a file, by default None
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        if fileName is not None:
            self.SetFileName(fileName)

        self.output = VtkToMesh(LoadVtuWithVTK(self.fileName))
        self.output.ConvertDataForNativeTreatment()

        return self.output

RegisterReaderClass(".vtu",VtuReader)

def CheckIntegrity():
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
