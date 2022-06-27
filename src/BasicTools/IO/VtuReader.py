# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Bridges.vtkBridge import VtkToMesh
from BasicTools.IO.IOFactory import RegisterReaderClass

class VtkReader(ReaderBase):
    def __init__(self,fileName = None):
        super(VtkReader,self).__init__(fileName=fileName)

    def Read(self, fileName=None):
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
      return self.output

RegisterReaderClass(".vtk",VtkReader)

def LoadVtuWithVTK(filename):
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
    def __init__(self,fileName = None):
        super(VtuReader,self).__init__(fileName =fileName)

    def Read(self, fileName=None,string=None,out=None):

      if fileName is not None:
          self.SetFileName(fileName)

      self.output = VtkToMesh(LoadVtuWithVTK(self.fileName))
      return self.output

RegisterReaderClass(".vtu",VtuReader)

def CheckIntegrity():
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
