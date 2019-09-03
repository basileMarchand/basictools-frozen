# -*- coding: utf-8 -*-
import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Containers.vtkBridge import VtkToMesh

def LoadVtuWithVTK(filename):
    import vtk
    readerVTU = vtk.vtkXMLUnstructuredGridReader()
    readerVTU.SetFileName(filename)
    readerVTU.Update()

    data = readerVTU.GetOutput()

    if data.GetNumberOfPoints() == 0:
        raise ValueError("No point data could be loaded from '" + filename)# pragma: no cover

    return data


class VtuReader(ReaderBase):
    def __init__(self,fileName = None):
        super(VtuReader,self).__init__()

    def Read(self, fileName=None,string=None,out=None):

      if fileName is not None:
          self.SetFileName(fileName)

      self.output = VtkToMesh(LoadVtuWithVTK(self.fileName))
      return self.output


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".vtu",VtuReader)

def CheckIntegrity():
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
