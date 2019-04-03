#from BasicTools.IO.UniversalReader import ReadMeshAndPopulateVtkObject as ReadMeshAndPopulateVtkObject
#ReadMeshAndPopulateVtkObject(filename,self.GetOutput(),TagsAsFields=TagsAsFields)


from BasicTools.IO.IOFactory import CreateReader
from BasicTools.IO.IOFactory import InitAllReaders
from BasicTools.Containers.vtkBridge import MeshToVtk

InitAllReaders()

#reader = CreateReader("."+filename.split(".")[-1])

reader = self.reader

if reader.canHandleTemporal :

   outInfo = self.GetOutputInformation(0)
   if outInfo.Has(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()):
     time = outInfo.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())
   else:
     time = 0
   reader.SetTimeToRead(time)

reader.SetFileName(filename)
mesh = reader.Read()



MeshToVtk(mesh, self.GetOutput(),TagsAsFields=TagsAsFields)



