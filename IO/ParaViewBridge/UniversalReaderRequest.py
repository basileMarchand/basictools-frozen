# to deal with time
# https://blog.kitware.com/defining-time-varying-sources-with-paraviews-programmable-source/

from BasicTools.IO.IOFactory import CreateReader
from BasicTools.IO.IOFactory import InitAllReaders

InitAllReaders()


reader = CreateReader("."+filename.split(".")[-1])
if reader.canHandleTemporal :
   reader.SetFileName(filename)
   self.metadata = reader.ReadMetaData()

   timeSteps = reader.GetAvilableTimes()

   outInfo = self.GetOutputInformation(0)
   timeRange = [timeSteps[0], timeSteps[-1]]
   outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), timeRange, 2)
   outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(), timeSteps, len(timeSteps))
