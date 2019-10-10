#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# to deal with time
# https://blog.kitware.com/defining-time-varying-sources-with-paraviews-programmable-source/

import os
path = self.GetPythonPath().strip("'")

import subprocess, os
pp = os.environ.get("PYTHONPATH","").split(":")
if not(path in pp):
   os.environ["PYTHONPATH"] = path + ":" + os.environ.get("PYTHONPATH","")


from BasicTools.IO.IOFactory import CreateReader
from BasicTools.IO.IOFactory import InitAllReaders

InitAllReaders()

try:
    self.reader
except:
    self.reader = CreateReader("."+filename.split(".")[-1])

reader = self.reader

if reader.canHandleTemporal :
   reader.SetFileName(filename)
   self.metadata = reader.ReadMetaData()

   timeSteps = reader.GetAvilableTimes()

   outInfo = self.GetOutputInformation(0)
   timeRange = [timeSteps[0], timeSteps[-1]]
   outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), timeRange, 2)
   outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(), timeSteps, len(timeSteps))
