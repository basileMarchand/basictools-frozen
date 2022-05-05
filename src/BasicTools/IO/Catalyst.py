# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import os

from BasicTools.IO.WriterBase import WriterBase as WriterBase
from BasicTools.Helpers.MPIInterface import MPIInterface

class Catalyst(WriterBase):

    bridge_initialize = False
    bridge_finalize = False
    def __init__(self):
        super(Catalyst,self).__init__()
        self.canHandleTemporal = True
        self.canHandleAppend = True
        self.canHandleMultidomain = True

        self.port = 22222
        self.hostname = "localhost"

        self.currentTime = 0
        self.currentStep = 0
        self.registrationName= "SimData"
        self.userscript = None

        self.mpiInterface= MPIInterface()
    def SetFileName(self,fileName):
        #possible options are ()
        #   None
        #   "/ignored/path/hostname:22222.catalyst"
        #   "22222.catalyst"
        #   "hostname:22222"
        #   "22222"

        if fileName is None:
            return

        if fileName.find(".catalyst") > -1 :
            fileName = fileName[0:fileName.find(".catalyst")]

        fileName = fileName.split(os.sep)[-1]

        if fileName.find(":") > -1:
            host,port = fileName.split(":")
            self.hostname = host
            self.port = int(port)
        elif len(fileName) >0 :
            self.port = int(fileName)

    def Open(self, filename=None):

        self.SetFileName(filename)

        from paraview.catalyst import bridge

        if self._isOpen:
            return

        self._isOpen = True
        if self.bridge_initialize == False:
            bridge.initialize()
            self.bridge_initialize = True

        if self.userscript is None:
            from BasicTools.Helpers.Tests import GetUniqueTempFile
            fileid, fileName  = GetUniqueTempFile(prefix='CatalystScript_',suffix=".py")
            filePointer = open(fileid,mode="w")
            filePointer.write(f"""
# File Generated by BasicTools.IO.Catalyst
from paraview.simple import *
from paraview import print_info
from paraview import catalyst

data = UnstructuredCellTypes(registrationName='{self.registrationName}')

options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = 'TimeStep'
options.CatalystLiveURL = "{self.hostname}:{self.port}"
            """)
            filePointer.close()
            bridge.add_pipeline(fileName, 2)
        else:
            # add analysis script
            bridge.add_pipeline(self.userscript, 2)


    def Close(self):
        if self.bridge_finalize == False:
            from paraview.catalyst import bridge
            bridge.finalize()
            self.bridge_finalize = True

    def Write(self,meshObject, PointFieldsNames=None, PointFields=None, CellFieldsNames=None, CellFields=None, GridFields=None, GridFieldsNames=None, Time= None, TimeStep = None):

        if not self._isOpen :
            self.Open()

        dt = 1
        if Time is not None:
            dt = Time - self.currentTime
        elif TimeStep is not None:
            dt = TimeStep
        self.currentTime += dt
        self.currentStep += 1

        from BasicTools.Containers.vtkBridge import MeshToVtk
        nodeFieldsBack = meshObject.nodeFields
        elemFieldsBack = meshObject.elemFields
        meshObject.nodeFields = { k:v for k,v in zip(PointFieldsNames,PointFields) }
        meshObject.elemFields = { k:v for k,v in zip(CellFieldsNames,CellFields) }
        vtkObj = MeshToVtk(meshObject,TagsAsFields=True)
        meshObject.nodeFields = nodeFieldsBack
        meshObject.elemFields = elemFieldsBack
        from paraview.catalyst import bridge
        bridge.coprocess(self.currentTime, self.currentStep, vtkObj, name=self.registrationName)#, wholeExtent=None)

    def __str__(self):
        res  = 'Catalyst : \n'
        res += '   hostname : '+ str(self.hostname) +'\n'
        res += '   port : '+ str(self.port) +'\n'
        return res

from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".catalyst",Catalyst)

def CheckIntegrity(GUI=False):

    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("CATALYST_NO_FAIL"): return "skip"

    numsteps = 50
    delay = 0.05
    writer = Catalyst()
    writer.Open("localhost:22221.catalyst")
    for step in range(numsteps):
        if delay > 0:
            import time
            time.sleep(delay)
        print("timestep: {0}/{1}".format(step+1, numsteps))
        # assume simulation time starts at 0
        time = step/float(numsteps)
        from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
        mesh = CreateCube(dimensions= [step+2,step+2,step+2],  origin=[-1.0,-1.0,-1.0], spacing=[1.,1.,1.], ofTetras=False)
        writer.Write(mesh, Time= time)
    writer.Close()
    return "ok"

if __name__ == "__main__":
    print(CheckIntegrity(GUI=True))