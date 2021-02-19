# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

# this files is inteded to be used inside paraview as a plugin
# compatible with paraview 5.7+

import copy

import numpy as np

from paraview.util.vtkAlgorithm import *
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

debug = False

try:
    from BasicTools.IO.IOFactory import InitAllReaders
    from BasicTools.IO.IOFactory import ReaderFactory
    from BasicTools.IO.IOFactory import InitAllWriters
    from BasicTools.IO.IOFactory import WriterFactory
    from BasicTools.Containers.vtkBridge import VtkToMesh, MeshToVtk


    paraview_plugin_name = "BasicTools ParaView Bridge"
    paraview_plugin_version = "5.7"

    #----------------------------- The Readers ------------------------------------
    InitAllReaders()


    for pext in ReaderFactory.keys():
        if pext==".PIPE":
            continue
        ext = pext[1:]
        readerClass = ReaderFactory.GetClass(pext)
        readerClassName = readerClass.__name__
        wrapperClassName = "BasicToolsPython"+ext.upper()+"Reader"

        def GetInit(readerClass):
            def myinit(self):
                VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkUnstructuredGrid')
                self._filename = None
                self.basicToolsReader = readerClass()
            return myinit

        @smproperty.stringvector(name="FileName")
        @smdomain.filelist()
        @smhint.filechooser(extensions=ext, file_description=ext + " files")
        def SetFileName(self, name):
            """Specify filename for the file to read."""
            if self._filename != name:
                self._filename = name
                self.Modified()
                self.basicToolsReader.SetFileName(name)
                if name is not None:
                    self.GetTimestepValues()

        @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
        def GetTimestepValues(self):
            if self._filename is None or self._filename == "None":
                return None

            if self.basicToolsReader.canHandleTemporal:
                self.basicToolsReader.SetFileName(self._filename)
                self.metadata = self.basicToolsReader.ReadMetaData()

                timeSteps = self.basicToolsReader.GetAvilableTimes()
                if len(timeSteps) == 0:
                    return None
                return timeSteps
            else:
                return None

        @smproperty.xml ("""
          <IntVectorProperty name="Tags As Fields"
                             command="SetTagsAsFields"
                             number_of_elements="1"
                             default_values="1">
              <BooleanDomain name="bool"/>
          </IntVectorProperty>""")
        def SetTagsAsFields(self, val):
            self.__TagsAsFields = bool(val)
            self.Modified()


        def RequestInformation(self, request, inInfoVec, outInfoVec):
            executive = self.GetExecutive()
            outInfo = outInfoVec.GetInformationObject(0)
            outInfo.Remove(executive.TIME_STEPS())
            outInfo.Remove(executive.TIME_RANGE())

            timesteps = self.GetTimestepValues()
            if timesteps is not None:
                for t in timesteps:
                    outInfo.Append(executive.TIME_STEPS(), t)
                outInfo.Append(executive.TIME_RANGE(), timesteps[0])
                outInfo.Append(executive.TIME_RANGE(), timesteps[-1])
            return 1


        def RequestData(self, request, inInfoVec, outInfoVec):
            reader =  self.basicToolsReader

            if reader.canHandleTemporal :

                outInfo = outInfoVec.GetInformationObject(0)
                executive = self.GetExecutive()
                #timesteps = self.GetTimestepValues()

                if outInfo.Has(executive.UPDATE_TIME_STEP()):
                    time = outInfo.Get(executive.UPDATE_TIME_STEP())
                else:
                    time = 0


                reader.SetTimeToRead(time= time,timeIndex=None)

            reader.SetFileName(self._filename)
            mesh = reader.Read()


            from BasicTools.Containers.vtkBridge import MeshToVtk

            MeshToVtk(mesh, vtkobject=vtkUnstructuredGrid.GetData(outInfoVec,0),TagsAsFields=self.__TagsAsFields)
            return 1

        obj = type(wrapperClassName,
                  (VTKPythonAlgorithmBase,),
                  {"__init__":GetInit(readerClass),
                   "SetFileName": SetFileName,
                   "GetTimestepValues": GetTimestepValues,
                   "SetTagsAsFields":SetTagsAsFields,
                   "RequestData":RequestData,
                   "RequestInformation":RequestInformation}
                   )


        obj2 = smproxy.reader(name=wrapperClassName, label="BasicTools Python-based "+ ext.upper() +" Reader",
                    extensions=ext, file_description=ext+" files [BasicTools Reader]")(obj)

        locals()[wrapperClassName] = obj2


    #----------------------------- The Writers ------------------------------------



    InitAllWriters()

    for pext in WriterFactory.keys():
        if pext==".PIPE":
            continue
        ext = pext[1:]

        writerClass = WriterFactory.GetClass(pext)
        writerClassName = writerClass.__name__
        wrapperClassName = "BasicToolsPython"+ext.upper()+"Writer"

        def GetInit(readerClass):
            def myinit(self):
                VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=0, inputType='vtkUnstructuredGrid')
                self._filename = None
                self._binary = True
                self.basicToolsReader = readerClass()
            return myinit


        @smproperty.stringvector(name="FileName", panel_visibility="never")
        @smdomain.filelist()
        def SetFileName(self, fname):
            """Specify filename for the file to write."""
            if self._filename != fname:
                self._filename = fname
                self.Modified()
                self.basicToolsReader.SetFileName(fname)


        @smproperty.xml ("""
          <IntVectorProperty name="In Binary (if avilable)"
                             command="SetBinary"
                             number_of_elements="1"
                             default_values="1">
              <BooleanDomain name="bool"/>
          </IntVectorProperty>""")
        def SetBinary(self, val):
            self._binary = bool(val)
            self.Modified()


        def RequestData(self, request, inInfoVec, outInfoVec):
            from vtkmodules.vtkCommonDataModel import vtkTable
            from vtkmodules.numpy_interface import dataset_adapter as dsa

            vtkdata = vtkUnstructuredGrid.GetData(inInfoVec[0], 0)
            mesh = VtkToMesh(vtkdata)

            from BasicTools.IO.IOFactory import CreateWriter


            writer = self.basicToolsReader

            writer.SetFileName(self._filename)
            writer.SetBinary(self._binary)
            writer.Open()

            PointFields = None
            PointFieldsNames = None
            if hasattr(mesh,"nodeFields"):
                PointFieldsNames = list(mesh.nodeFields.keys())
                PointFields = list(mesh.nodeFields.values())

            CellFields = None
            CellFieldsNames = None
            if hasattr(mesh,"elemFields"):
                CellFieldsNames = list(mesh.elemFields.keys())
                CellFields = list(mesh.elemFields.values())

            writer.Write(mesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields )
            writer.Close()

            return 1

        def Write(self):
            self.Modified()
            self.Update()

        writerInstanse = writerClass()
        props = {"__init__":GetInit(writerClass),
                   "SetFileName": SetFileName,
                   "RequestData":RequestData,
                   "Write":Write}

        if writerInstanse.canHandleBinaryChange:
            props["SetBinary"] = SetBinary


        obj = type(wrapperClassName,
                  (VTKPythonAlgorithmBase,),
                  props)

        #print("Python"+classname+"Reader", obj )

        obj1 = smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)(obj)
        obj2 = smproperty.input(name="Input", port_index=0)(obj1)
        obj3 = smproxy.writer(extensions=ext, file_description=ext+" files [BasicTools Writer]", support_reload=False)(obj2)

        locals()[wrapperClassName] = obj3


    #----------------------------- The Filters ------------------------------------
    # this is experimental
    @smproxy.filter(name="Center DataSet")
    @smhint.xml("""<ShowInMenu category="BasicTools Filters" />""")
    @smproperty.input(name="Target mesh", port_index=1)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    @smproperty.input(name="Data to move", port_index=0)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    class CenterSecondObjectOnTheFirst(VTKPythonAlgorithmBase):
        def __init__(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1, outputType="vtkUnstructuredGrid")

        def FillInputPortInformation(self, port, info):
            if port == 0:
                info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
                info.Append(self.INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData")
            else:
                info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
                info.Append(self.INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData")
            return 1

        def RequestData(self, request, inInfoVec, outInfoVec):
            from vtkmodules.vtkCommonDataModel import  vtkDataSet
            input0 = VtkToMesh(vtkDataSet.GetData(inInfoVec[0], 0))
            input1 = VtkToMesh(vtkDataSet.GetData(inInfoVec[1], 0))

            mean0 = np.sum(input0.nodes,axis=0)/input0.GetNumberOfNodes()
            mean1 = np.sum(input1.nodes,axis=0)/input1.GetNumberOfNodes()

            output = vtkUnstructuredGrid.GetData(outInfoVec, 0)

            # the user must not modify the inputs
            outputmesh = copy.copy(input0)
            outputmesh.nodes = input0.nodes + (mean1 - mean0)
            MeshToVtk(outputmesh, output)

            return 1


    print("BasicTools ParaView Plugin Loaded")
except:
    print("Error loading BasicTools ParaView Plugin")
    print("BasicTools in the PYTHONPATH ??? ")
    if debug:
        raise
