# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

# this files is inteded to be used inside paraview as a plugin
# compatible with paraview 5.7+

import copy
import time

_startTime = time.time()
debug = False
if debug :
    def PrintDebug(mes):
        import time
        print(mes,time.time() - _startTime)
else:
    def PrintDebug(mes):
        pass

import numpy as np

from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain, smhint
from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

try:
    PrintDebug("Loading libs")
    from BasicTools.Containers.vtkBridge import GetInputVtk, GetInputBasicTools,  SetOutputBasicTools
    from BasicTools.IO.IOFactory import InitAllReaders
    from BasicTools.IO.IOFactory import ReaderFactory
    from BasicTools.IO.IOFactory import InitAllWriters
    from BasicTools.IO.IOFactory import WriterFactory
    from BasicTools.Containers.vtkBridge import VtkToMesh, VtkToMeshOnlyMeta
    import BasicTools.Containers.ElementNames as EN
    PrintDebug("Loading")

    try :
        PrintDebug("loading meshio interface Done")
        import BasicTools.Containers.MeshIOBridge as MeshIOBridge
        MeshIOBridge.InitAllReaders()
        MeshIOBridge.InitAllWriters()
        MeshIOBridge.AddReadersToBasicToolsFactory()
        MeshIOBridge.AddWritersToBasicToolsFactory()
        PrintDebug("loading meshio interface Done")
    except:
        print("Error loading meshio interface ")
        if debug:
            raise

    paraview_plugin_name = "BasicTools ParaView Bridge"
    paraview_plugin_version = "5.7"

    #----------------------------- The Readers ------------------------------------
    PrintDebug("Init readers ")
    InitAllReaders()
    PrintDebug("Init readers Done")

    for pext,readerClass,_ in ReaderFactory.AllEntries():
        if pext==".PIPE":
            continue
        ext = pext[1:]
        #readerClass = ReaderFactory.GetClass(pext)
        readerClassName = readerClass.__name__
        wrapperClassName = "BasicToolsPython"+ext.upper()+"Reader_"+readerClassName

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

                timeSteps = self.basicToolsReader.GetAvailableTimes()
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
        <Documentation>
          This property indicates if tags (points/cells) must be conerted to
          a 0/1 char field
        </Documentation>
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
            SetOutputBasicTools(request,inInfoVec,outInfoVec,mesh,TagsAsFields=self.__TagsAsFields)
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

        name = " files [BasicTools Reader]"
        if readerClassName.find("MeshIO") >=0 :
            name = " files [BasicTools MeshIO Reader]"

        obj2 = smproxy.reader(name=wrapperClassName, label="BasicTools Python-based "+ ext.upper() +" Reader",
                    extensions=ext, file_description=ext+name)(obj)

        locals()[wrapperClassName] = obj2


    #----------------------------- The Writers ------------------------------------


    PrintDebug("InitAllWriters")
    InitAllWriters()
    PrintDebug('InitAllWriters Done')
    for pext in WriterFactory.keys():
        if pext==".PIPE":
            continue
        ext = pext[1:]

        writerClass = WriterFactory.GetClass(pext)
        writerClassName = writerClass.__name__
        wrapperClassName = "BasicToolsPython"+ext.upper()+"Writer_"+writerClassName

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
        <Documentation>
          This property indicates if a binary version of the format must be used.
        </Documentation>
          </IntVectorProperty>""")
        def SetBinary(self, val):
            self._binary = bool(val)
            self.Modified()


        def RequestData(self, request, inInfoVec, outInfoVec):
            mesh = GetInputBasicTools(request,inInfoVec,outInfoVec,FieldsAsTags=True,connectino=0,port=0)

            from BasicTools.IO.IOFactory import CreateWriter
            writer = self.basicToolsReader

            writer.SetFileName(self._filename)

            if writer.canHandleBinaryChange:
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
        name = ext+" files [BasicTools Writer]"
        if writerClassName.find("MeshIO") >=0 :
            name = ext+" files [BasicTools MeshIO Writer]"

        obj1 = smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)(obj)
        obj2 = smproperty.input(name="Input", port_index=0)(obj1)
        obj3 = smproxy.writer(extensions=ext, file_description=name, support_reload=False)(obj2)

        locals()[wrapperClassName] = obj3

    #----------------------------- The Filters ------------------------------------
    # this is experimental
    @smproxy.filter(name="Prefix Filter")
    @smhint.xml("""<ShowInMenu category="BasicTools" />""")
    @smproperty.input(name="Input", port_index=0)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    class PrefixFilter(VTKPythonAlgorithmBase):
        def __init__(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
            self.prefix = ""
            self.invert = False
            self.onCellTags = False
            self.onPointTags = False
            self.onCellFields = False
            self.onPointFields = False

        def RequestData(self, request, inInfoVec, outInfoVec):
            inMesh = GetInputBasicTools(request, inInfoVec, outInfoVec, True)
#for elem in :
#    print elem
            if self.onPointTags:
                tagsToDelete = list(filter(lambda x: x.find(self.prefix) == -1, inMesh.nodesTags.keys()))
                print("tagsToDelete",tagsToDelete)
                inMesh.nodesTags.DeleteTags(tagsToDelete)

            if self.onCellTags:
                for elemtype, data in inMesh.elements.items():
                    tagsToDelete = list(filter(lambda x: x.find(self.prefix) == -1, data.tags.keys()))
                    print(elemtype +" tagsToDelete",tagsToDelete)
                    data.tags.DeleteTags(tagsToDelete)


            if self.onPointFields:
                inMesh.nodeFields =dict(filter(lambda x: x[0].find(self.prefix) == 0, inMesh.nodeFields.items()))

            if self.onCellFields:
                inMesh.elemFields =dict(filter(lambda x: x[0].find(self.prefix) == 0, inMesh.elemFields.items()))


                #inMesh.nodeFields = {(k if (k.find(self.prefix) == 0) ):v for k,v in inMesh.nodeFields.items()  }
            SetOutputBasicTools(request,inInfoVec,outInfoVec,inMesh)
            return 1

        @smproperty.xml("""<IntVectorProperty name="OnCellTags" command="SetOnCellTags"
                             number_of_elements="1" default_values="0">
                             <BooleanDomain name="bool"/>
        <Documentation>
          This property indicates if the prefix must be apply to CellTags
        </Documentation>
                             </IntVectorProperty>""")
        def SetOnCellTags(self, val):
            val = int(val)
            if self.onCellTags != val :
                self.onCellTags = val
                self.Modified()

        @smproperty.xml("""<IntVectorProperty name="OnPointTags" command="SetOnPointTags"
                             number_of_elements="1" default_values="0">
                             <BooleanDomain name="bool"/>
        <Documentation>
          This property indicates if the prefix must be apply to PointTags
        </Documentation>
                             </IntVectorProperty>""")
        def SetOnPointTags(self, val):
            val = int(val)
            if self.onPointTags != val :
                self.onPointTags = val
                self.Modified()

        @smproperty.xml("""<IntVectorProperty name="OnCellFields" command="SetOnCellFields"
                             number_of_elements="1" default_values="0">
                             <BooleanDomain name="bool"/>
        <Documentation>
          This property indicates if the prefix must be apply to CellFields
        </Documentation>
                             </IntVectorProperty>""")
        def SetOnCellFields(self, val):
            val = int(val)
            if self.onCellFields != val :
                self.onCellFields = val
                self.Modified()

        @smproperty.xml("""<IntVectorProperty name="OnPointFields" command="SetOnPointFields"
                             number_of_elements="1" default_values="0">
                             <BooleanDomain name="bool"/>
        <Documentation>
          This property indicates if the prefix must be apply to PointFields
        </Documentation>
                             </IntVectorProperty>""")
        def SetOnPointFields(self, val):
            val = int(val)
            if self.onPointFields != val :
                self.onPointFields = val
                self.Modified()

        @smproperty.xml ("""
        <StringVectorProperty name="Prefix"
                             command="SetPrefix"
                             number_of_elements="1"
                             default_values="">
        <Documentation>
          This property indicates if the prefix to be used by the filter
        </Documentation>
        </StringVectorProperty>""")
        def SetPrefix(self,val ):
            self.prefix = str(val)
            self.Modified()

    @smproxy.filter(name="Center DataSet")
    @smhint.xml("""<ShowInMenu category="BasicTools" />""")
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
            input0 = GetInputBasicTools(request, inInfoVec, outInfoVec,FieldsAsTags=True,connection=0)
            input1 = GetInputBasicTools(request, inInfoVec, outInfoVec,FieldsAsTags=False,connection=1)

            mean0 = np.sum(input0.nodes,axis=0)/input0.GetNumberOfNodes()
            mean1 = np.sum(input1.nodes,axis=0)/input1.GetNumberOfNodes()
            # the user must not modify the inputs
            outputmesh = copy.copy(input0)
            outputmesh.nodes = input0.nodes + (mean1 - mean0)
            SetOutputBasicTools(request, inInfoVec, outInfoVec,outputmesh )
            return 1


    @smproxy.filter(name="Mesh Filter")
    @smhint.xml("""<ShowInMenu category="BasicTools" />""")
    @smproperty.input(name="Mesh To Filter", port_index=0)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    class BasicToolsMeshFilter(VTKPythonAlgorithmBase):
        def __init__(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
            self.__dimensionality = None
            self.__elemTypeFilter = {x:0 for x in EN.geoSupport}
            self.__purgeTags = 1
            self.__elemTagsFilter = {}

        def RequestInformation(self, request, inInfoVec, outInfoVec):
            metadata = VtkToMeshOnlyMeta(GetInputVtk(request, inInfoVec, outInfoVec),FieldsAsTags=True)
            tagnames = metadata.elemTags
            self.__elemTagsFilter = {n:self.__elemTagsFilter.get(n,0) for n in tagnames}
            #possible improvement. not ready yet
            #self.__elemTypeFilter = {x:self.__elemTypeFilter.get(n,0) for x in inmesh.elements.keys()}
            return 1

        def RequestData(self, request, inInfoVec, outInfoVec):
            inMesh = GetInputBasicTools(request, inInfoVec, outInfoVec, True)
            from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByElementFilter
            from BasicTools.Containers.UnstructuredMeshFieldOperations import CopyFieldsFromOriginalMeshToTargetMesh
            from BasicTools.Containers.Filters import ElementFilter
            elementTypes = [ k for k,v in self.__elemTypeFilter.items() if v != 0]
            tags = [ k for k,v in self.__elemTagsFilter.items() if v != 0]
            ef  = ElementFilter(mesh=inMesh,dimensionality=self.__dimensionality,elementTypes=elementTypes, tags=tags)

            outMesh = ExtractElementsByElementFilter(inMesh, ef)
            if self.__purgeTags:
                outMesh.nodesTags.RemoveEmptyTags()
                for name,data in outMesh.elements.items():
                    data.tags.RemoveEmptyTags()
            CopyFieldsFromOriginalMeshToTargetMesh(inMesh,outMesh)
            SetOutputBasicTools(request, inInfoVec, outInfoVec, outMesh, TagsAsFields=True)
            return 1

        @smproperty.xml("""<IntVectorProperty
                            name="Dimensionality Filter"
                            command="SetDimensionalityFilter"
                            number_of_elements="1"
                            default_values="-100">
            <EnumerationDomain name="enum">
              <Entry value="-100" text="no dimensionaly filter"/>
              <Entry value="-3" text="Eliminate 3D elements"/>
              <Entry value="-2" text="Eliminate 2D elements"/>
              <Entry value="-1" text="Eliminate 1D elements"/>
              <Entry value="0"  text="Keep only 0D elements"/>
              <Entry value="1"  text="Keep only 1D elements"/>
              <Entry value="2"  text="Keep only 2D elements"/>
              <Entry value="3"  text="Keep only 3D elements"/>
            </EnumerationDomain>
            <Documentation>
              This property indicates which type of elements to keep or filter
            </Documentation>
             </IntVectorProperty>""")
        def SetDimensionalityFilter(self, val):
            val = int(val)
            if val == -100:
                val = None
            elif val <-3 or val > 3:
                raise(Exception("dimensionality must be between -3 and 3"))
            if self.__dimensionality != val :
                self.__dimensionality  = val
                self.Modified()

        @smproperty.xml("""<IntVectorProperty
                             name="PurgeEmptyTags"
                             command="SetPurgeEmptyTags"
                             number_of_elements="1"
                             default_values="1">
                             <BooleanDomain name="bool"/>
            <Documentation>
              This property indicates if tags with zero elements/point must be
              eliminated from the output
            </Documentation>
                             </IntVectorProperty>""")
        def SetPurgeEmptyTags(self, val):
            val = int(val)
            if self.__purgeTags != val :
                self.__purgeTags = val
                self.Modified()

        @smproperty.xml("""<StringVectorProperty information_only="1"
                                name="ElementTypesArrayInfo">
            <ArraySelectionInformationHelper attribute_name="ElementTypes" />
          </StringVectorProperty>
          <StringVectorProperty command="SetElementTypesArrayStatus"
                                element_types="2 0"
                                information_property="ElementTypesArrayInfo"
                                label="Element Types Filter"
                                name="ElementTypesArrayStatus"
                                number_of_elements="0"
                                number_of_elements_per_command="2"
                                repeat_command="1">
            <ArraySelectionDomain name="array_list">
              <RequiredProperties>
                <Property function="ArrayList"
                          name="ElementTypesArrayInfo" />
              </RequiredProperties>
            </ArraySelectionDomain>
            <Documentation>This property lists which ElementTypes are used to filter
            the mesh.</Documentation>
          </StringVectorProperty>""")
        def SetElementTypesArrayStatus(self,key, val):
            if self.__elemTypeFilter[key] != int(val):
                self.__elemTypeFilter[key] = int(val)
                self.Modified()
        def GetNumberOfElementTypesArrays(self):
            return len(EN.geoSupport)

        def GetElementTypesArrayName(self,index):
            return   list(self.__elemTypeFilter.keys())[index]

        def GetElementTypesArrayStatus(self,key):
            return self.__elemTypeFilter[key]

        @smproperty.xml("""
          <StringVectorProperty command="SetElementTagsArrayStatus"
                                element_types="2 0"
                                information_property="ElementTagsArrayInfo"
                                label="Element Tags Filter"
                                name="ElementTagsArrayStatus"
                                number_of_elements="0"
                                number_of_elements_per_command="2"
                                repeat_command="1">
            <ArraySelectionDomain name="array_list">
              <RequiredProperties>
                <Property function="ArrayList"
                          name="ElementTagsArrayInfo" />
              </RequiredProperties>
            </ArraySelectionDomain>
            <Documentation>This property lists which ElementTags are used to filter
            the mesh.</Documentation>
          </StringVectorProperty>""")
        def SetElementTagsArrayStatus(self,key, val):
            if self.__elemTagsFilter[key] != int(val):
                self.__elemTagsFilter[key] = int(val)
                self.Modified()

        @smproperty.xml("""<StringVectorProperty information_only="1"
                                name="ElementTagsArrayInfo">
            <ArraySelectionInformationHelper attribute_name="ElementTags" />
          </StringVectorProperty>""")
        def GetNumberOfElementTagsArrays(self):
            #print("int GetNumberOfElementTagsArrays")
            return len(self.__elemTagsFilter)

        def GetElementTagsArrayName(self,index):
            return list(self.__elemTagsFilter.keys())[index]

        def GetElementTagsArrayStatus(self,key):
            return self.__elemTagsFilter[key]

    PrintDebug("BasicTools ParaView Plugin Loaded")
except:
    print("Error loading BasicTools ParaView Plugin")
    print("BasicTools in the PYTHONPATH ??? ")
    if debug:
        raise
