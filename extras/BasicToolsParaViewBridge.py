# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

# this files is inteded to be used inside paraview as a plugin
# compatible with paraview 5.7+

import copy
import time
import locale

from BasicTools.Containers.Filters import ElementFilter
_startTime = time.time()
debug = False
if debug :
    def PrintDebug(mes):
        import time
        print(mes,time.time() - _startTime)
else:
    def PrintDebug(mes):
        pass
try:

    import numpy as np

    from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain, smhint
    from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

    PrintDebug("Loading libs")
    from BasicTools.Bridges.vtkBridge import GetInputVtk, GetOutputVtk, GetInputBasicTools,  SetOutputBasicTools, VtkFieldToNumpyFieldWithDims
    from BasicTools.IO.IOFactory import InitAllReaders
    from BasicTools.IO.IOFactory import ReaderFactory
    from BasicTools.IO.IOFactory import InitAllWriters
    from BasicTools.IO.IOFactory import WriterFactory
    from BasicTools.Bridges.vtkBridge import VtkToMesh, VtkToMeshOnlyMeta
    import BasicTools.Containers.ElementNames as EN
    PrintDebug("Loading")

    paraview_plugin_name = "BasicTools ParaView Bridge"
    paraview_plugin_version = "5.7"

    #----------------------------- The Readers ------------------------------------
    PrintDebug("Init readers ")
    InitAllReaders()
    PrintDebug("Init readers Done")
    try :
        PrintDebug("loading meshio readers")

        import BasicTools.Bridges.MeshIOBridge as MeshIOBridge
        MeshIOBridge.InitAllReaders()
        MeshIOBridge.AddReadersToBasicToolsFactory()
        PrintDebug("loading meshio readers Done")
    except:
        print("Error in the registration of meshio readers")
        if debug:
            raise


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

        @smproperty.xml("""<IntVectorProperty
                            name="Encoding"
                            command="SetEncoding"
                            panel_visibility="advanced"
                            number_of_elements="1"
                            default_values="-1">
            <EnumerationDomain name="enum">
              <Entry value="-1" text="locale"/>
              <Entry value="0" text="utf-8"/>
              <Entry value="1" text="ascii"/>
              <Entry value="2" text="latin-1"/>
            </EnumerationDomain>
            <Documentation>
              This property indicates with which encoding to use to open the file
            </Documentation>
             </IntVectorProperty>""")
        def SetEncoding(self, val):
            reader =  self.basicToolsReader
            val = int(val)
            data = {-1:locale.getpreferredencoding(False),
                     0:"utf-8",
                     1:"ascii",
                     2:"latin-1"}
            val = data[val]

            if reader.encoding == val:
                return
            reader.encoding = val
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
                   "SetEncoding":SetEncoding,
                   "RequestData":RequestData,
                   "RequestInformation":RequestInformation}
                   )

        name = " files [BasicTools Reader]"
        extraname = ""
        if readerClassName.find("MeshIO") >=0 :
            extraname = "MeshIO "
            name = " files [BasicTools "+extraname+"Reader]"

        obj2 = smproxy.reader(name=wrapperClassName, label="BasicTools "+extraname+"Python-based "+ ext.upper() +" Reader",
                    extensions=ext, file_description=ext+name)(obj)

        locals()[wrapperClassName] = obj2


    #----------------------------- The Writers ------------------------------------


    PrintDebug("InitAllWriters")
    InitAllWriters()
    PrintDebug('InitAllWriters Done')

    try :
        PrintDebug("loading meshio  Done")
        import BasicTools.Bridges.MeshIOBridge as MeshIOBridge
        MeshIOBridge.InitAllWriters()
        MeshIOBridge.AddWritersToBasicToolsFactory()
        PrintDebug("loading meshio writers Done")
    except:
        print("Error in the registration of meshio writers")
        if debug:
            raise

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
            mesh = GetInputBasicTools(request,inInfoVec,outInfoVec,FieldsAsTags=True,connection=0,port=0)

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

    def ConstructSkeletonForTensorVectorFields(data, sep="_",sufix = [{"x":0,"y":1,"z":2},{"X":0,"Y":1,"Z":2},{"XX":0, "YY":1, "ZZ":2, "XY":3, "XZ":4, "YZ":5},{"xx":0, "yy":1, "zz":2, "xy":3, "xz":4, "yz":5}]):
        from collections import defaultdict
        res = defaultdict(lambda : dict())
        for name in data:
            done = False
            for family in sufix:
                for i,(s,pos) in enumerate(family.items()):
                    ss = sep+s
                    index= name.find(ss)
                    if index > -1 and index == len(name)-len(ss):
                        newname = name[0:len(name)-len(ss)]
                        if pos in res[newname]:
                            raise(Exception(f"field {name} confict with {res[newname][pos]}"))
                        res[newname][pos] = name

                        done = True
            if not done:
                if 0 in res[name]:
                    raise(Exception(f"field {name} confict with {res[newname][0]}"))
                res[name][0] = name
        return res


    @smproxy.filter(name="Scalar To Tensor/Vector")
    @smhint.xml("""<ShowInMenu category="BasicTools" />""")
    @smproperty.input(name="Input", port_index=0)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    class ScalarToTensorVector(VTKPythonAlgorithmBase):
        def __init__(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
            self.onCellData =  True
            self.onPointData =  True

        def RequestData(self, request, inInfoVec, outInfoVec):

            def TreatDataSet(inputData,outputData,nbtuples):
                nbarrays = outputData.GetNumberOfArrays()
                arraynames = []
                for x in range(nbarrays):
                    array = inputData.GetArray(x)
                    if array.GetNumberOfComponents() == 1:
                        arraynames.append(array.GetName() )

                res = ConstructSkeletonForTensorVectorFields(arraynames)

                for name, data in res.items():
                    nbcomponents = max(data.keys())+1
                    npdata = np.zeros((nbtuples,nbcomponents), dtype=float)
                    for index, fn in data.items():
                        npdata[:,index] = VtkFieldToNumpyFieldWithDims(inputData.GetArray(fn) )[1]
                        outputData.RemoveArray(fn)
                    newdata = numpy_support.numpy_to_vtk(num_array=npdata, deep=True)
                    newdata.SetNumberOfComponents(nbcomponents)
                    newdata.SetName(name)
                    if nbcomponents == 3:
                        for i in range(3):
                            newdata.SetComponentName(i,["X","Y","Z"][i])
                    elif nbcomponents == 6:
                        for i in range(6):
                            newdata.SetComponentName(i,["XX","YY","ZZ","XY","XZ","YZ"][i])
                    outputData.AddArray(newdata)

            from vtkmodules.util import numpy_support
            inputMesh = GetInputVtk(request, inInfoVec, outInfoVec)
            outputMesh = GetOutputVtk(request, inInfoVec, outInfoVec, copyAttr = False, outputNumber= 0)
            outputMesh.ShallowCopy(inputMesh)
            from vtkmodules.vtkCommonDataModel import vtkPointData

            if self.onCellData:
                TreatDataSet(inputMesh.GetCellData(), outputMesh.GetCellData(), inputMesh.GetNumberOfCells())

            if self.onPointData:
                TreatDataSet(inputMesh.GetPointData(), outputMesh.GetPointData(), inputMesh.GetNumberOfPoints())

            return 1

        @smproperty.xml("""<IntVectorProperty name="OnCellData" command="SetOnCellData"
                             number_of_elements="1" default_values="1">
                             <BooleanDomain name="bool"/>
        <Documentation>
          This property indicates if the cell data must be treated
        </Documentation>
                             </IntVectorProperty>""")
        def SetOnCellData(self, val):
            val = int(val)
            if self.onCellData != val :
                self.onCellData = val
                self.Modified()

        @smproperty.xml("""<IntVectorProperty name="OnPointData" command="SetOnPointData"
                             number_of_elements="1" default_values="1">
                             <BooleanDomain name="bool"/>
        <Documentation>
          This property indicates if the point data must be treated
        </Documentation>
                             </IntVectorProperty>""")
        def SetOnPointData(self, val):
            val = int(val)
            if self.onPointData != val :
                self.onPointData = val
                self.Modified()

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

    @smproxy.filter(name="Transfer Data")
    @smhint.xml("""<ShowInMenu category="BasicTools" />""")
    @smproperty.input(name="Target Points", port_index=1)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    @smproperty.input(name="Source of Data", port_index=0)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    class TransferData(VTKPythonAlgorithmBase):
        def __init__(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1, outputType="vtkUnstructuredGrid")
            self.__method = 2

        def FillInputPortInformation(self, port, info):
            if port == 0:
                info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
                info.Append(self.INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData")
            else:
                info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
                info.Append(self.INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData")
            return 1

        def RequestData(self, request, inInfoVec, outInfoVec):
            possibleMethods =["Interp/Nearest","Nearest/Nearest","Interp/Clamp","Interp/Extrap","Interp/ZeroFill"]

            from vtkmodules.vtkCommonDataModel import  vtkDataSet
            input0 = GetInputBasicTools(request, inInfoVec, outInfoVec,FieldsAsTags=True,connection=0)
            input1 = GetInputBasicTools(request, inInfoVec, outInfoVec,FieldsAsTags=False,connection=1)
            input0.ConvertDataForNativeTreatment()
            input1.ConvertDataForNativeTreatment()

            from BasicTools.FE.Fields.FEField import FEField
            from BasicTools.FE.FETools import PrepareFEComputation
            from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
            space, numberings, offset, NGauss = PrepareFEComputation(input0,numberOfComponents=1)
            field = FEField("",mesh=input0,space=space,numbering=numberings[0])

            op,status = GetFieldTransferOp(field,input1.nodes,method=possibleMethods[self.__method], verbose= False,elementFilter=ElementFilter(input0))

            # the user must not modify the inputs
            outputmesh = copy.copy(input1)
            outputmesh.nodeFields["status"] = status
            for n,v in input0.nodeFields.items():
                outputmesh.nodeFields[n] = op.dot(v)
            outputmesh.nodeFields["SourceCoords"] = op.dot(input0.nodes)
            SetOutputBasicTools(request, inInfoVec, outInfoVec,outputmesh )
            return 1

        @smproperty.xml("""<IntVectorProperty
                            name="Method inside/outside"
                            command="SetMethod"
                            number_of_elements="1"
                            default_values="2">
            <EnumerationDomain name="enum">
              <Entry value="0" text="Interp/Nearest"/>
              <Entry value="1" text="Nearest/Nearest"/>
              <Entry value="2" text="Interp/Clamp"/>
              <Entry value="3" text="Interp/Extrap"/>
              <Entry value="4"  text="Interp/ZeroFill"/>

            </EnumerationDomain>
            <Documentation>
              This property indicates the type of method to use when a point is inside/outside of the Source Data
            </Documentation>
             </IntVectorProperty>""")
        def SetMethod(self, val):
            val = int(val)
            if val < 0 or val > 4:
                raise(Exception("Method must be between 0 and 4"))
            if self.__method != val :
                self.__method  = val
                self.Modified()

    @smproxy.filter(name="Tensor Rotation")
    @smhint.xml("""<ShowInMenu category="BasicTools" />""")
    @smproperty.input(name="Input", port_index=0)
    @smdomain.datatype(dataTypes=["vtkUnstructuredGrid","vtkPolyData"], composite_data_supported=False)
    class TensorRotation(VTKPythonAlgorithmBase):
        def __init__(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
            self.V1Name = None
            self.V2Name = None
            self.__firstStressArray = None
            self.WorkOnCells = False
            self.__inverse = False

        def RequestData(self, request, inInfoVec, outInfoVec):
            from vtkmodules.vtkCommonDataModel import  vtkDataSet
            input0 = GetInputBasicTools(request, inInfoVec, outInfoVec,FieldsAsTags=True,connection=0)

            mean0 = np.sum(input0.nodes,axis=0)/input0.GetNumberOfNodes()
            #mean1 = np.sum(input1.nodes,axis=0)/input1.GetNumberOfNodes()
            # the user must not modify the inputs
            outputmesh = copy.copy(input0)
            if self.WorkOnCells:
                data = copy.copy(input0.elemFields)
            else:
                data = copy.copy(input0.nodeFields)
            from BasicTools.Containers.UnstructuredMeshFieldOperations import  ApplyRotationMatrixTensorField

            db = {}
            db["XX"] =[["XX","XY","XZ"],["XY","YY","YZ"],["XZ","YZ","ZZ"] ]
            db["xx"] =[["xx","xy","xz"],["xy","yy","yz"],["xz","yz","zz"] ]
            db["00"] =[["00","01","02"],["01","11","12"],["02","12","22"] ]
            db["11"] =[["11","12","13"],["12","22","23"],["13","23","33"] ]

            part1 = self.__firstStressArray[:-2]
            part2 = db[self.__firstStressArray[-2:]]
            fieldstoTreat = [ [ part1+x for x in p] for p in part2]

            newfields = ApplyRotationMatrixTensorField(data, fieldstoTreat, baseNames=[self.V1Name,self.V2Name],inplace=False,prefix="new_",inverse=self.__inverse)

            if self.WorkOnCells:
                data.update(newfields)
                outputmesh.elemFields = data
            else:
                data = input0.nodeFields
                outputmesh.elemFields = data

            SetOutputBasicTools(request, inInfoVec, outInfoVec, outputmesh, TagsAsFields=True )
            return 1

        @smproperty.xml("""<StringVectorProperty command="SetV1Array"
                            default_values="1"
                            element_types="0 0 0 0 2"
                            name="OrientationV1Array"
                            number_of_elements="5">
        <ArrayListDomain attribute_type="Vectors"
                         input_domain_name="vector_array"
                         name="array_list"
                         none_string="V1 not selected">
          <RequiredProperties>
            <Property function="Input"
                      name="Input" />
          </RequiredProperties>
        </ArrayListDomain>
        <Documentation>
Select the input array to use for orienting the glyphs.
        </Documentation>
      </StringVectorProperty>""")
        def SetV1Array(self, *val):
            #PointData            vals  (1, 0, 0, 0, 'Normals')
            #CellData             vals  (1, 0, 0, 1, 'Normals')

            if self.V1Name != val[4] or self.WorkOnCells != val[3]:
                self.V1Name = val[4]
                self.WorkOnCells = val[3]
                self.Modified()

        @smproperty.xml("""<StringVectorProperty command="SetV2Array"
                            default_values="1"
                            element_types="0 0 0 0 2"
                            name="OrientationV2Array"
                            number_of_elements="5">
        <ArrayListDomain attribute_type="Vectors"
                         input_domain_name="vector_array"
                         name="array_list"
                         none_string="V2 not selected">
          <RequiredProperties>
            <Property function="Input"
                      name="Input" />
          </RequiredProperties>
        </ArrayListDomain>
        <Documentation>
Select the input array to use for orienting the glyphs.
        </Documentation>
      </StringVectorProperty>""")
        def SetV2Array(self, *val):
            if self.V2Name != val[4] or self.WorkOnCells != val[3] :
                self.V2Name = val[4]
                self.WorkOnCells = val[3]
                self.Modified()

        @smproperty.xml("""<StringVectorProperty command="SetStressArray"
                            default_values="1"
                            element_types="0 0 0 0 2"
                            name="FirstStressArray"
                            number_of_elements="5">
        <ArrayListDomain attribute_type="Vectors"
                         input_domain_name="vector_array"
                         name="array_list"
                         none_string="V2 not selected">
          <RequiredProperties>
            <Property function="Input"
                      name="Input" />
          </RequiredProperties>
        </ArrayListDomain>
        <Documentation>
Select the input array to use for orienting the glyphs.
        </Documentation>
      </StringVectorProperty>""")
        def SetStressArray(self, *val):
            if self.__firstStressArray != val[4] or self.WorkOnCells != val[3] :
                self.__firstStressArray = val[4]
                self.WorkOnCells = val[3]
                self.Modified()

        @smproperty.xml("""<IntVectorProperty
                             name="Inverse"
                             command="SetInverse"
                             number_of_elements="1"
                             default_values="0">
                             <BooleanDomain name="bool"/>
            <Documentation>
              This property indicates if we need to apply the invers of the transformation.
            </Documentation>
            </IntVectorProperty>""")
        def SetInverse(self, val):
            val = int(val)
            if self.__inverse != val :
                self.__inverse = val
                self.Modified()

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
            if (key in self.__elemTagsFilter) and (self.__elemTagsFilter[key] != int(val)) :
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


###Generic algorithm interface
# This interface is intended to facilitate the interface of single function BasicTools algorithm.
# this means a function with the signature :
#    Func(UnstructuredMesh inputA, UnstructuredMesh inputB, ...) -> (UnstructuredMesh outputA, ...)

from BasicTools.Bridges.vtkBridge import MeshToVtk

def WrapBasicToolsFunctionToVTK(function, inputs, outputs, options, description=""):
    if isinstance(inputs,str):
        inputs = (inputs,)

    if isinstance(outputs,str):
        outputs = (outputs,)

    def GetInit(function, nInputPorts=1, nOutputPorts=1, ops={}):

        def myFunctionAsFilterInit(self):
            VTKPythonAlgorithmBase.__init__(self, nInputPorts=nInputPorts, nOutputPorts=nOutputPorts, outputType='vtkUnstructuredGrid')
            self.function = function
            self.nInputPorts = nInputPorts
            self.nOutputPorts = nOutputPorts
            self.__TagsAsFields = True
            self.ops = dict(ops)

        return myFunctionAsFilterInit

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

    def RequestData(self, request, inInfoVec, outInfoVec):

        inputmesh = []
        for i in range(self.nInputPorts):
            mesh = GetInputBasicTools(request, inInfoVec, outInfoVec,FieldsAsTags=True, connection=0, port=i)
            mesh.ConvertDataForNativeTreatment()
            inputmesh.append(mesh)


        meshs = self.function(*inputmesh,**self.ops)

        if self.nOutputPorts > 1:
            for i in range(self.nOutputPorts):
                vtkoutputi = GetOutputVtk(request,inInfoVec, outInfoVec,False,i)
                MeshToVtk(meshs[i], vtkoutputi,TagsAsFields=self.__TagsAsFields)
        else:
            SetOutputBasicTools(request, inInfoVec, outInfoVec, meshs,TagsAsFields=self.__TagsAsFields)

        return 1

    data = {}
    ops = {}
    def get_default_args(func):
        import inspect
        signature = inspect.signature(func)
        return {
            k: v.default
            for k, v in signature.parameters.items()
            if v.default is not inspect.Parameter.empty
        }

    signature = get_default_args(function)
    for name, otype in options:

        default = signature.get(name,otype())
        if otype == int:
            @smproperty.xml (f"""
            <IntVectorProperty name="Set {name}"
                                command="Set{name}"
                                number_of_elements="1"
                                default_values="{default}">
            <Documentation>
                This is an automatic property setter
            </Documentation>
            </IntVectorProperty>""")
            def SetParam(self, val):
                if val != self.ops[name]:
                    self.ops[name] = val
                    self.Modified()
        elif otype == bool:
            @smproperty.xml (f"""
            <IntVectorProperty name="Set {name}"
                                command="Set{name}"
                                number_of_elements="1"
                                default_values="{default}">
                <BooleanDomain name="bool"/>
            <Documentation>
                This is an automatic property setter
            </Documentation>
            </IntVectorProperty>""")
            def SetParam(self, val):
                if val != self.ops[name]:
                    self.ops[name] = val
                    self.Modified()

        data[f"Set{name}"] = SetParam
        ops[name] = default

    wrapperClassName = "BasicToolsWrapped" + function.__name__
    data.update({"__init__":GetInit(function, len(inputs), len(outputs),ops ),
                "RequestData":RequestData,
                "SetTagsAsFields":SetTagsAsFields})
    obj = type(wrapperClassName,
               (VTKPythonAlgorithmBase,),
                data
    )

    for i,inputtype in enumerate(inputs):
        obj2_tmp = smdomain.datatype(dataTypes=[inputtype], composite_data_supported=False)(obj)
        obj = smproperty.input(name=f"Input Mesh {i}", port_index=i)(obj2_tmp)
    obj2 = smhint.xml("""<ShowInMenu category="BasicTools" />""")(obj)
    obj3 = smproxy.filter(name=description)(obj2)

    PrintDebug(f"{wrapperClassName} registered")
    obj3.wrappedNameBasicTools = wrapperClassName
    return wrapperClassName, obj3

################ Stl to mesh ##########################
from BasicTools.Bridges.gmshBridge import StlToMesh
name,cl = WrapBasicToolsFunctionToVTK( StlToMesh,("vtkPolyData",),("vtkUnstructuredGrid",),(), description= "2D Mesh To 3D Mesh (gmsh)" )
locals()[name] = cl

from BasicTools.Containers.UnstructuredMeshCreationTools import SubDivideMesh
name,cl = WrapBasicToolsFunctionToVTK( SubDivideMesh,("vtkUnstructuredGrid",),("vtkUnstructuredGrid",),( ("level",int),  ), description= "SubDivide a Mesh " )
locals()[cl.wrappedNameBasicTools] = cl

