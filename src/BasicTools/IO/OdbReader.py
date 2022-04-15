# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import os
import numpy as np
from collections import defaultdict

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.IO.AbaqusTools import InpNameToBasicTools, permutation

abaqus_EXEC = os.environ.get("ABAQUS_EXEC","abaqus")

class OdbReader(BaseOutputObject):
    def __init__(self):
        super(OdbReader,self).__init__()
        self.canHandleTemporal = True

        self.fieldsNamesToRead = []
        self.intanceNumberToRead = 0
        self.timeToRead = -1
        self.filename =None

        self.time = None
        self.stepData = None
        self.__VariablesToPush = ["fieldsNamesToRead","filename","timeToRead"]
        self.__FunctionToBypass = ["SetFileName","SetFieldNameToRead","SetTimeToRead"]
        self.proc = None
        self.client = None
        self.odb = None
        self.output = None
        self.encoding = None



    def Reset(self):
        self.fieldsNamesToRead = []
        self.intanceNumberToRead = 0
        self.timeToRead = None
        self.time = None
        self.stepData =None

    def GetAvailableTimes(self):
           return self.time

    def __getattribute__(self, name):
        attr = object.__getattribute__(self, name)
        if hasattr(attr, '__call__'):
            if name in self.__FunctionToBypass:
                return attr
            try:
                import odbAccess as OA
                return attr
            except:
                if self.proc == None:
                    from BasicTools.IO.Wormhole import WormholeClient,GetPipeWormholeScript
                    from BasicTools.Helpers.Tests import WriteTempFile
                    import subprocess
                    script = GetPipeWormholeScript()
                    fn = WriteTempFile("WormholeServer.py",script)
                    self.proc = subprocess.Popen([abaqus_EXEC,"python","-B",fn], stdout=subprocess.PIPE,stdin=subprocess.PIPE)
                    self.client = WormholeClient(proc=self.proc)
                    self.client.RemoteExec("from BasicTools.IO.OdbReader import OdbReader")
                    self.client.RemoteExec("reader = OdbReader()")

                def newfunc(*args, **kwargs):
                    for var in self.__VariablesToPush:
                        val = object.__getattribute__(self, var)
                        self.client.SendData(var,object.__getattribute__(self, var))
                        if type(val) == str:
                            self.client.RemoteExec("{0} = str({0})".format(var) )

                        self.client.RemoteExec("reader.{0} = {0}".format(var) )
                    self.client.SendData("args",args)
                    self.client.SendData("kwargs",kwargs)
                    self.client.RemoteExec("res = reader.{0}(*args, **kwargs)".format(name) )
                    res = self.client.RetrieveData("res")
                    return res

                return newfunc
        else:
            return attr

    def SetFileName(self, fn):
        self.filename = fn

    def SetIntanceNumberToRead(self, i):
        self.intanceNumberToRead = i

    def SetFieldsNamesToRead(self, val):
        self.fieldsNamesToRead = val

    def SetTimeToRead(self, time, timeIndex=None):
        if time is not None:
            self.timeToRead = time
        elif timeIndex is not None:
            self.timeToRead = self.time[timeIndex]

    def ReadMetaData(self,odb = None):
        if not(self.stepData is None):
            return self.time, self.stepData

        if odb is None:
            odb = self.Open()
        self.stepData = []
        cpt = 0
        time =0
        self.time = []
        # loop over the steps
        for k,step in odb.steps.items():
            fcpt = -1
            for frame in step.frames:
                fcpt += 1
                if cpt != 0:
                    # if the frameValue is 0 then the frame is droped because
                    # is the same as the last frame (last frame of the previous
                    # step)
                    if frame.frameValue == 0:
                        continue
                    time += frame.frameValue
                cpt += 1
                self.time.append(float(time))
                self.stepData.append( (str(k),int(fcpt) ))
        self.time = np.array(self.time)
        return self.time, self.stepData

    def ConvertInstanceToBasicTools(self,instance):
        res = UnstructuredMesh()
        nodes = instance.nodes

        nbnodes = len(nodes)
        res.nodes = np.empty((nbnodes,3),dtype=float)
        res.originalIDNodes = np.empty((nbnodes),dtype=int)
        abaToMeshNode = {}
        cpt = 0
        self.PrintDebug("Reading Nodes")
        for i in nodes:
            res.nodes[cpt,:] = i.coordinates
            res.originalIDNodes[cpt] = i.label
            abaToMeshNode[i.label] = cpt
            cpt += 1

        self.PrintDebug("Reading Nodes Keys")
        res.PrepareForOutput()
        nSets = instance.nodeSets
        for nSetK in nSets.keys():
            nSet = nSets[nSetK]
            name = nSet.name
            tag = res.nodesTags.CreateTag(name,False)
            buffer = []
            for node in nSet.nodes:
                buffer.append(node.label)
            enum = [abaToMeshNode[x] for x in buffer ]
            tag.AddToTag(enum)
            tag.RemoveDoubles()

        elements = instance.elements
        self.PrintDebug("Reading Elements")
        elemToMeshElem = {}

        for elem in elements:
            eltype = InpNameToBasicTools[elem.type]
            conn = [abaToMeshNode[n] for n in elem.connectivity ]
            elems = res.GetElementsOfType(eltype)
            per = permutation.get(elem.type,None)
            enum = elems.AddNewElement(conn,elem.label)-1
            elemToMeshElem[elem.label] = (eltype,enum)

        for name,data in res.elements.items():
            per = permutation.get(name,None)
            if per is not None :
                data.connectivity = data.connectivity[:,per]

        self.PrintDebug("Reading Elements Keys")
        res.PrepareForOutput()
        eSets = instance.elementSets
        for eSetK in eSets.keys():
            eSet =eSets[eSetK]
            self.PrintDebug("Reading set {eSetK}".format(eSetK=eSetK))
            name = eSet.name
            cpt = 0
            temptagdata = defaultdict(lambda : list() )
            for elem in eSet.elements:
                temptagdata[elem.type].append(elem.label)

            for elemstype,label in  temptagdata.items():
                res.GetElementsOfType(InpNameToBasicTools[elemstype]).GetTag(name).SetIds( [ elemToMeshElem[x][1] for x in label ] )

        for name,data in res.elements.items():
            for tag in data.tags:
               tag.RemoveDoubles()

        res.PrepareForOutput()
        return res, abaToMeshNode, elemToMeshElem

    def Read(self):
        if self.output == None:

            odb = self.Open()
            self.ReadMetaData(odb)

            instance = odb.rootAssembly.instances
            instancename = list(instance.keys())[self.intanceNumberToRead]
            self.__currentInstance = instance[instancename]
            res, abaToMeshNode, elemToMeshElem = self.ConvertInstanceToBasicTools(instance[instancename])

            self.abatomesh = abaToMeshNode
            self.elemMap = elemToMeshElem
            self.output = res

        self.PrintDebug("Reading Fields")
        self.output.nodeFields,self.output.elemFields = self.ReadFields(self.abatomesh,self.elemMap,self.output)

        return self.output

    def GetActiveFrame(self):
        if self.timeToRead == -1.:
            timeIndex = len(self.time)-1
        else:
            timeIndex = np.argmin(abs(self.time - self.timeToRead ))
        name, val = self.stepData[timeIndex]
        odb = self.Open()
        frame = odb.steps[name].frames[val]
        return frame

    def Open(self):
        if not(self.odb is None):
            return self.odb
        import odbAccess as OA
        self.odb = OA.openOdb(self.filename,readOnly=True)
        self.ReadMetaData()
        return self.odb

    def ReadFields(self,nodeMap,elemMap,res):
        frame = self.GetActiveFrame()
        import odbAccess as OA
        nodalFields = {}
        elemFields = {}
        s1 = 0
        s2 = 1
        for name,data in frame.fieldOutputs.items():

            if len(self.fieldsNamesToRead) != 0  and name not in self.fieldsNamesToRead:
                continue

            if len(self.fieldsNamesToRead) == 0:
                if name in [ "COORD" ]:
                    self.PrintDebug("skip field {}".format(name))
                    continue

            self.PrintDebug("Reading {name} (type:{type})".format(name=name,type=str(data.type)))
            if data.type == OA.SCALAR:
              s2 = 1
            elif data.type == OA.VECTOR:
              s2 = 3
            elif data.type == OA.TENSOR_3D_FULL:
              s2 = 6
            else:
                self.PrintDebug("do not how to treat {}".format(data.type))
                raise(Exception("error"))

            if data.locations[0].position == OA.CENTROID:
                s1 = len(elemMap)
                storage = elemFields
                storage[name] = self.ReadFieldWithMapElement(elemMap,data,s1,s2,res)

            elif data.locations[0].position == OA.NODAL:
                s1 = len(nodeMap)
                storage = nodalFields
                storage[name] = self.ReadFieldWithMapNode(nodeMap,data,s1,s2)
            elif data.locations[0].position == OA.INTEGRATION_POINT:
                sdata = data.getSubset(position=OA.CENTROID)
                s1 = len(elemMap)
                storage = elemFields
                storage[name] = self.ReadFieldWithMapElement(elemMap,data,s1,s2,res)
            else:
                sdata = data.getSubset(position=OA.NODAL)
                s1 = len(nodeMap)
                storage = nodalFields
                storage[name] = self.ReadFieldWithMapNode(nodeMap,sdata,s1,s2)
        self.PrintDebug("Done Reading fields")
        return nodalFields, elemFields

    def ReadFieldWithMapNode(self,entityMap,field,s1,s2):
        res = np.zeros((s1,s2))
        fieldValues = field.values
        for v in fieldValues :
            if self.__currentInstance != v.instance:
                continue
            try:
                nid = entityMap[v.nodeLabel]
            except:
                continue
            res[nid,:] = v.data
        return res

    def ReadFieldWithMapElement(self,entityMap,field,s1,s2,mesh):
        res = np.zeros((s1,s2))
        fieldValues = field.values
        for v in fieldValues :
            if self.__currentInstance != v.instance:
                continue
            eltype, localnumb = entityMap[v.elementLabel]
            nid = localnumb + mesh.elements[eltype].globaloffset
            res[nid,:] = v.data
        return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".odb",OdbReader)


def CheckIntegrity(GUI = False):
    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("ABAQUS_NO_FAIL"): return "ok"

    import time as tt

    at = tt.time()
    reader = OdbReader()
    # no .odb in the database for test for the moment
    return "ok"

    reader.SetFileName("path/Job-1.odb")

    time, stepData = reader.ReadMetaData()
    print("time")
    print(time)

    print(stepData)
    reader.timeToRead = 2.0
    reader.SetFieldsNamesToRead(["U"])

    print(reader.Read())
    print(tt.time() - at)
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))
