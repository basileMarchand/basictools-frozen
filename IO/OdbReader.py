# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
import BasicTools.Containers.ElementNames as EN

eltype = {}
eltype["S3"] = EN.Triangle_3
eltype["C3D4"] = EN.Tetrahedron_4
eltype["C3D8"] = EN.Hexaedron_8

permutaion = {}
#permutaion["C3D4"] = [0, 1, 3, 2]

class OdbReader(object):
    def __init__(self):
        super(OdbReader,self).__init__()
        self.canHandleTemporal = True

        self.fieldNameToRead = None
        self.timeToRead = -1
        self.filename =None

        self.time = None
        self.stepData = None
        self.__VariablesToPush = ["fieldNameToRead","filename","timeToRead"]
        self.__VariablesToPull = ["time","stepData"]
        self.__FunctionToBypass = ["SetFileName","SetFieldNameToRead","SetTimeToRead"]
        self.proc = None
        self.client = None
        self.odb = None
        self.output = None

    def __del__(self):
        if not( self.proc is None):
            self.client.Close()
            self.proc = None

    def Reset(self):
        self.fieldNameToRead = None
        self.timeToRead = None
        self.time = None
        self.stepData =None

    def GetAvilableTimes(self):
           return self.time


    def __getattribute__(self,name):
        attr = object.__getattribute__(self, name)
        if hasattr(attr, '__call__'):
            if name in self.__FunctionToBypass:
                return attr
            try:
                import odbAccess as OA
                return attr
            except:
                if self.proc == None:
                    from BasicTools.IO.Wormhole import WormholeServer,WormholeClient,GetPipeWrormholeScript
                    from BasicTools.Helpers.Tests import WriteTempFile
                    import subprocess
                    script = GetPipeWrormholeScript()
                    fn = WriteTempFile("WormholeServer.py",script)
                    self.proc = subprocess.Popen(["abaqus","python",fn], stdout=subprocess.PIPE,stdin=subprocess.PIPE)
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

                    for var in self.__VariablesToPull:
                        print(var)
                        print(object.__getattribute__(self, var))
                        #object.__getattribute__(self, var)  = self.client.RetrieveData("reader.{0}".format(var) )
                    return res

                return newfunc
        else:
            return attr

    def OLD__getattribute__(self,name):
        attr = object.__getattribute__(self, name)
        if hasattr(attr, '__call__'):
            if name in self.__FunctionToBypass:
                return attr
            print(name)
            try:
                import odbAccess as OA
                return attr
            except:
                print("Abaqus API not avilable unsing external call")

                # return a function wrapped in an external call
                from BasicTools.Helpers.Tests import GetUniqueTempFile
                script = """
from BasicTools.Helpers.PrintBypass import PrintBypass
with PrintBypass() as f:
    f.ToDisk("{4}","{5}")

    from BasicTools.IO.OdbReader import OdbReader
    from BasicTools.IO.PipeIO import PipeWriter

    reader = OdbReader()
    reader.SetFileName("{0}")
    if "{1}" == "None":
        reader.fieldNameToRead = None
    else:
        reader.fieldNameToRead = "{1}"
    reader.timeToRead = {2}
    reader.ReadMetaData()
    res = reader.{3}()
    wr = PipeWriter()
    wr.outbuffer = f.stdout_
    wr.Write(res)

""".format(self.filename,self.fieldNameToRead,self.timeToRead,name,GetUniqueTempFile(".log","output_")[1],GetUniqueTempFile(".log2","output_")[1])

                def newfunc(*args, **kwargs):
                    print('before calling %s' %attr.__name__)
                    import os
                    from BasicTools.IO.PipeIO import PipeReader
                    from BasicTools.Helpers.Tests import WriteTempFile
                    import subprocess

                    r, w = os.pipe()
                    os.close(w)

                    fn = WriteTempFile("AbaqusReader.py",script)
                    pr = PipeReader()

                    proc = subprocess.Popen(["abaqus", "python",fn], stdout=subprocess.PIPE)
                    pr.inbuffer = proc.stdout
                    result = pr.Read()
                    proc.wait()
                    #result = attr(*args, **kwargs)
                    print('done calling %s' %attr.__name__)
                    return result

                return newfunc
        else:
            return attr

    def SetFileName(self,fn):
        self.filename = fn

    def pprint(self,vals):
        if self.silent == True:
            print(vals)

    def SetFieldNameToRead(self,val):
        self.fieldNameToRead = val

    def SetTimeToRead(self,time):
        if time is not None:
            self.timeToRead = time

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
                #print(step.totalTime)
        self.time = np.array(self.time)
        return self.time, self.stepData



    def Read(self):
        if self.output == None:
            res = UnstructuredMesh()
            odb = self.Open()
            self.ReadMetaData(odb)

            instance = odb.rootAssembly.instances
            instancename = list(instance.keys())[0]

            nodes = instance[instancename].nodes


            nbnodes = len(nodes)
            res.nodes = np.empty((nbnodes,3),dtype=float)
            res.originalIDNodes = np.empty((nbnodes),dtype=int)
            abatomesh = {}
            cpt = 0
            print("Reading Nodes")
            for i in nodes:
                res.nodes[cpt,:] = i.coordinates
                res.originalIDNodes[cpt] = i.label
                abatomesh[i.label] = cpt
                cpt += 1

            print("Reading Nodes Keys")
            res.PrepareForOutput()
            nSets = instance[instancename].nodeSets
            for nSetK in nSets.keys():
                print(nSetK)
                print("This does not work for the moment")
                nSet = nSets[nSetK]
                name = nSet.name
                tag = res.nodesTags.CreateTag(name,False)
                for node in nSet.nodes:
                    enum = abatomesh[node.label]
                    tag.AddToTag(enum)
                tag.RemoveDoubles()



            elements = instance[instancename].elements
            print("ReadingElements")
            elemMap = {}
            #print(elements.bulkDataBlocks)

            for elem in elements:
                conn = [abatomesh[n] for n in elem.connectivity ]
                elems = res.GetElementsOfType(eltype[elem.type])
                per = permutaion.get(elem.type,None)
                if per is None:
                    enum = elems.AddNewElement(conn,elem.label) - 1
                else:
                    enum = elems.AddNewElement([conn[x] for x in per],elem.label) - 1
                elemMap[elem.label] = (eltype[elem.type],enum)

            print("Reading Keys")
            res.PrepareForOutput()
            eSets = instance[instancename].elementSets
            for eSetK in eSets.keys():
                eSet =eSets[eSetK]
                name = eSet.name
                for elem in eSet.elements:
                    elems = res.GetElementsOfType(eltype[elem.type])
                    enum = elemMap[elem.label][1]
                    elems.GetTag(elem.instanceName).AddToTag(enum)
                    elems.GetTag(name).AddToTag(enum)

            for name,data in res.elements.items():
                for tag in data.tags:
                   tag.RemoveDoubles()

            print("Reading fields")
            res.PrepareForOutput()
            self.abatomesh = abatomesh
            self.elemMap = elemMap
            self.output = res

        self.output.nodeFields,self.output.elemFields = self.ReadFields(self.abatomesh,self.elemMap,self.output)

        return self.output

    def GetActiveFrame(self):

        if self.timeToRead == -1.:
            timeIndex = len(self.time)-1
        else:
            timeIndex = np.argmin(abs(self.time - self.timeToRead ))
        name, val = self.stepData[timeIndex]
        odb = self.Open()
        frame = odb.steps[name].getFrame(val)
        return frame


    def Open(self):
        if not(self.odb is None):
            return self.odb
        print("opening")
        import odbAccess as OA
        print(self.filename)
        self.odb = OA.openOdb(self.filename,readOnly=True)
        print("done")
        self.ReadMetaData()
        return self.odb


    def ReadField(self,fieldname=None,time=None,odb=None):
        if odb is None:
            odb = self.Open()

        if fieldname != None:
            self.SetFieldNameToRead(fieldname)

        if self.fieldNameToRead is None:
            raise(Exception("need a fieldNameToRead to read"))

        if self.timeToRead == -1.:
            timeIndex = len(self.time)-1
        else:
            timeIndex = np.argmin(abs(self.time - self.timeToRead ))
        name, val = self.stepData[timeIndex]

        frame = odb.steps[name].getFrame(val)
        field = frame.fieldOutputs[self.fieldNameToRead]

        # for the moment only nodal is supported
        return np.array(field.bulkDataBlocks[0].data)

    def ReadFields(self,nodeMap,elemMap,res):
        frame = self.GetActiveFrame()
        import odbAccess as OA
        nodalFields = {}
        elemFields = {}
        s1 = 0
        s2 = 1
        for name,data in frame.fieldOutputs.items():
            if self.fieldNameToRead and name != self.fieldNameToRead:
                continue

            if data.type == OA.SCALAR:
              s2 = 1
            elif data.type == OA.VECTOR:
              s2 = 3
            elif data.type == OA.TENSOR_3D_FULL:
              s2 = 6
            else:
                print("do not how to treat {}".format(data.type))
                raise(Exception("error"))

            if data.locations[0].position == OA.CENTROID:
                s1 = len(elemMap)
                storage = elemFields
                storage[name] = self.ReadFieldWithMapElement(elemMap,data,s1,s2,res)

            elif data.locations[0].position == OA.NODAL:
                s1 = len(nodeMap)
                storage = nodalFields
                storage[name] = self.ReadFieldWithMapNode(nodeMap,data,s1,s2)
        return nodalFields, elemFields

    def ReadFieldWithMapNode(self,entityMap,field,s1,s2):
        res = np.zeros((s1,s2))
        fieldValues = field.values
        for v in fieldValues :
            nid = entityMap[v.nodeLabel]
            res[nid,:] = v.data
        return res

    def ReadFieldWithMapElement(self,entityMap,field,s1,s2,mesh):
        res = np.zeros((s1,s2))
        fieldValues = field.values
        for v in fieldValues :
            nid = entityMap[v.elementLabel][1] + mesh.elements[entityMap[v.elementLabel][0]].globaloffset
            res[nid,:] = v.data
        return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".odb",OdbReader)


def CheckIntegrity(GUI = False):
    import time as tt

    at = tt.time()
    reader = OdbReader()
    # no .odb in the database for test for the moment
    return "ok"

    reader.SetFileName("/home/fbordeu-weld/test/odbReader/Job-1.odb")

    time, stepData = reader.ReadMetaData()
    print("time")
    print(time)

    print(stepData)
    reader.timeToRead = 2.0
    reader.SetFieldNameToRead("U")

    print(reader.ReadField())
    print(reader.Read())
    print(tt.time() - at)
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))
