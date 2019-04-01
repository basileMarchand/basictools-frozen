# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
import BasicTools.Containers.ElementNames as EN

eltype = {}
eltype["S3"] = EN.Triangle_3
eltype["C3D4"] = EN.Tetrahedron_4

permutaion = {}
#permutaion["C3D4"] = [0, 1, 3, 2]

class OdbReader(object):
    def __init__(self):
        super(OdbReader,self).__init__()
        self.fieldNameToRead = None
        self.timeToRead = -1
        self.time = None
        self.stepData = None
        self.filename =None
        self.canHandleTemporal = True

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
            if name in ["SetFileName","SetFieldNameToRead","SetTimeToRead"]:
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

        if odb is None:
            odb = self.Open()
        print("coucou")
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
        for i in nodes:
            res.nodes[cpt,:] = i.coordinates
            res.originalIDNodes[cpt] = i.label
            abatomesh[i.label] = cpt
            cpt += 1

        elements = instance[instancename].elements

        elemMap = {}
        for elem in elements:
            #print(elem.label)
            conn = [abatomesh[n] for n in elem.connectivity ]
            elems = res.GetElementsOfType(eltype[elem.type])
            per = permutaion.get(elem.type,None)
            if per is None:
                enum = elems.AddNewElement(conn,elem.label) - 1
            else:
                enum = elems.AddNewElement([conn[x] for x in per],elem.label) - 1
            elemMap[elem.label] = (eltype[elem.type],enum)

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

        if self.timeToRead == -1.:
            timeIndex = len(self.timeSteps)-1
        else:
            timeIndex = np.argmin(abs(self.time - self.timeToRead ))

        name, val = self.stepData[timeIndex]
        frame = odb.steps[name].getFrame(val)

        res.PrepareForOutput()
        res.nodeFields,res.elemFields = self.ReadFields(frame,abatomesh,elemMap,res)

        return res

    def Open(self):
        print("opening")
        import odbAccess as OA
        print(self.filename)
        odb = OA.openOdb(self.filename,readOnly=True)
        print("done")
        self.ReadMetaData(odb)
        return odb


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

    def ReadFields(self,frame,nodeMap,elemMap,res):
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
    reader = OdbReader()
    # no .odb in the database for test for the moment
    return "OK"
    reader.SetFileName("/home/fbordeu-weld/test/odbReader/Job-1.odb")

    time, stepData = reader.ReadMetaData()
    print(time)
    print(stepData)
    reader.timeToRead = 2.0
    reader.SetFieldNameToRead("U")

    print(reader.ReadField())
    print(reader.Read())

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))