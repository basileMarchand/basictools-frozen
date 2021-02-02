# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Containers.Filters import ElementFilter
from BasicTools.FE.KR.KRBase import KRBaseVector
from BasicTools.Containers.Octree import Octree



class KRConformalTieVector(KRBaseVector):
    def __init__(self):
        super(KRConformalTieVector,self).__init__()
        self.tol = 0.0001
        self.onII =[]
        self.argsII = []

    def AddArgII(self,name):
        self.argsII.append(name)
        return self

    def From(self,offset=None,first=None,second=None):
        if offset is not None:
            self.originSystem.SetOffset(offset)
        if first is not None:
            self.originSystem.SetFirst(first)
        if second is not None:
            self.originSystem.SetSecond(second)
        return self

    def To(self,offset=None,first=None,second=None):
        if offset is not None:
            self.targetSystem.SetOffset(offset)
        if first is not None:
            self.targetSystem.SetFirst(first)
        if second is not None:
            self.targetSystem.SetSecond(second)
        return self

    def SideI(self,zone):
        return self.On(zone)

    def SideII(self,zone):
        if type(zone) is list:
            self.onII.extend(zone)
        else:
            self.onII.append(zone)
        self.onII = list(set(self.onII))
        return self


    def GenerateEquations(self,meshI=None,fields=None,CH=None,meshII=None,fieldsII=None):

        CH = self._GetConstraintHolder(CH)

        if meshII is None:
            meshII = meshI
            fieldsII = fields

        if len(self.argsII) == 0:
            argsII = self.args
        else:
            argsII = self.argsII

        #submesh1 = ExtractElementByTags(meshI,self.on,cleanLonelyNodes=False)
        #submesh2 = ExtractElementByTags(meshII,self.onII,cleanLonelyNodes=False)

        meshI.ComputeBoundingBox()
        meshII.ComputeBoundingBox()

        bmin1 = meshI.boundingMin
        bmax1 = meshI.boundingMax
        bmin2 = meshII.boundingMin
        bmax2 = meshII.boundingMax

        fieldDicI = {f.name:f for f in fields }
        fieldDicII = {f.name:f for f in fields }
        #self.__fieldOffsetsI = { }
        #self.__fieldOffsetsII = { }

        offsets  , fieldOffsetsI  = self._ComputeOffsets(fields)
        offsetsII, fieldOffsetsII = self._ComputeOffsets(fieldsII)

        oc1 = Octree(bmax1[0],bmax1[1],bmax1[2],bmin1[0],bmin1[1],bmin1[2])
        oc2 = Octree(bmax2[0],bmax2[1],bmax2[2],bmin2[0],bmin2[1],bmin2[2])


        def FillOctree(mesh,oc, tags, base ):
            usedNodes = set()

            ef = ElementFilter(mesh,tags=tags)
            for name,data, ids in ef:
                nids = data.GetNodesIdFor(ids)
                usedNodes.update(nids)

            usedNodes = list(usedNodes)
            nodes = base.ApplyTransform(mesh.nodes[usedNodes,:])
            for cpt,nid in enumerate(usedNodes):
                oc.add_item(nid,nodes[cpt])

            return (usedNodes, nodes)

        usedNodesMeshI, nodesI  = FillOctree(meshI,oc1, self.on, self.originSystem.GetOrthoNormalBase())
        usedNodesMeshII, nodesII = FillOctree(meshII,oc2,self.onII, self.targetSystem.GetOrthoNormalBase())


        #print("usedNodesMeshI")
        #print(usedNodesMeshI)
        #print("usedNodesMeshII")
        #print(usedNodesMeshII)

        ffI = []
        usedOffsetsI = []
        ffII = []
        usedOffsetsII = []

        for arg, argII in zip(self.args,argsII):

            if arg in fieldDicI.keys():
                  dim = 1
                  ffI.append(fieldDicI[arg])
                  usedOffsetsI.append(fieldOffsetsI[arg])
                  ffII.append(fieldDicII[argII])
                  usedOffsetsII.append(fieldOffsetsII[argII])
            else:
                  dim =0
                  #field = fieldDic[arg+"_0"]
                  for i in range(3):
                      if arg+"_"+str(i) in fieldDicI:
                         dim += 1
                         ffI.append(fieldDicI[arg+"_"+str(i)])
                         usedOffsetsI.append(fieldOffsetsI[arg+"_"+str(i)])
                         ffII.append(fieldDicII[argII+"_"+str(i)])
                         usedOffsetsII.append(fieldOffsetsI[argII+"_"+str(i)])

                      else:
                         break


        for cpt,nidI in enumerate(usedNodesMeshI):
                posI = nodesI[cpt,:]
                entries = oc2.find_within_range_cube(nodesI[cpt], self.tol)
                for entry in entries:
                    posII, nidII = entry
                    print(posII,nidII)
                    dist = np.linalg.norm(posII - posI)
                    if dist > self.tol:
                        continue



                    for i in range(len(ffI)):
                        firstOff = usedOffsetsI[i]
                        firstNumbering = ffI[i].numbering.GetDofOfPoint(nidI)+firstOff
                        secondOff = usedOffsetsII[i]
                        secondNumbering = ffII[i].numbering.GetDofOfPoint(nidI) + secondOff
                        CH.AddFactor(firstNumbering,1)
                        CH.AddFactor(secondNumbering,-1)
                        CH.NextEquation()
                        #print("Adding ",(firstNumbering,secondNumbering))

        return CH

class KRConformalTieScalar(KRConformalTieVector):
    def __init__(self):
        super(KRConformalTieScalar,self).__init__()



def CheckIntegrityKRConformalTieScalar(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    from BasicTools.FE.FETools import PrepareFEComputation

    mesh = CreateCube()
    space, numberings, offset, _ = PrepareFEComputation(mesh, numberOfComponents=3)

    print(mesh)
    obj = KRConformalTieVector()
    obj.From([-1,0,0])
    obj.To([0,0,0])
    obj.AddArg("temp")
    obj.SideI("X0")
    obj.SideII("X1")


    from BasicTools.FE.Fields.FEField import FEField
    fields =  []
    fields.append(FEField("temp",mesh=mesh,space=space, numbering=numberings[0]) )

    CH = obj.GenerateEquations(mesh,fields)
    CH.SetNumberOfDofs(numberings[0]["size"])
    mat, dofs = CH.ToSparse()
    #print(CH.cols)
    #print(CH.rows)
    #print(CH.vals)
    return 'ok'


def CheckIntegrityKRConformalTieVector(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    from BasicTools.FE.FETools import PrepareFEComputation

    mesh = CreateCube()
    space, numberings, offset, _ = PrepareFEComputation(mesh, numberOfComponents=3)

    print(mesh)
    obj = KRConformalTieVector()
    obj.From([-1,0,0])
    obj.To([0,0,0])
    obj.AddArg("temp")
    obj.SideI("X0")
    obj.SideII("X1")


    from BasicTools.FE.Fields.FEField import FEField
    fields =  []
    for x in range(3):
        fields.append(FEField("temp_"+str(x),mesh=mesh,space=space, numbering=numberings[x]) )

    CH = obj.GenerateEquations(mesh,fields)
    CH.SetNumberOfDofs(numberings[0]["size"]*3)
    mat, dofs = CH.ToSparse()
    #print(CH.cols)
    #print(CH.rows)
    #print(CH.vals)
    print(dofs)
    print(mat.toarray())
    return 'ok'


def CheckIntegrity(GUI=False):
    totest = [CheckIntegrityKRConformalTieScalar,
              CheckIntegrityKRConformalTieVector]

    for f in totest:
        print(str(f))
        res = f(GUI)
        if res.lower() != "ok":
            return res

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
