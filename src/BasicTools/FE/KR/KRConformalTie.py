# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
from scipy.spatial import KDTree

from BasicTools.Containers.Filters import ElementFilter
from BasicTools.FE.KR.KRBase import KRBaseVector
import BasicTools.Helpers.CPU  as CPU


class KRConformalTieVector(KRBaseVector):
    def __init__(self):
        super(KRConformalTieVector,self).__init__()
        self.tol = 0.000001
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

        offsets  , fieldOffsetsI, totalNumberOfDofsI  = self._ComputeOffsets(fields)
        offsetsII, fieldOffsetsII, totalNumberOfDofsII = self._ComputeOffsets(fieldsII)

        def FillOctree(mesh, tags, base ):

            usedNodes = set()

            ef = ElementFilter(mesh,tags=tags)
            for name,data, ids in ef:
                nids = data.GetNodesIdFor(ids)
                usedNodes.update(nids)

            usedNodes = list(usedNodes)
            nodes = base.ApplyTransform(mesh.nodes[usedNodes,:])

            return ( KDTree(nodes), usedNodes, nodes)

        kdTree1, usedNodesMeshI, nodesI  = FillOctree(meshI, self.on, self.originSystem.GetOrthoNormalBase())
        kdTree2, usedNodesMeshII, nodesII = FillOctree(meshII, self.onII, self.targetSystem.GetOrthoNormalBase())

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

        neighboors = kdTree2.query_ball_point(x=nodesI[:len(usedNodesMeshI)], r=self.tol, workers=CPU.GetNumberOfAvailableCpus())
        for cpt,nidI in enumerate(usedNodesMeshI):
            posI = nodesI[cpt,:]
            entries = neighboors[cpt]

            for entry in entries:
                nidII = entry
                posII = usedNodesMeshII[entry]
                dist = np.linalg.norm(posII - posI)
                if dist > self.tol:
                    continue

                for i in range(len(ffI)):
                    firstOff = usedOffsetsI[i]
                    firstNumbering = ffI[i].numbering.GetDofOfPoint(nidI)+firstOff
                    secondOff = usedOffsetsII[i]
                    secondNumbering = ffII[i].numbering.GetDofOfPoint(nidII) + secondOff
                    CH.AddFactor(firstNumbering,1)
                    CH.AddFactor(secondNumbering,-1)
                    CH.NextEquation()
                break

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
