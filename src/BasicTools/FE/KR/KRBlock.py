# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from BasicTools.FE.KR.KRBase import KRBaseScalar, KRBaseVector
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

class KRBlockScalar(KRBaseScalar):
    def __init__(self):
        super(KRBlockScalar,self).__init__()
        self.type = "BlockScalar"


    def GenerateEquations(self,mesh,fields,CH=None):
        CH = self._GetConstraintHolder(CH)

        offsets, fieldOffsets, totalNumberOfDofs  = self._ComputeOffsets(fields)

        fieldDic = {f.name:f for f in fields }

        ef = ElementFilter(mesh=mesh, tags=self.on )

        for name, data, elids in ef:
            geoSpace = LagrangeSpaceGeo[name]

            for arg in self.args:

                if arg in fieldDic.keys():
                    dim = 1
                    field = fieldDic[arg]
                    offsetUsed = fieldOffsets[arg]
                else:
                    raise(Exception("Field ("+str(arg)+") not found "))

                sp = field.space[name]
                nbsf = sp.GetNumberOfShapeFunctions()
                geoSpace.SetIntegrationRule(sp.posN,np.ones(nbsf) )

                treated = np.zeros(field.numbering["size"])
                for elid in elids:
                    for i in range(nbsf):
                        dofid = field.numbering[name][elid,i]
                        if treated[dofid]:
                           continue
                        treated[dofid] = True
                        parampos = sp.posN[i,:]
                        valN = geoSpace.GetShapeFunc(parampos)
                        xcoor = mesh.nodes[data.connectivity[elid,:],:]
                        pos = np.dot(valN ,xcoor).T
                        if dim == 1:
                            CH.AddFactor(dofid+offsetUsed,1)
                            CH.AddConstant(self.value(pos))
                            CH.NextEquation()
        return CH

    def __str__(self):
        res = self.type + " "
        if len(self.args) > 1:
            res += "("
        res += " and ".join(str(x) for x in self.args)
        if len(self.args) > 1:
            res += ")"

        return res

class KRBlockVector(KRBaseVector):

    def __init__(self):
        super(KRBlockVector,self).__init__()
        self.type = "BlockVector"

    def ComputeDisplacement(self,pos):
        disp = self.targetSystem.ApplyInvTransform(self.originSystem.ApplyTransform(pos))
        return disp



    def From(self,offset=None,first=None):
        if offset is not None:
            self.originSystem.SetOffset(offset)
        if first is not None:
            self.originSystem.SetFirst(first)
        return self

    def To(self,offset=None,first=None):
        if offset is not None:
            self.targetSystem.SetOffset(offset)
        if first is not None:
            self.targetSystem.SetFirst(first)
        return self

    def GenerateEquations(self,mesh,fields,CH=None):
        CH = self._GetConstraintHolder(CH)
        offsets, fieldOffsets, totalNumberOfDofs  = self._ComputeOffsets(fields)

        # for the moment the aproximation spaces of a vector must be homogenious
        # (the same in every direction)
        fieldDic = {f.name:f for f in fields }

        ef = ElementFilter(mesh=mesh, tags=self.on )


        for name, data, elids in ef:

        #for zone in self.on:
        #  for name,data in mesh.elements.items():
            geoSpace = LagrangeSpaceGeo[name]
            for arg in self.args:
                offsetUsed = []
                for i in range(3):
                    if arg + "_" + str(i) in fieldOffsets:
                        offsetUsed.append(fieldOffsets[arg + "_" + str(i)])
                    else:
                        break

                field = fieldDic[arg+"_0"]

                sp = field.space[name]
                nbsf = sp.GetNumberOfShapeFunctions()
                geoSpace.SetIntegrationRule(sp.posN,np.ones(nbsf) )

              #if zone in data.tags:
                #elids = data.tags[zone].GetIds()
                treated = np.zeros(field.numbering["size"])
                for elid in elids:
                    for i in range(nbsf):
                        dofid = field.numbering[name][elid,i]
                        if treated[dofid]:
                           continue
                        treated[dofid] = True
                        parampos = sp.posN[i,:]
                        valN = geoSpace.GetShapeFunc(parampos)
                        xcoor = mesh.nodes[data.connectivity[elid,:],:]
                        pos = np.dot(valN ,xcoor).T
                        if len(pos) == 2:
                            pos = [pos[0], pos[1], 0]
                        disp =  self.ComputeDisplacement(pos) - pos
                        if self.constraintDiretions == "Local":
                            Jack, Jdet, Jinv = geoSpace.GetJackAndDetI(i,xcoor)
                            normal = self.geoSpace.GetNormal(Jack)

                            dirs = []
                            for x,y in zip([0,1,2],self.blockDirections):
                                vec1 = np.lialg.cross([ 1,0,0 ],normal ) + np.lialg.cross([ 0,1,0 ],normal )
                                vec1 = vec1/np.linalg.norm(vec1)
                                if y and x== 0:
                                    dirs.append(normal)
                                elif y and x == 1:
                                    dirs.append(vec1)
                                elif y and x == 2:
                                    vec2 = np.linal.cross(normal,vec1)
                                    vec2 = vec2/np.linalg.norm(vec2)
                                    dirs.append(vec2)
                        else:
                            dirs = self.GetConstrainedDirections(pos,disp)

                        for dirToBlock in dirs:
                            dim = len(offsetUsed)
                            CH.AddEquationSparse(dofid+offsetUsed,dirToBlock[0:dim],np.dot(disp,dirToBlock)  )

        return CH

    def __str__(self):
        res = self.type + " "
        if len(self.args) > 1:
            res += "("
        res += " and ".join(str(x) for x in self.args)
        if len(self.args) > 1:
            res += ")"
        if self.constraintDiretions == "Global":
            res += "_BlockDir("+"".join([x if y else "" for x,y in zip(["x","y","z"],self.blockDirections )]) + ")"
        else :
            res += "_BlockDir("+"".join([str(self.GetDirections(x)) if y else "" for x,y in zip([0,1,2],self.blockDirections )]) + ")"
        res += "_On(" + ",".join(self.on) + ")"
        res += "_" + str(self.constraintDiretions)+ "\n"
        res += "---- Origin ---------"
        res += str(self.originSystem) + "\n"
        res += "---- Target ----------"
        res += str(self.targetSystem)
        return res

def CheckIntegrity(GUI=False):
    obj = KRBlockScalar()
    obj.AddArg("u").On("Z0")
    obj.SetValue(1.0)
    print(obj)

    obj = KRBlockVector()
    obj.AddArg("u").On("Z0").Fix0(True).Fix2(True)
    print(obj)



    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
