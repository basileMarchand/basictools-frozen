# -*- coding: utf-8 -*-

#import numpy as np
from BasicTools.FE.KR.KRBase import KRBase
from BasicTools.FE.ProblemData import Transform
import numpy as np

class KRBlock(KRBase):

    def __init__(self):
        super(KRBlock,self).__init__()
        self.type = "Block"
        self.constraintDiretions = "Final" # ["Global","Local","Original","Final"]

        self.OriginSystem = Transform()
        self.FinalSystem = Transform()

    def CompluteDisplacement(self,pos):
        disp = self.FinalSystem.ApplyInvTransform(self.OriginSystem.ApplyTransform(pos))
        return disp


    def GetConstrainedDirections(self,pos=None,direction=None):
        res = []

        for x,y in zip([0,1,2],self.blockDirections ):
            if y:
                res.append(self.GetDirections(x,pos,direction))
        return   res

    def GetDirections(self,i,pos=None,direction=None):
        if self.constraintDiretions == "Global":
            res = np.zeros(3)
            res[i] = 1
            return res
        elif self.constraintDiretions == "Original":
            return self.OriginSystem.GetOrthoNormalBase().GetDirection(i,pos,direction)
        elif self.constraintDiretions == "Final":
            return self.FinalSystem.GetOrthoNormalBase().GetDirection(i,pos,direction)
        elif self.constraintDiretions == "Local":
            if i == 0:
                return direction / np.linalg.norm(direction)
            #Create a local base based on the direction vector
            # need a second vector to generate a connsistent base(???)

            raise
        else:
            raise

    def From(self,offset=None,first=None):
        if offset is not None:
            self.OriginSystem.SetOffset(offset)
        if first is not None:
            self.OriginSystem.SetFirst(first)
        return self

    def To(self,offset=None,first=None):
        if offset is not None:
            self.FinalSystem.SetOffset(offset)
        if first is not None:
            self.FinalSystem.SetFirst(first)
        return self



    def GenerateEquations(self,mesh,fields,CH=None):

        if CH is None:
            from BasicTools.FE.ConstraintsHolder import ConstraintsHolder
            CH = ConstraintsHolder()

        # for the moment the aproximation spaces of a vector must be homogenious
        # (the same in every direction)

        #compute the offsets
        offsets = []


        fieldDic = {f.name:f for f in fields }
        fieldOffsets = { }

        cpt = 0
        for field in fields:
            offsets.append(cpt)
            fieldOffsets[field.name] = cpt
            cpt += field.numbering["size"]

        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

        for zone in self.on:
          for name,data in mesh.elements.items():
            geoSpace = LagrangeSpaceGeo[name]
            for arg in self.args:
              #field = usedFields[0]
              if arg in fieldDic.keys():
                  dim = 1
                  field = fieldDic[arg]
                  raise
              else:
                  dim = 3
                  field = fieldDic[arg+"_0"]

              sp = field.space[name]
              nbsf = sp.GetNumberOfShapeFunctions()

              if zone in data.tags:
                elids = data.tags[zone].GetIds()
                for elid in elids:
                    for i in range(nbsf):
                        dofid = field.numbering[name][elid,i]
                        if dim == 1:
                            CH.AddFactor(dofid+offsets[0],1)
                            CH.AddConstant()
                            raise
                            CH.NextEquation()
                        else:
                            parampos = sp.posN[i,:]
                            valN = geoSpace.GetShapeFunc(parampos)
                            pos = np.dot(valN ,mesh.nodes[data.connectivity[elid,:],:]).T
                            disp =  self.CompluteDisplacement(pos) -pos
                            dirs = self.GetConstrainedDirections(pos,disp)
                            for dirToBlock in dirs:
                                CH.AddEquationSparse(dofid+offsets,dirToBlock,np.dot(disp,dirToBlock)  )

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
        res += str(self.OriginSystem) + "\n"
        res += "---- Final ----------"
        res += str(self.FinalSystem)
        return res

