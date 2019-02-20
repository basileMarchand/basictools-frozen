# -*- coding: utf-8 -*-

import numpy as np
from BasicTools.FE.KR import KRBase

class KRMasterSlaves(KRBase):
    def __init__(self):
        super(KRMasterSlaves,self).__init__()
        self.type = "KRMasterSlaves"
        from BasicTools.FE.ProblemData import Transform
#        self.OriginSystem = Transform()
        self.FinalSystem = Transform()
        self.constraintDiretions = "Final" # ["Global","Local","Original","Final"]

    def MasterNode(self,tag=None):
        self.master = tag


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

        #compute the ofsets
        offsets = []
        cpt = 0

        fieldDic = {f.name:f for f in fields }
        fieldOffsets = { }

        for field in fields:
            offsets.append(cpt)
            cpt += field.numbering["size"]
            offsets.append(cpt)
            fieldOffsets[field.name] = cpt

        #recover the dofs of the master no de
        self.master
        nids = mesh.nodestags[self.master].GetIds()
        if len(nids) > 1 :
            raise Exception("master tag has more than 1 nodes")

        initPos = mesh.nodes[nids[0],:]

        usedFields = []

        if self.master in fieldDic.keys() :
            # we have a scalar field (like temperature)
            masterDofs = [ np.array([ fieldDic[self.master].numbering['almanac'][('P',x,None)] for x in nids])+fieldOffsets[self.master] ]
            usedFields = [ fieldDic[self.master] ]
        else :
            # we work with a tensor quantity (displacement
            masterDofs = []
            for sufix in ["_0", "_1" ,"_2" ]:
                if self.master + sufix in fieldDic.keys() :
                    masterDofs.append([ np.array([ fieldDic[self.master + sufix].numbering['almanac'][('P',x,None)] for x in nids])+fieldOffsets[self.master + sufix] ] )
                    usedFields.append(fieldDic[self.master + sufix])
                    if usedFields[0] != usedFields[-1]:
                        raise Exception("Cant treat kinematic Relation using different approximation spaces")
                else:
                    break

        #from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
        #geoSpace = LagrangeSpaceGeo

        dim = len(usedFields)
        for zone in self.on:
          for name,data in mesh.elements.items():
            for arg in self.args:
              field = usedFields[0]
              sp = field.space[name]
              nbsf = sp.GetNumberOfShapeFunctions()

              if zone in data.tags:
                elids = data.tags[zone].GetIds()
                for elid in elids:
                    for i in range(nbsf):

                        dofid = field.numbering[name][elid,i]

                        if dim == 1:
                            CH.AddFactor(dofid+offsets[0],1)
                            CH.AddFactor(masterDofs[0],-1)
                            CH.NextEquation()
                        else:
                            # the posicion of
                            pos = sp.GetPosOfShapeFunction(i,mesh.nodes[data.connectivity[elid,:],:] )
                            # vector from the master to the final point
                            disp = pos - initPos
                            dirToBlock = self.GetConstrainedDirections(pos,disp)
                            for dtb in dirToBlock:
                                for d in range(dim):
                                    CH.AddFactor(dofid+offsets[d],dtb[d])
                                    CH.AddFactor(masterDofs[d],-dtb[d])
                                CH.NextEquation()
                            cpt+=1
        return CH


    def __str__(self):
        res = self.type + " "
        if len(self.arg) > 1:
            res += "("
        res += " and ".join(str(x)+"."+str(y) for x,y in self.arg)

        if len(self.arg) > 1:
            res += ")"
        res += "_On(" + ",".join(self.on) + ")"
        res += "_To(" + ",".join(self.to) + ")"

        return res

def CheckIntegrity(GUI=False):
    obj = KRMasterSlaves()
    obj.AddArg("u").On("Z0").Fix0().Fix1(False).Fix2(True)

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
