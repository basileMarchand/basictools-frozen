# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.FE.KR.KRBase import KRBaseScalar, KRBaseVector
from BasicTools.Containers.Filters import ElementFilter


class KRMasterNodeToSlaveScalar(KRBaseScalar):
    def __init__(self):
        super(KRMasterNodeToSlaveScalar,self).__init__()
        self.type = "KRMasterNodeToSlaveScalar"

    def MasterNode(self,tag=None):
        self.master = tag

    def GenerateEquations(self,mesh,fields,CH=None):
        CH = self._GetConstraintHolder(CH)
        offsets, fieldOffsets  = self._ComputeOffsets(fields)

        fieldDic = {f.name:f for f in fields }

        nids = mesh.nodesTags[self.master].GetIds()
        if len(nids) > 1 :
            raise Exception("master tag has more than 1 nodes")

        if len(self.args) > 1:
            raise(Exception("Cant treat only cases with one arg"))

        fieldname = self.args[0]
        field = fieldDic[fieldname]


        masterDofs = np.array([ field.numbering.GetDofOfPoint(x)  for x in nids])+fieldOffsets[fieldname]

        ef = ElementFilter(mesh=mesh, tags=self.on )
        for name, data, elids in ef:
            numbering = np.unique(field.numbering[name][elids,:].flatten())
            dofs = numbering+fieldOffsets[fieldname]
            for dof in dofs :
                if dof == masterDofs :
                    continue
                CH.AddFactor(dof,1)
                CH.AddFactor(masterDofs[0],-1)
                CH.NextEquation()

        return CH

class KRMasterNodeToSlaveVector(KRBaseVector):
    def __init__(self):
        super(KRMasterNodeToSlaveVector,self).__init__()
        self.type = "KRMasterNodeToSlaveVector"

    def MasterNode(self,tag=None):
        self.master = tag

    def GenerateEquations(self,mesh,fields,CH=None):
        CH = self._GetConstraintHolder(CH)
        offsets, fieldOffsets  = self._ComputeOffsets(fields)

        fieldDic = {f.name:f for f in fields }

        nids = mesh.nodesTags[self.master].GetIds()
        if len(nids) > 1 :
            raise Exception("master tag has more than 1 nodes")

        if len(self.args) > 1:
            raise(Exception("Cant treat only cases with one arg"))

        fieldname = self.args[0]

        masterDofs = []
        usedFields = []

        for sufix in range(3):
            cfieldname = fieldname + "_"+ str(sufix)

            masterDofs.append( fieldDic[cfieldname].numbering.GetDofOfPoint(nids[0]) + fieldOffsets[cfieldname]  )

            usedFields.append(fieldDic[cfieldname])


        initPos = mesh.nodes[nids[0],:]

        ef = ElementFilter(mesh=mesh, tags=self.on )

        field = usedFields[0]
        print(masterDofs)
        for name, data, elids in ef:
            print(name)
            sp = field.space[name]
            sp.Create()
            nbsf = sp.GetNumberOfShapeFunctions()
            nu = [x.numbering[name] for x in usedFields]
            for elid in elids:
                print(elid)
                for i in range(nbsf):
                    dofids = [x[elid,i] for x in nu]

                    # the posicion of
                    pos = sp.GetPosOfShapeFunction(i,mesh.nodes[data.connectivity[elid,:],:] )
                    # vector from the master to the final point
                    disp = pos - initPos
                    dirToBlock = self.GetConstrainedDirections(pos,disp)
                    for dtb in dirToBlock:
                        for d in range(3):
                            CH.AddFactor(dofids[d]+offsets[d],dtb[d])
                            CH.AddFactor(masterDofs[d],-dtb[d])
                        CH.NextEquation()
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

def CheckIntegrityKRMasterNodeToSlaveScalar(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    from BasicTools.FE.FETools import PrepareFEComputation

    mesh = CreateSquare()
    space, numberings, offset, _ = PrepareFEComputation(mesh, numberOfComponents=1)

    obj = KRMasterNodeToSlaveScalar()
    obj.MasterNode("x0y0")
    obj.AddArg("temp")
    obj.On("X1")


    from BasicTools.FE.Fields.FEField import FEField
    temp = FEField("temp",mesh=mesh,space=space, numbering=numberings[0])

    CH = obj.GenerateEquations(mesh,[temp])
    CH.SetNumberOfDofs(numberings[0]["size"])
    mat, dofs = CH.ToSparse()
    print(dofs)
    print(mat.toarray())

    #obj.AddArg("u").On("Z0").Fix0().Fix1(False).Fix2(True)

    return "ok"

def CheckIntegrityKRMasterNodeToSlaveVector(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    from BasicTools.FE.FETools import PrepareFEComputation

    mesh = CreateSquare()
    space, numberings, offset, _ = PrepareFEComputation(mesh, numberOfComponents=3)

    obj = KRMasterNodeToSlaveVector()
    obj.MasterNode("x0y0")
    obj.AddArg("temp")
    obj.On("X1").Fix0()


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


    return "ok"

def CheckIntegrity(GUI=False):
    totest = [CheckIntegrityKRMasterNodeToSlaveScalar,
              CheckIntegrityKRMasterNodeToSlaveVector]

    for f in totest:
        print(str(f))
        res = f(GUI)
        if res.lower() != "ok":
            return res

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
