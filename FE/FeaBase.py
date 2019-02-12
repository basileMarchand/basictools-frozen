# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Linalg.LinearSolver import LinearProblem
from  BasicTools.Containers.UnstructuredMesh import AllElements


class FeaBase(BaseOutputObject):
    def __init__(self,spaceDim=3, size= 1):
        super(FeaBase,self).__init__()
        self.mesh = None

        self.sol = None
        self.solver = LinearProblem()
        self.spaceDim = spaceDim
        self.totalNumberOfDof = 0

    def SetMesh(self,mesh):
        self.mesh = mesh

    def ComputeDofNumbering(self,tag=AllElements):
        from BasicTools.FE.DofNumbering import ComputeDofNumbering
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

        if self.space is  LagrangeSpaceGeo and tag is None:
            # fast generation of the numbering based on the physical Geo space
            # warning !!!!!!
            # will add numbering for lonely nodes also
            self.numbering = ComputeDofNumbering(self.mesh,self.space,fromConnectivity =True,dofs=self.numbering)
        else :
            self.numbering = ComputeDofNumbering(self.mesh,self.space,fromConnectivity =False,tag=tag,dofs=self.numbering)

    #def AddDirichlet(self,zone,dofs,val):
    #    self.dirichletBC.append((zone,dofs,val))

#    def ApplyBC(self,K,rhs):
#        #for every bc we generate the list of dofs
#
#        self.dofs = np.ones(len(rhs), dtype=bool)
#        if self.sol is None or len(rhs) != len(self.sol) :
#            self.sol = np.zeros(len(rhs),dtype=np.float)
#
#        for zone,dof,val in self.dirichletBC:
#            offset = 0
#            for i in range(dof):
#                offset += self.numbering["size"]
#
#            if zone in self.mesh.nodesTags:
#                 nids = self.mesh.nodesTags[zone].GetIds()
#                 nids = np.array([ self.numbering['almanac'][('P',x,None)] for x in nids])
#                 dofsids = nids + offset
#
#                 self.dofs[dofsids] = False
#                 self.sol[dofsids] = val
#            else :
#                for name,data in self.mesh.elements.items():
#                    if zone in data.tags:
#                        elids = data.tags[zone].GetIds()
#                        dofsids = np.unique(self.numbering[name][elids,:].ravel()) + offset
#
#                        self.dofs[dofsids] = False
#                        self.sol[dofsids] = val
#
#        [cleanK, self.rhsfixed] = deleterowcol(K, np.logical_not(self.dofs), np.logical_not(self.dofs), self.sol)
#
#        return cleanK, self.ApplyBCF(rhs)
#
#    def ApplyBCF(self,rhs):
#        return rhs[self.dofs]-self.rhsfixed[self.dofs]

    def Reset(self):
        self.sol = None
        self.solver.u = None
        pass

    def ComputeConstraintsEquations(self):
        self.solver.constraints.ComputeConstraintsEquations(self.mesh,self.unkownFields)

    def Solve(self,cleanK,cleanff):
        #self.PrintDebug(cleanK.tocsc())
        #self.PrintDebug(cleanff)
        self.solver.SetOp(cleanK.tocsc())
        self.Resolve(cleanff)

    def Resolve(self,cleanff):
        res = self.solver.Solve(cleanff)
        self.sol = res
        #self.sol[self.dofs] = res

    def PushVectorToUnkownFields(self):
        offset =0
        for f in self.unkownFields:
              f.data = self.sol[offset:offset+f.numbering["size"]]
              offset += f.numbering["size"]


    def PushVectorToMesh(self,onNodes,fieldValues,name,numberings):
        # vectorSize = 1 for scalar field
        # vectorSize = 3 for 3D vector field (example dep in 3D )
        F = fieldValues.view()

        ncomp = len(numberings)
        if onNodes:
            res = np.zeros((self.mesh.GetNumberOfNodes(), ncomp),dtype=float)
            cpt = 0
            for i in range(ncomp):
                if len(numberings[i]["doftopointLeft"]) == 0:
                    print("Warning : PushVectorToMesh()  transfert vector is empty")
                res[numberings[i]["doftopointLeft"],i] =  F[cpt+numberings[i]["doftopointRight"]]
                cpt += numberings[i]["size"]
            self.mesh.nodeFields[name] = res
        else:
            res = np.zeros((self.mesh.GetNumberOfElements(), ncomp),dtype=float)
            cpt = 0
            for i in range(ncomp):
                if len(numberings[i]["doftocellLeft"]) == 0:
                    print("Warning : PushVectorToMesh()  transfert vector is empty")
                res[numberings[i]["doftocellLeft"],i] =  F[cpt+numberings[i]["doftocellRight"]]
                cpt += numberings[i]["size"]
            self.mesh.elemFields[name] = res

def CheckIntegrity(GUI=False):
    FeaBase()

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
