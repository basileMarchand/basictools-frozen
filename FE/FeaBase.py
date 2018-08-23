# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.FE.LinearSolver import LinearProblem


class FeaBase(BaseOutputObject):
    def __init__(self,dim=3, size= 1):
        super(FeaBase,self).__init__()
        self.mesh = None
        self.space = None
        self.numbering = None



        # tuple (zone,dof,val)
        self.dirichletBC = []
        self.rhsfixed = None
        self.sol = None
        self.solver = LinearProblem()

    def SetMesh(self,mesh):
        self.mesh = mesh

    def ComputeDofNumbering(self,tag=None):
        from BasicTools.FE.DofNumbering import ComputeDofNumbering
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

        if self.space is  LagrangeSpaceGeo and tag is None:
            # fast generation of the numbering based on the physical Geo space
            # warning !!!!!!
            # will add numbering for lonely nodes also
            self.numbering = ComputeDofNumbering(self.mesh,self.space,fromConnectivity =True,dofs=self.numbering)
        else :
            self.numbering = ComputeDofNumbering(self.mesh,self.space,fromConnectivity =False,tag=tag,dofs=self.numbering)

    def AddDirichlet(self,zone,dofs,val):
        self.dirichletBC.append((zone,dofs,val))

    def ApplyBC(self,K,rhs):
        #for every bc we generate the list of dofs

        self.dofs = np.ones(len(rhs), dtype=bool)
        if self.sol is None or len(rhs) != len(self.sol) :
            self.sol = np.zeros(len(rhs),dtype=np.float)

        for zone,dof,val in self.dirichletBC:
            offset = 0
            for i in range(dof):
                offset += self.numbering["size"]

            if zone in self.mesh.nodesTags:
                 nids = self.mesh.nodesTags[zone].GetIds()
                 dofsids = nids + offset

                 self.dofs[dofsids] = False
                 self.sol[dofsids] = val
            else :
                for name,data in self.mesh.elements.items():
                    if zone in data.tags:
                        elids = data.tags[zone].GetIds()
                        dofsids = np.unique(self.numbering[name][elids,:].ravel()) + offset

                        self.dofs[dofsids] = False
                        self.sol[dofsids] = val

        [cleanK, self.rhsfixed] = deleterowcol(K, np.logical_not(self.dofs), np.logical_not(self.dofs), self.sol)

        return cleanK, self.ApplyBCF(rhs)

    def ApplyBCF(self,rhs):
        return rhs[self.dofs]-self.rhsfixed[self.dofs]

    def Reset(self):
        self.sol = None
        self.solver.u = None
        pass

    def Solve(self,cleanK,cleanff):
        #self.PrintDebug(cleanK.tocsc())
        #self.PrintDebug(cleanff)
        self.solver.SetOp(cleanK.tocsc())
        self.Resolve(cleanff)

    def Resolve(self,cleanff):
        res = self.solver.Solve(cleanff)
        self.sol[self.dofs] = res

    def PushVectorToMesh(self,onNodes,field,name):
        F = field.view()
        if onNodes:
            if self.mesh.GetNumberOfNodes() == F.size:
                self.mesh.nodeFields[name] = np.zeros((self.mesh.GetNumberOfNodes()),dtype=float)
                self.mesh.nodeFields[name][self.numbering["doftopointLeft"]] =  F[self.numbering["doftopointRight"]]
                return
            if self.mesh.GetNumberOfNodes()*3 != F.size :
                raise Exception("incompatible field")
            F.shape = (3,F.size/3)
            F = F.T

            self.mesh.nodeFields[name] = np.zeros((self.mesh.GetNumberOfNodes(),3),dtype=float)
            self.mesh.nodeFields[name][self.numbering["doftopointLeft"],:] =  F[self.numbering["doftopointRight"],:]
            #print(self.numbering["doftopointLeft"])
            #print(self.numbering["doftopointRight"])
            #for cpt in range(len(self.numbering["doftopointLeft"])):
                #i = self.numbering["doftopointLeft"][cpt]
                #j = self.numbering["doftopointRight"][cpt]
                #print(i),
                #print(j)
                #self.mesh.nodeFields[name][i,:] =  F[j,:]
        else:
            raise
            if self.mesh.GetNumberOfElements(3)*3 != F.size :
                raise Exception("incompatible field")
            F.shape = (3,self.mesh.GetNumberOfElements(3))
            F = F.T
            self.mesh.elemFields[name] = np.zeros((self.mesh.GetNumberOfElements(),3),dtype=float)
            self.mesh.elemFields[name][self.numbering["doftopointLeft"],:] =  F[self.numbering["doftopointRight"],:]



def deleterowcol(A, delrow, delcol, fixedValues ):
    # Assumes that matrix is in symmetric csc form !

    rhs = A.dot(fixedValues)
    #keep = np.delete (np.arange(0, m), delrow)
    A = A[np.logical_not(delrow) , :]
    #keep = np.delete (np.arange(0, m), delcol)
    A = A[:, np.logical_not(delcol)]

    return [A, rhs]



def CheckIntegrity(GUI=False):
    FeaBase()
    from scipy.sparse import csr_matrix

    fv = np.array([1,2]).T
    fv = np.array([[1,2],]).T
    mask = np.zeros(2)
    mask[1] = True

    for sp in [True,False]:
        if sp:
            K = csr_matrix([[1, 2], [3, 4]])
        else:
            K = np.array([[1, 2], [3, 4]]);

        A,rhs = deleterowcol(K, mask, mask, fv )
        print("using Sparse : " + ("True" if sp else "False" ) )
        print("Vals")
        print(A)
        print("--")
        print(rhs)
        print("Types")
        print(type(A))
        print(type(rhs))
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
