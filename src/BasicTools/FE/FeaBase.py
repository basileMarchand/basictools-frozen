#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Linalg.LinearSolver import LinearProblem
from  BasicTools.Containers.UnstructuredMesh import AllElements


class FeaBase(BaseOutputObject):
    """
    Base class for a finit element solver, this class is experimental

    normaly a finit element solver has a mesh (slef.mesh), a solution vector
    (self.sol), a linear solver (self.solver), the dimensionality of the physical
    space (1D,2D,3D) (self.spaceDim), and the number of dofs to allocate the
    objects. All the other parts (asembly operator, IO). must be defined in the
    derived class

    """

    def __init__(self,spaceDim=3, size= 1):
        super(FeaBase,self).__init__()
        self.mesh = None

        self.sol = None
        self.solver = LinearProblem()
        self.spaceDim = spaceDim
        self.totalNumberOfDof = 0

    def SetMesh(self,mesh):
        """
        To set the mesh
        """
        self.mesh = mesh

    def ComputeDofNumbering(self,tag=AllElements):
        """
        This fuction must be eliminated (it uses self.space).
        """
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
        """
        To eliminate the solution vector and to reset the linear solver
        """
        self.sol = None
        self.solver.u = None
        pass

    def ComputeConstraintsEquations(self):
        """
        To computhe the cinematic relation in terms of dofs.

        The the cinematic relations are stored in the solver
        """
        self.solver.constraints.ComputeConstraintsEquations(self.mesh,self.unkownFields)

    def Solve(self,cleanK,cleanff):
        """
        Solve a linear system using the internal solver with the cinematic
        reations calculated previously
        """
        self.solver.SetOp(cleanK.tocsc())
        self.Resolve(cleanff)

    def Resolve(self,cleanff):
        """
        To solve a problem with the same tangent operator but with a different
        RHS term
        """
        res = self.solver.Solve(cleanff)
        self.sol = res


    def PushVectorToUnkownFields(self):
        """
        Function to extract fields from the solution vector and to put it into
        fields data
        """
        offset =0
        for f in self.unkownFields:
              f.data = self.sol[offset:offset+f.numbering["size"]]
              offset += f.numbering["size"]

    def PushVectorToMesh(self,onNodes,fieldValues,name,numberings,justReturn=False):
        """
        Function to push the data from the fields into a vector homogeneous to
        the mesh (for visualition for example). Entities with no dofs are filled
        with zeros. Use the justReturn to get the vector but do not put it in the mesh.
        """
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
            if justReturn :
                return res
            self.mesh.nodeFields[name] = res
        else:
            res = np.zeros((self.mesh.GetNumberOfElements(), ncomp),dtype=float)
            cpt = 0
            for i in range(ncomp):
                if len(numberings[i]["doftocellLeft"]) == 0:
                    print("Warning : PushVectorToMesh()  transfert vector is empty")
                res[numberings[i]["doftocellLeft"],i] =  F[cpt+numberings[i]["doftocellRight"]]
                cpt += numberings[i]["size"]
            if justReturn :
                return res
            self.mesh.elemFields[name] = res

def CheckIntegrity(GUI=False):
    FeaBase()

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
