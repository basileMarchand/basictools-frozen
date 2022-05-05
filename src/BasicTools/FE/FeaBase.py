# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Linalg.LinearSolver import LinearProblem

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

    def ComputeDofNumbering(self, elementFilter=None):
        """
        This fuction must be eliminated (it uses self.space).
        """


        from BasicTools.FE.DofNumbering import ComputeDofNumbering
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

        if self.space is  LagrangeSpaceGeo and elementFilter is None:
            # fast generation of the numbering based on the physical Geo space
            # warning !!!!!!
            # will add numbering for lonely nodes also
            self.numbering = ComputeDofNumbering(self.mesh,self.space,fromConnectivity =True,dofs=self.numbering)
        else :
            self.numbering = ComputeDofNumbering(self.mesh,self.space,fromConnectivity =False,dofs=self.numbering,elementFilter=elementFilter)

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


    def PushSolutionVectorToUnkownFields(self):
        """
        Function to extract fields from the solution vector and to put it into
        fields data
        """
        from BasicTools.FE.Fields.FieldTools import VectorToFEFieldsData
        VectorToFEFieldsData(self.sol,self.unkownFields)

    def PushUnkownFieldsToSolutionVector(self):
        """
        Function to extract from the Unkown fields a solution vector
        """
        from BasicTools.FE.Fields.FieldTools import FEFieldsDataToVector
        self.sol = FEFieldsDataToVector(self.unkownFields)
        self.solver.u = self.sol

def CheckIntegrity(GUI=False):
    FeaBase()

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
