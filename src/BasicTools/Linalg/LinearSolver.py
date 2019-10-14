# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

__author__ = "Felipe Bordeu"

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spslin

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from BasicTools.Helpers.TextFormatHelper import TFormat as TF

from BasicTools.Linalg.ConstraintsHolder import ConstraintsHolder
# you must set SetAlgo before setting the Op

class LinearProblem(BOO):
    def __init__(self):
        super(LinearProblem,self).__init__()
        self.op = None   # the operator to solve
        self.solver = None
        self.type = "CG"
        self.u = None
        self.tol = 1.e-6

        self.constraints = ConstraintsHolder()

    def HasConstraints(self):
        return self.constraints.numberOfEquations > 0

    def GetNumberOfDofs(self):
        return self.op.shape[0]

    def ComputeProjector(self,mesh,fields):
        self.constraints.ComputeConstraintsEquations(mesh,fields)
        ndofs = np.sum([f.numbering["size"] for f in fields])
        print("Number Of dofs" , ndofs)

        self.constraints.SetNumberOfDofs(ndofs)
        self.constraints.ComputeProjector()
        self.constraints.SetOp( self.originalOp   )
        print("Number Of dofs " , self.originalOp.shape)
        self.op = self.constraints.GetReducedOp()


    def SetOp(self, op):
        self.PrintVerbose('In SetOp (type:' +str(self.type) + ')')

        self.originalOp = op

        if self.HasConstraints():
            self.PrintVerbose('With constraints (' +str(self.constraints.numberOfEquations) + ')')
            self.constraints.SetNumberOfDofs(op.shape[1])
            self.constraints.SetOp( op  )
            self.constraints.ComputeProjector()
            op = self.constraints.GetReducedOp()
            self.PrintVerbose('Constraints treatememnt Done ')

        self.op = op


        if self.type == "Direct":
            self.PrintDebug('Starting Factorisaton')
            self.solver = sps.linalg.factorized(op)
        elif self.type == "CG":
            self.solver = None
            pass
        elif self.type == "lstsq":
            self.solver = None
            pass
        elif self.type == "gmres":
            self.solver = None
            pass
        elif self.type == "cholesky":
            from sksparse.cholmod import cholesky
            self.solver = cholesky(op)
        elif len(self.type) >= 5 and  self.type[0:5] == "Eigen":
            self.PrintVerbose('Loading Eigensolver ')

            import BasicTools.FE.EigenSolver as EigenSolver
            self.PrintVerbose('Loading Eigensolver Done')
            self.solver = EigenSolver.EigenSolvers()
            #print(self.type[5:])
            if self.type[5:] =="CG":
                self.solver.SetSolverType(1)
            elif self.type[5:] =="LU":
                self.solver.SetSolverType(2)
            else:
                raise(Exception("Solver not found"))
            self.solver.SetOp(op)
        else:
            raise(Exception("Error"))

        self.PrintDebug('In SetOp Done')

    def SetAlgo(self, algoType):
        if algoType in  ["Direct" ,"CG", "lstsq","cholesky" ,"EigenCG","EigenLU","gmres"] :
            self.type = algoType
        else:
            raise(ValueError(TF.InRed("Error : ") + "Type not allowed ("+algoType+")"  ) )#pragma: no cover

    def Solve(self, rhs):
        if self.HasConstraints():
            rhs = self.constraints.GetReducedRhs(rhs)

        self.PrintDebug('In LinearProblem Solve')
        if self.type == "Direct":
            self.u = self.solver(rhs)
        elif len(self.type) >= 5 and  self.type[0:5] == "Eigen":
            #print("solve using eigen")
            self.u = self.solver.Solve(rhs)
        elif self.type == "CG":
            diag = self.op.diagonal()
            diag[diag == 0] = 1.0
            M = sps.dia_matrix((1./diag,0), shape=self.op.shape)

            norm = np.linalg.norm(rhs)
            if self.u is None:
                res = spslin.cg(self.op, rhs/norm, M = M,tol = self.tol)
            else:
                res = spslin.cg(self.op, rhs/norm, M = M, x0 = self.u/norm,tol = self.tol)
            self.u = res[0][np.newaxis].T*norm
            self.u = self.u[:,0]

            if res[1] > 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Convergence to tolerance not achieved"))) #pragma: no cover
            if res[1] < 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Illegal input or breakdown"))) #pragma: no cover
        elif self.type == "gmres":
            self.u = np.linalg.gmres(self.op, rhs,tol = self.tol)[0]
        elif self.type == "lstsq":
            self.u = np.linalg.lstsq(self.op, rhs)[0]
        elif self.type == "cholesky":
            self.u = self.solver(rhs)

        self.PrintDebug("Done Linear solver")

        if self.HasConstraints():
            res = np.empty(self.constraints.nbdof, dtype=float)
            res[self.constraints.masters] = self.u
            res[self.constraints.slaves] = self.constraints.slaveValues
            self.u = res

        return self.u

def CheckIntegrity(GUI=False):

    LS = LinearProblem ()

    LS.SetAlgo('CG')
    LS.SetGlobalDebugMode()

    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]])) )
    sol = LS.Solve(np.array([[1],[2]]))
    # second run
    sol = LS.Solve(np.array([[1],[2]]))
    if abs(sol[0] - 2.) >1e-15  :
        print(sol[0]-2) #pragma: no cover
        raise Exception() #pragma: no cover
    if abs(sol[1] - 2.) > 1e-15 : raise Exception()

    LS.SetAlgo("Direct")

    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]])))
    sol = LS.Solve(np.array([[1],[2]]))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()

    LS.SetAlgo("cholesky")

    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]])))
    sol = LS.Solve(np.array([[1],[2]]))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()


    b = np.array([1.,2.])
    matrix = sps.csc_matrix(np.array([[0.5,0],[0,1.]]))

    LS.SetAlgo("EigenLU")
    LS.SetOp(matrix)
    sol = LS.Solve(b)
    #print("sol " +str(sol))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()

    LS.SetAlgo("EigenCG")
    LS.SetOp(matrix)
    sol = LS.Solve(b)
    #print("sol " +str(sol))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()
#    if sol[0] != 2. : raise Exception()
#    if sol[1] != 2. : raise Exception()

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
