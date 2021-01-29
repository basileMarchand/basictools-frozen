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


def _available_algorithms():
    result = ["Direct", "CG", "lsqr", "gmres","lgmres"]
    try:
        import sksparse.cholmod
        result.append("cholesky")
    except ImportError:
        pass

    try:
        import pyamg
        result.append("AMG")
    except ImportError:
        pass

    try:
        import BasicTools.Linalg.EigenSolver
        result.extend(("EigenCG", "EigenLU","EigenBiCGSTAB","EigenSPQR"))
    except ImportError:
        pass

    return tuple(result)


class LinearProblem(BOO):
    def __init__(self):
        super(LinearProblem,self).__init__()
        self.op = None   # the operator to solve
        self.solver = None
        self.type = "EigenCG"
        self.u = None
        self.tol = 1.e-6

        self.constraints = ConstraintsHolder()
        self.SetAlgo(self.type)

    def SetTolerance(self,tol):
        self.tol = tol

    def HasConstraints(self):
        return self.constraints.numberOfEquations > 0

    def GetNumberOfDofs(self):
        return self.op.shape[0]

    def ComputeProjector(self,mesh,fields):
        self.constraints.ComputeConstraintsEquations(mesh,fields)
        ndofs = np.sum([f.numbering["size"] for f in fields])
        print("Number Of dofs" , ndofs)

        self.constraints.SetNumberOfDofs(ndofs)
        self.constraints.UpdateCSystem()
        self.constraints.SetOp( self.originalOp   )
        print("Number Of dofs " , self.originalOp.shape)
        self.op = self.constraints.GetReducedOp()

    # you must set SetAlgo before setting the Op
    def SetOp(self, op):
        self.PrintVerbose('In SetOp (type:' +str(self.type) + ')')

        self.originalOp = op

        if self.HasConstraints():
            self.PrintVerbose(" With constraints (" +str(self.constraints.numberOfEquations) + ")")
            self.PrintVerbose(" Treating constraints using "+ str(type(self.constraints.method)) )
            self.constraints.SetNumberOfDofs(op.shape[1])
            self.PrintDebug(" Setting Op" )
            self.constraints.SetOp( op  )
            self.PrintDebug(" UpdateCSystem()" )
            self.constraints.UpdateCSystem()
            self.PrintDebug(" GetCOp ")
            op = self.constraints.GetCOp()
            self.PrintVerbose('Constraints treatememnt Done ')

        self.op = op

        if self.type in [ "CG", "gmres","lsqr", "lgmres" ]:
            self.solver = None
        elif self.type == "Direct":
            self.PrintDebug('Starting Factorisaton')
            self.solver = sps.linalg.factorized(op)
        elif self.type == "AMG":
            import pyamg
            self.solver = pyamg.ruge_stuben_solver(op.tocsr())
        elif self.type == "cholesky":
            from sksparse.cholmod import cholesky
            self.solver = cholesky(op)
        elif len(self.type) >= 5 and  self.type[0:5] == "Eigen":
            import BasicTools.Linalg.EigenSolver as EigenSolver
            self.solver = EigenSolver.EigenSolvers()
            self.solver.SetSolverType(self.type[5:])
            self.solver.SetTolerance(self.tol)
            self.solver.SetOp(op)
        else:
            raise(Exception("Error setting the operator"))

        self.PrintDebug('In SetOp Done')

    def SetAlgo(self, algoType, withErrorIfNotFound=False):
        if algoType in self.GetAvailableAlgorithms():
            self.type = algoType
        else:
            if withErrorIfNotFound:
                raise ValueError(TF.InRed("Error : ") + "Type not allowed ("+algoType+")") #pragma: no cover
            else:
                defaultIfError = "CG"
                print(f"Solver {algoType} unavailable, falling back to solver {defaultIfError} instead.")
                self.SetAlgo(defaultIfError)

    def Solve(self, rhs):
        u0 = None
        if self.HasConstraints():
            rhs = self.constraints.GetCRhs(rhs)
            if self.u is not None:
                u0 =  self.constraints.RestrictSolution(self.u)
        else:
            u0 = self.u

        if u0 is None:
           u0 = np.zeros_like(rhs)

        rhs = rhs.squeeze()
        u0 = u0.squeeze()

        self.PrintDebug('In LinearProblem Solve')
        if self.type in ["Direct", "cholesky"]:
            self.u = self.solver(rhs)
        elif self.type == "AMG":
            self.u = self.solver.solve(rhs,x0=u0,tol=self.tol)
        elif len(self.type) >= 5 and  self.type[0:5] == "Eigen":
            self.u = self.solver.Solve(rhs,u0)
        elif self.type == "CG":
            diag = self.op.diagonal()
            diag[diag == 0] = 1.0
            M = sps.dia_matrix((1./diag,0), shape=self.op.shape)

            norm = np.linalg.norm(rhs)
            if self.u is None:
                res = spslin.cg(self.op, rhs/norm, M = M,tol = self.tol)
            else:
                res = spslin.cg(self.op, rhs/norm, M = M, x0 = u0/norm,tol = self.tol)
            self.u = res[0][np.newaxis].T*norm
            self.u = self.u[:,0]

            if res[1] > 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Convergence to tolerance not achieved"))) #pragma: no cover
            if res[1] < 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Illegal input or breakdown"))) #pragma: no cover
        elif self.type == "gmres":
            diag = self.op.diagonal()
            diag[diag == 0] = 1.0
            M = sps.dia_matrix((1./diag,0), shape=self.op.shape)
            self.u = spslin.gmres(self.op, rhs, x0 = u0,M = M,tol = self.tol, atol= self.tol)[0]
        elif self.type == "lgmres":
            diag = self.op.diagonal()
            diag[diag == 0] = 1.0
            M = sps.dia_matrix((1./diag,0), shape=self.op.shape)
            self.u = spslin.lgmres(self.op, rhs, x0 = u0,M = M,tol = self.tol, atol= self.tol)[0]
        elif self.type == "lsqr":
            self.u = spslin.lsqr(self.op, rhs, atol=self.tol, btol=self.tol, x0 = u0)[0]
        else:
            raise(Exception("Error solving system"))


        self.PrintDebug("Done Linear solver")

        if self.HasConstraints():
            res = np.empty(self.constraints.nbdof, dtype=float)
            self.u = self.constraints.RestoreSolution(self.u)

        return self.u

    __available_algorithms = _available_algorithms()

    @classmethod
    def GetAvailableAlgorithms(cls):
        return cls.__available_algorithms

def CheckSolver(GUI,solver):
    print("Solver "+ str(solver))
    LS = LinearProblem ()
    LS.SetGlobalDebugMode()
    LS.SetAlgo(solver)

    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]])) )
    sol = LS.Solve(np.array([[1.],[2.]]))
    # second run
    sol = LS.Solve(np.array([[1.],[2.]]))

    if abs(sol[0] - 2.) >1e-15  :
        print(sol[0]-2) #pragma: no cover
        raise Exception() #pragma: no cover
    if abs(sol[1] - 2.) > 1e-15 : raise Exception()

    if solver == "EigenSPQR" :

        A = sps.csc_matrix(np.array([[0,0,0,1,1],
                                     [0,0,0,1,1],
                                     [0.51,0.5,0.5,0,10],
                                     [0.5,0.5,0.5,0,10],
                                     [0,1,0,0,5],
                                     ]))

        def QRTest(A):
            LS.SetTolerance(1e-5)
            print(A.toarray())
            LS.SetOp(A)
            rank = LS.solver.GetSPQRRank()
            print("Rank",rank )
            Q = LS.solver.GetSPQR_Q()
            print("Q")
            print(Q.toarray())
            R = LS.solver.GetSPQR_R()
            print("R")
            print(R.toarray())
            P = LS.solver.GetSPQR_P()
            print("P" ,P)
            print("-------------------------------------------------")
            print("AP")
            print(A.toarray()[:,P])
            print("QR")
            print(Q.tocsr().dot(R.tocsr().toarray() ) )
            print("norm A[:,P]-QR")
            print(np.linalg.norm(A.toarray()[:,P] - Q.tocsr()[:,0:rank].dot(R.tocsr()[0:rank,:].toarray() ) ) )
            print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            print("Q réduit")
            print(Q.tocsr()[:,0:rank].toarray())
            print("R réduit")
            print(R.tocsr()[0:rank,0:rank].toarray())
            print("QR reduit")
            print(Q.tocsr()[:,0:rank].dot(R.tocsr()[0:rank,0:rank].toarray() ) )
            print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            error = np.linalg.norm(A.toarray()[:,P[0:rank]] - Q.tocsr()[:,0:rank].dot(R.tocsr()[0:rank,0:rank].toarray() ) )
            if abs(error) > 1e-10:
                raise
            print(error)


            print ("OK EigenSPQR"  )

        QRTest(A.T)



def CheckIntegrity(GUI=False):

    solvers = _available_algorithms()

    for s in solvers:
        CheckSolver(GUI,s)

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
