# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spslin

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from BasicTools.Helpers.TextFormatHelper import TFormat as TF
from BasicTools.Linalg.ConstraintsHolder import ConstraintsHolder
from BasicTools.Helpers.Factory import Factory

class SolverFactory(Factory):
    _Catalog = {}
    _SetCatalog = set()
    def __init__(self):
        super().__init__()

def GetAvailableSolvers():
    return list(SolverFactory._Catalog.keys())

def RegisterSolverClass(name, classtype, constructor=None, withError = True):
    return SolverFactory.RegisterClass(name,classtype, constructor=constructor, withError = withError )

def RegisterSolverClassUsingName(cls):
    RegisterSolverClass(cls().name, cls)

class LinearSolverBase(BOO):
    def __init__(self):
        super().__init__()
        self.op = None   # the operator to solve
        self.originalOp = None
        self.solver = None
        self.name = ''
        self.u = None
        self.constraints = ConstraintsHolder()

        self._can_use_u0 = False

    def GetConstraints(self):
        return self.constraints

    def SetConstraints(self,constraints):
        self.constraints = constraints

    def HasConstraints(self):
        return self.constraints.numberOfEquations > 0

    def GetNumberOfDofs(self):
        return self.op.shape[0]

    def ComputeProjector(self, op ):
        self.PrintVerbose(" With constraints (" +str(self.constraints.numberOfEquations) + ")")
        self.PrintVerbose(" Treating constraints using "+ str(type(self.constraints.method)) )
        self.constraints.SetNumberOfDofs(op.shape[1])
        self.PrintDebug(" Setting Op" )
        self.constraints.SetOp( op  )
        self.PrintDebug(" UpdateCSystem()" )
        self.constraints.UpdateCSystem()
        self.PrintDebug(" GetCOp ")
        op = self.constraints.GetCOp()
        self.PrintVerbose('Constraints treatment Done ')
        return op

    def SetOp(self, op):
        self.PrintVerbose('In SetOp (type:' +str(self.name) + ')')
        self.originalOp = op

        if self.HasConstraints():
            op = self.ComputeProjector(op)

        self.op = op

        self.PrintDebug('In LinearSolver.SetOp(...) ')
        self._setop_imp(op)
        self.PrintDebug('In LinearSolver.SetOp(...) Done')

    def Solve(self, rhs, u0=None):
        if self.HasConstraints():
            rhs = self.constraints.GetCRhs(rhs.squeeze())
        rhs = np.atleast_1d(rhs.squeeze())

        if self.u is not None and self.originalOp is not None and len(self.u) != self.originalOp.shape[1]:
            self.u = None

        self.PrintDebug('In LinearProblem Solve ' + self.name)
        if self._can_use_u0 :
            if self.HasConstraints():
                if u0 is not None:
                    u0 =  self.constraints.RestrictSolution(u0)
            else:
                if self.u is not None:
                    u0 =  self.constraints.RestrictSolution(self.u)
            if u0 is None:
                u0 = np.zeros_like(rhs)
            u0 = np.atleast_1d(u0.squeeze())
            u = self._solve_imp(rhs,u0=u0)
        else:
            if u0 is not None and self.__print_warning_u0_ignored:
                print("u0 ignored for direct solvers")
                self.__print_warning_u0_ignored = False

            u = self._solve_imp(rhs,u0=None)
        self.PrintDebug("Done Linear solver "+str(u.shape))

        if self.HasConstraints():
            self.u = self.constraints.RestoreSolution(u)
        else:
            self.u = u

        return self.u

    def _setop_imp(self,op):
        pass


class LinearSolverIterativeBase(LinearSolverBase):
    def __init__(self):
        super().__init__()
        self.tol = 1.e-6
        self._can_use_u0 = True

    def SetTolerance(self,tol):
        self.tol = tol

class LinearSolverCG(LinearSolverIterativeBase):
    def __init__(self):
        super().__init__()
        self.name = "CG"

    def _solve_imp(self, rhs, u0):
        diag = self.op.diagonal()
        diag[diag == 0] = 1.0
        M = sps.dia_matrix((1./diag,0), shape=self.op.shape)

        norm = np.linalg.norm(rhs)
        if u0 is None:
            res = spslin.cg(self.op, rhs/norm, M = M,tol = self.tol, atol = self.tol)
        else:
            res = spslin.cg(self.op, rhs/norm, M = M, x0 = u0/norm,tol = self.tol, atol = self.tol)
        u = res[0][np.newaxis].T*norm
        u = u[:,0]

        if res[1] > 0 :
            self.Print(TF.InYellowBackGround(TF.InRed("Convergence to tolerance not achieved"))) #pragma: no cover
        if res[1] < 0 :
            self.Print(TF.InYellowBackGround(TF.InRed("Illegal input or breakdown"))) #pragma: no cover

        return u

RegisterSolverClassUsingName(LinearSolverCG)

class LinearSolvergmres(LinearSolverIterativeBase):
    def __init__(self):
        super().__init__()
        self.name = "gmres"

    def _solve_imp(self, rhs, u0):
        diag = self.op.diagonal()
        diag[diag == 0] = 1.0
        M = sps.dia_matrix((1./diag,0), shape=self.op.shape)
        return  spslin.gmres(self.op, rhs, x0 = u0,M = M,tol = self.tol, atol= self.tol)[0]

RegisterSolverClassUsingName(LinearSolvergmres)

class LinearSolverlsqr(LinearSolverIterativeBase):
    def __init__(self):
        super().__init__()
        self.name = "lsqr"

    def _solve_imp(self, rhs, u0):
        return spslin.lsqr(self.op, rhs, atol=self.tol, btol=self.tol, x0 = u0)[0]

RegisterSolverClassUsingName(LinearSolverlsqr)

class LinearSolverlgmres(LinearSolverIterativeBase):
    def __init__(self):
        super().__init__()
        self.name = "lgmres"

    def _solve_imp(self, rhs, u0):
        diag = self.op.diagonal()
        diag[diag == 0] = 1.0
        M = sps.dia_matrix((1./diag,0), shape=self.op.shape)
        return spslin.lgmres(self.op, rhs, x0 = u0,M = M,tol = self.tol, atol= self.tol)[0]

RegisterSolverClassUsingName(LinearSolverlgmres)

class LinearSolverlAMG(LinearSolverIterativeBase):
    def __init__(self):
        super().__init__()
        self.name = "AMG"

    def _setop_imp(self,op):
        import pyamg
        self._internal_solver = pyamg.ruge_stuben_solver(op.tocsr())

    def _solve_imp(self, rhs, u0):
        return self._internal_solver.solve(rhs,x0=u0,tol=self.tol)

try:
    import pyamg
    RegisterSolverClassUsingName(LinearSolverlAMG)
except:
    pass

class LinearSolverDirect(LinearSolverIterativeBase):
    def __init__(self):
        super().__init__()
        self.name = "Direct"

    def _setop_imp(self,op):
        self._internal_solver = sps.linalg.factorized(op)

    def _solve_imp(self, rhs, u0=None):
        return self._internal_solver(rhs)

RegisterSolverClassUsingName(LinearSolverDirect)

class LinearSolverCholesky(LinearSolverDirect):
    def __init__(self):
        super().__init__()
        self.name = "cholesky"

    def _setop_imp(self,op):
        from sksparse.cholmod import cholesky
        self._internal_solver = cholesky(op)

    def _solve_imp(self, rhs, u0= None):
        return self._internal_solver(rhs)

try:
    from sksparse.cholmod import cholesky
    RegisterSolverClassUsingName(LinearSolverCholesky)
except:
    pass

class LinearSolverEigen(LinearSolverIterativeBase):
    def __init__(self,subtype):
        super().__init__()
        self.SetSolver(subtype)
        import BasicTools.Linalg.NativeEigenSolver as NativeEigenSolver
        self.solver = NativeEigenSolver.CEigenSolvers()
        from BasicTools.Helpers.CPU import GetNumberOfAvailableCpus
        self.solver.ForceNumberOfThreads(GetNumberOfAvailableCpus())

    def SetSolver(self, subtype):
        self.name = "Eigen"+subtype
        self.subtype = subtype
        self.solver = None

    def _setop_imp(self,op):
        self.solver.SetSolverType(self.subtype)
        self.solver.SetTolerance(self.tol)
        self.solver.SetOp(op)

    def _solve_imp(self, rhs, u0=None):
        # for the eigen solver we must allocate on the python side
        if u0 is None:
            u0 = np.zeros_like(rhs)
        return self.solver.Solve(rhs,u0)

    def GetSPQRRank(self):
        return self.solver.GetSPQRRank()

    def GetSPQR_Q(self):
        return self.solver.GetSPQR_Q()

    def GetSPQR_R(self):
        return self.solver.GetSPQR_R()

    def GetSPQR_P(self):
        return self.solver.GetSPQR_P()

    @classmethod
    def GetAvailableSolvers(cls):
        return ['CG','LU','BiCGSTAB', 'SPQR']

try:
    import BasicTools.Linalg.NativeEigenSolver as NativeEigenSolver
    for eigenSubTypes  in LinearSolverEigen.GetAvailableSolvers():
        def GenerateEigenConstructor(type):
            return lambda x :  LinearSolverEigen(type)
        RegisterSolverClass("Eigen"+eigenSubTypes,LinearSolverEigen, GenerateEigenConstructor(eigenSubTypes) )
    defaultIfError = "EigenCG"
except:
    print("WARNING! Error loading Eigen linear solver using CG as default")
    defaultIfError = "CG"

###############################################################################################################################
class LinearProblem(BOO):
    def __init__(self):
        super().__init__()
        self.realsolver = None
        self.SetAlgo(defaultIfError)

    def SetTolerance(self,tol):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        self.realsolver.SetTolerance(tol)

    def HasConstraints(self):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        return self.realsolver.HasConstraints()

    def GetNumberOfDofs(self):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        return self.realsolver.GetNumberOfDofs()

    def ComputeProjector(self,mesh,fields):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        self.realsolver.ComputeProjector(mesh,fields)

    # you must set SetAlgo before setting the Op
    def SetOp(self, op):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        self.realsolver.SetOp(op)

    def SetAlgo(self, name, ops=None, withErrorIfNotFound=False):
        try:
            if self.realsolver is not None:
                constraints = self.realsolver.GetConstraints()
                self.realsolver = SolverFactory.Create(name,ops=ops)
                self.realsolver.SetConstraints(constraints)
            else:
                self.realsolver = SolverFactory.Create(name)
        except:
            if withErrorIfNotFound: #pragma: no cover
                raise
            else:
                print(f"Solver {name} unavailable, falling back to solver {defaultIfError} instead.")
                self.SetAlgo(defaultIfError)

    def Solve(self, rhs):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        return self.realsolver.Solve(rhs)

    @property
    def constraints(self):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        return self.realsolver.GetConstraints()

    @constraints.setter
    def constraints(self, constraints):
        if self.realsolver == None: #pragma: no cover
            raise(Exception("Please set the solver type first"))
        return self.realsolver.SetConstraints(constraints)

def CheckSolver(GUI,solver):
    print("Solver "+ str(solver))
    LS = LinearProblem ()
    if GUI:
        LS.SetGlobalDebugMode()
    LS.SetAlgo(solver)

    print("Number of dofs")
    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]]),dtype=PBasicFloatType))
    print("HasConstraints : ", LS.HasConstraints() )
    print("Number of dofs : ", LS.GetNumberOfDofs() )

    sol = LS.Solve(np.array([[1.],[2.]]))
    # second run
    sol = LS.Solve(np.array([[1.],[2.]]))

    if abs(sol[0] - 2.) >1e-15  :
        print(sol[0]-2) #pragma: no cover
        print("sol : ",sol) #pragma: no cover
        raise Exception() #pragma: no cover
    if abs(sol[1] - 2.) > 1e-15 : raise Exception()


def CheckSPQR(GUI):
    realsolver = LinearSolverEigen("SPQR")


    A = sps.csc_matrix(np.array([[0,0,0,1,1],
                                    [0,0,0,1,1],
                                    [0.51,0.5,0.5,0,10],
                                    [0.5,0.5,0.5,0,10],
                                    [0,1,0,0,5],
                                    ]),dtype=PBasicFloatType)

    def QRTest(A):
        realsolver.SetTolerance(1e-5)
        print(A.toarray())
        realsolver.SetOp(A)
        rank = realsolver.GetSPQRRank()
        print("Rank",rank )
        Q = realsolver.GetSPQR_Q()
        print("Q")
        print(Q.toarray())
        R = realsolver.GetSPQR_R()
        print("R")
        print(R.toarray())
        P = realsolver.GetSPQR_P()
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
    obj =  SolverFactory()
    solvers = GetAvailableSolvers()

    for s in solvers:
        CheckSolver(GUI,s)

    CheckSPQR(GUI)

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
