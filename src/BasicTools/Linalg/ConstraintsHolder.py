# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Dict, Any

import numpy as np
from scipy.sparse import linalg as splinalg

from scipy import sparse
import scipy.linalg as sp_linalg
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

methodFactory : Dict[str,Any] = {}

class ConstraintsHolder(BOO):
    """
    Constrained linear problem: Class to store and to enforce constrains to a
    quadratic system on the form:

        minimize 1/2 u.K.u - u.f
        with the constrain :
        A.u = b

    This class can store, and manipulate the system Au=b and enforced on the Ku=f
    by different methods (penalisation, subtitution, lagrange multiplier)

    the first set of function is used to fill this class with the system Au=b:

        NextEquation must be called after the call of a set of : AddFactor, AddConstant
        One per constraint

        AddEquationSparse, AddEquation, AddEquationsFromIJV to add a full
        constraint at once

        Compact will clean (zeros,and duplicate) entries


    """

    def __init__(self,nbdof=0):
        super(ConstraintsHolder,self).__init__()
        self.SetNumberOfDofs(nbdof)
        self.constraints = []
        self.Reset()

        self.SetConstraintsMethod("Projection")
        self.tol = 1e-8 # tolerance for the detection of redundant constraints

#----------------------- for population of the class --------------------------

    def Reset(self):
        self.numberOfEquations=0
        self.rows = []
        self.cols = []
        self.vals = []
        self.op = None
        self.rhs =[]

    def AddConstraint(self,constraints):
        self.constraints.append(constraints)

    def ComputeConstraintsEquations(self, mesh, unkownFields):
        self.Reset()
        for cons in self.constraints:
            cons.GenerateEquations(mesh,unkownFields,self)
        self.PrintVerbose(self.nbdof)
        self.PrintVerbose(len(self.constraints) )
        self.PrintVerbose(self.numberOfEquations )

    def SetNumberOfDofs(self,nbdof):
        self.nbdof = nbdof

    def AddFactor(self,ddl,factor):
        if(factor == 0 ): return
        self.rows.append(self.numberOfEquations)
        self.cols.append(ddl)
        self.vals.append(factor)

    def AddConstant(self,constant):
        if constant == 0 : return

        self.rows.append(self.numberOfEquations)
        self.cols.append(self.nbdof)
        self.vals.append(constant)

    def AddEquation(self,vals,const=0):
        for j,val in enumerate(vals):
            self.AddFactor(j,val)
        if const :
            self.AddConstant(const)
        self.NextEquation()

    def AddEquationSparse(self,index,vals,const):
        for j,v in zip(index,vals):
            self.AddFactor(j,v)
        if const :
            self.AddConstant(const)
        self.NextEquation()

    def AddEquationsFromIJV(self,ei,ej,ev):
        ei = np.require(ei, dtype=PBasicIndexType)
        ej = np.require(ej, dtype=PBasicIndexType)
        ev = np.require(ev, dtype=PBasicFloatType)
        s = np.unique(ei)

        for i in s:
            mask = np.equal(ei,i)
            self.AddEquationSparse(ej[mask],ev[mask],0)

    def NextEquation(self):
        self.numberOfEquations += 1

    def Compact(self):
        r = self.numberOfEquations
        c = self.nbdof+1

        res = coo_matrix((self.vals, (self.rows, self.cols)), shape=((r, c)), copy=False, dtype=PBasicFloatType )

        self.rows = res.row
        self.cols = res.col
        self.vals = res.data


    def AddEquationsFromOperatorAAndb(self, A, b = None):
        self.rows.extend(A.row+self.numberOfEquations)
        self.cols.extend(A.col)
        self.vals.extend(A.data)

        if b is not None:
            self.rows.extend(np.arange(A.shape[0],dtype=PBasicIndexType)+self.numberOfEquations )
            self.cols.extend(np.ones(A.shape[0],dtype=PBasicIndexType)*self.nbdof,dtype=PBasicIndexType)
            self.vals.extend(b)
            raise

        self.numberOfEquations += A.shape[0]

    def SetOp(self,op,rhs=None):
        self.op = sparse.csr_matrix(op)
        self.rhs = rhs
        self.nbdof = self.op.shape[1]
#----------------------- To recover the system of equations -------------------

    def ToSparse(self):
        """
        Function to convert to the constraints system Ax=B to a sparse matrix

        this function returns:
            res : a sparse matrix of [A'|B] , A' has only the non zero columns of A
            usedDofs : the indices of the columns present in A'
        """
        return self.CleanEmptyColumns(self.ToSparseFull())

    def ToSparseFull(self):
        """
        Function to convert to the constraints system Ax=B to a sparse matrix

        this function returns:
            res : a sparse matrix of [A|B] ,
        """
        mat = coo_matrix((self.vals, (self.rows, self.cols)), shape=((self.numberOfEquations, self.nbdof+1 )), copy= False )
        return mat

    def CleanEmptyColumns(self,mat):
        """
        Function to clean empty columns for a matrix,

        this function returns:
            res : a sparse matrix with only non zero columns of mat (the last
            column is always present, event if it is full of zeros)
            usedDofs : the indices of the columns prensent in res
        """

        mat = mat.tocoo()
        nbdof = mat.shape[1]-1

        usedDofs = np.unique(mat.col).astype(dtype=PBasicIndexType)
        if len(usedDofs) == 0:
            return coo_matrix(([], ([], [])), shape=((0, 0 )), copy= True ), usedDofs

        if usedDofs[-1] != nbdof :
            #we add the independent term in the case is all zeros (last column = 0)
            usedDofs= np.append(usedDofs, [nbdof])

        nbUsedDofs = len(usedDofs)
        # map fron global to the reduced matrix
        mapGToR = np.zeros(mat.shape[1],dtype=PBasicIndexType)
        mapGToR[usedDofs] = range(nbUsedDofs)
        col = mapGToR[mat.col]
        res = coo_matrix((mat.data, (mat.row, col)), shape=((mat.shape[0], nbUsedDofs )), copy= True )

        return res, usedDofs

    def GetNumberOfConstraints(self):
        return self.numberOfEquations

    def GetNumberOfDofsOnOriginalSystem(self):
        return self.nbdof

    def GetCleanConstraintBase(self, purePython=False):
        """
        This function compute the QR decompostion of the augmentd system [A|b]
        Il will use the self.tol value to detect redundants constrants (rank revealing)

        this functions returns:
            res : a parse matrix q' with only the relevant part of the QR
            decomposition.
            Q' containt only the non zero columns of [A|b] and always the last
            column (the b)
            usedDofs : the indices of the columns present in Q

        option:
            purePython: to force scipy for computing the QR decomposition
        """
        self.PrintDebug("GetCleanConstraintBase ")

        if purePython == False:
            try:
                import BasicTools.Linalg.NativeEigenSolver as NativeEigenSolver
                LS = NativeEigenSolver.CEigenSolvers()
            except: # pragma: no cover
                purePython = True
                print("Warning!!! Eigen SPQR decomposition not available using slower algorithm (scipy.linalg.qr)")

        if purePython == True:
            from scipy.linalg  import qr
            class WrapedScipyQr():
                def __init__(self):
                    self.tol = 1.e-6
                def SetTolerance(self,tol):
                    self.tol = tol
                def SetSolverType(self,st):
                    if st != "SPQR": # pragma: no cover
                        raise(Exception("SolverType '"+str(st)+"' not allowed "))
                def SetOp(self,op):
                    Q,r,P= qr(op.toarray(),mode="full", pivoting=True,check_finite=False)
                    diag = r.diagonal()
                    diag = diag/diag[0]
                    mask = (abs(diag)<self.tol)
                    self.rank = mask.argmax()
                    if self.rank == 0:
                        self.rank = op.shape[1]
                    self.Q = sparse.csr_matrix(Q)
                def GetSPQRRank(self):
                    return self.rank
                def GetSPQR_Q(self):
                    return self.Q
            LS = WrapedScipyQr()
        LS.SetTolerance(self.tol)
        LS.SetSolverType("SPQR")

        M, usedDofs = self.ToSparse()

        Mcsc = M.tocsc()
        am = abs(Mcsc[:,:-1])

        self.PrintDebug("compute of adjacency_matrix ")
        adjacency_matrix = am.dot(am.T)
        self.PrintDebug("compute of adjacency_matrix done")
        del am

        Mcsr = M.tocsr()

        self.PrintDebug("connected_components ")
        ncc, cc = connected_components(adjacency_matrix,directed=False)
        self.PrintDebug("connected_components Done")

        res_vals = []
        res_rows = []
        res_cols = []
        rankcpt = 0
        self.PrintDebug(f"treating sub matrices ({ncc})" )

        # shortcut to check if we have independent 1 dof equations
        # this means all dofs are slaves
        if ncc == M.shape[0] :
            self.PrintVerbose(f"using shortcut to treat KCs" )
            res_vals = M.data
            res_rows = M.row
            res_cols = M.col
            rankcpt = M.shape[0]
            slaves = np.arange(M.shape[1]-1,dtype=PBasicIndexType)
        else:
            atonce = np.zeros( M.shape[0],dtype=bool)
            slaves = []
            for incc in range(ncc):
                mask = (cc == incc)
                nl = np.count_nonzero(mask)
                if nl == 1:
                    # we have one independet equation
                    np.logical_or(atonce,mask,out=atonce)
                    continue
                # extraction sof the equation to treat (rows) and filter the matrix
                # to have only the used dofs (cols)
                subsubM, subusedDofs =  self.CleanEmptyColumns(Mcsr[mask,:])
                LS.SetOp(subsubM.T)
                rank = LS.GetSPQRRank()
                slaves.extend(subusedDofs[0:rank])
                Q = LS.GetSPQR_Q()
                realQ = Q.tocsr()[:,0:rank].tocoo().T

                res_vals.extend(realQ.data)
                res_rows.extend(realQ.row+rankcpt)
                res_cols.extend(subusedDofs[realQ.col])
                rankcpt +=rank
            # we treat all the single equation constraint at once
            nb_single_constraint = np.count_nonzero(atonce)
            if nb_single_constraint:
                subM = Mcsr[atonce,:].tocoo()
                res_vals.extend(subM.data)
                res_rows.extend(subM.row+rankcpt)
                res_cols.extend(subM.col)
                rankcpt += nb_single_constraint
                cols = np.unique(subM.col)
                if cols[-1] == M.shape[1]-1:
                    slaves.extend(cols[0:-1] )
                else:
                    slaves.extend(cols )

        # normalisation
        res = coo_matrix((res_vals,(res_rows,res_cols)), shape=(rankcpt,len(usedDofs) ) )
        rescsr = res.tocsr()
        norm = 1/splinalg.norm(rescsr[:,:-1],axis=1)
        res.data *= norm[res.row]

        self.PrintDebug("treating sub matrices Done ")

        return res, usedDofs, slaves

#-----------------------  External API ------------------

    def SetConstraintsMethod(self,method):
        MethodClass =  methodFactory.get(method.lower(),None)
        if MethodClass is None:# pragma: no cover
            raise(Exception(f"Method '{method}' to avaible: possible options are {list(methodFactory.keys())}"))

        self.method = MethodClass()

    def UpdateCSystem(self):
        self.method.UpdateSystem(self)

    def GetNumberOfDofsOnCSystem(self):
        return self.method.GetNumberOfDofs()

    def GetCOp(self):
        return self.method.GetCOp(self.op)

    def matvec(self,arg):
        return self.method.matvec(self.op,arg)

    def rmatvec(self,arg):
        return self.method.rmatvec(self.op,arg)

    def GetCRhs(self,rhs=None):
        if rhs is not None:
            self.rhs = rhs
        return self.method.GetCRhs(self.op,self.rhs)

    def RestoreSolution(self,arg):
        return self.method.RestoreSolution(self.op, self.rhs, arg)

    def RestrictSolution(self,arg):
        return self.method.RestrictSolution(self.op, self.rhs, arg)


    def __str__(self):
        res = "Constraints Holder: \n"
        res = "   Number of contraints: "+str(self.GetNumberOfConstraints()) +"\n"
        res += str(self.method)
        return res

def ExpandMatrix(op,mattoglobal,nbdofs,  treatrows=True, treatcols=True):
    op = op.tocoo()
    data = op.data
    if treatrows :
        rows = mattoglobal[op.row]
        rs = nbdofs[0]
    else:
        rows = op.row
        rs = op.shape[0]

    if treatcols :
        cols = mattoglobal[op.col]
        cs = nbdofs[1]
    else:
        cols = op.col
        cs = op.shape[1]

    return  sparse.coo_matrix((data,(rows,cols)),shape=(rs,cs)).tocsr()

#--------------------------------------------------------
class Ainsworth(BOO):
    """
    Method from https://doi.org/10.1016/S0045-7825(01)00236-5
    Because the CleanConstraint(C) matrix is orthonormal then C.dot(C.T) = I
    this make the P Q R expression simpler

    """
    def __init__(self):
        super(Ainsworth, self).__init__()
        self.P = None
        self.Q = None
        self.R = None
        self.g = None
        self.mattoglobal = None

    def UpdateSystem(self,CH):
        self.nbdof = CH.GetNumberOfDofsOnOriginalSystem()
        mat, self.mattoglobal, slaves = CH.GetCleanConstraintBase()
        matcsc = mat.tocsc()
        C = matcsc[:,0:-1]
        nbdofs = (self.nbdof,self.nbdof)
        self.g = matcsc[:,-1]
        CCT_1 = splinalg.inv(C.dot(C.T))
        r = C.T.dot(CCT_1)
        self.R = ExpandMatrix(r,self.mattoglobal,nbdofs,treatcols=False).tocsc()
        p = r.dot(C)
        self.P = ExpandMatrix(p,self.mattoglobal,nbdofs).tocsc()
        self.Q = (sparse.eye(self.nbdof) - self.P).tocsc()

    def GetCOp(self,op):
        return self.P + self.Q.T.dot(op.dot(self.Q))

    def GetCRhs(self, op,rhs):

        a = self.R.dot(self.g).toarray()[:,0]
        b = self.Q.T.dot(rhs - op.dot(self.R.dot(self.g).toarray()[:,0]))
        return (a +b).flatten()

    def GetNumberOfDofs(self):
        return self.nbdof

    def matvec(self,op,arg):
        return self.P.dot(arg) + self.Q.T.dot(op.dot(self.Q.dot(arg)))

    def RestoreSolution(self,op,rhs,arg):
        return arg

    def RestrictSolution(self,op, rhs, arg):
        return arg

methodFactory["Ainsworth".lower()] = Ainsworth

class Penalisation(BOO):
    def __init__(self):
        super(Penalisation, self).__init__()
        self.factor = 1e8
        self.maxdiag =1.

    def UpdateSystem(self,CH):
        self.nbdof = CH.GetNumberOfDofsOnOriginalSystem()

        mat, self.mattoglobal, slaves = CH.GetCleanConstraintBase()

        matcsc = mat.tocsc()
        A = matcsc[:,0:-1]
        b = matcsc[:,-1]

        op =  A.T.dot(A).tocoo()  #  A^T A

        data = op.data
        rows = self.mattoglobal[op.row]
        cols = self.mattoglobal[op.col]
        self.penalOp =  sparse.coo_matrix((data,(rows,cols)),shape=(self.nbdof,self.nbdof)).tocsr()

        self.penalRhs = np.zeros(self.nbdof)
        rhs = A.T.dot(b).toarray().ravel()
        self.penalRhs[self.mattoglobal[:-1]] = rhs

    def GetCOp(self,op):
        self.maxdiag = max(op.diagonal())
        return op + self.maxdiag*self.factor * self.penalOp

    def GetCRhs(self, op,rhs):
        return rhs + self.maxdiag*self.factor*self.penalRhs

    def GetNumberOfDofs(self):
        return self.nbdof

    def matvec(self,op,arg):
        return op.dot(arg) +self.factor*self.maxdiag*self.penalOp.dot(arg)

    def RestoreSolution(self,op,rhs,arg):
        return arg

    def RestrictSolution(self,op, rhs, arg):
        return arg

methodFactory["Penalisation".lower()] = Penalisation

class Lagrange(BOO):
    def __init__(self):
        super(Lagrange,self).__init__()

    def __str__(self):
        res = "Lagrange method:\n"
        return res

    def UpdateSystem(self,CH):
        mat, self.mattoglobal, slaves = CH.GetCleanConstraintBase()
        self.mat = mat.tocsr()[:,0:-1].tocoo()
        self.rhs = mat.tocsr()[:,-1].toarray().ravel()
        self.nbdof = CH.GetNumberOfDofsOnOriginalSystem()

    def GetCOp(self,op):
        nbdof = op.shape[0]
        nbcons = self.mat.shape[0]
        op= op.tocoo()

        data = np.hstack((op.data,self.mat.data,self.mat.data))
        rows = np.hstack((op.row, self.mat.row+nbdof,self.mattoglobal[self.mat.col]))
        cols = np.hstack((op.col, self.mattoglobal[self.mat.col], self.mat.row+nbdof))

        return sparse.coo_matrix((data,(rows,cols)),shape=(nbdof+nbcons,nbdof+nbcons))

    def GetCRhs(self, op,rhs):
        return np.hstack((rhs.ravel(),self.rhs.ravel()))

    def GetNumberOfDofs(self):
        nbcons = self.mat.shape[0]
        return self.nbdof+nbcons

    def matvec(self,op,arg):
        u = np.zeros(self.GetNumberOfDofs())
        u[0:self.nbdof] += op.dot(arg[0:self.nbdof])
        u[self.nbdof:]  += self.mat.dot(arg[self.mattoglobal[:-1]])
        u[self.mattoglobal[:-1]] += self.mat.T.dot(arg[self.nbdof:])
        return u

    def RestoreSolution(self,op,rhs,arg):
        return arg[0:self.nbdof]

    def RestrictSolution(self, op, rhs, arg):
        res = np.zeros(self.GetNumberOfDofs())
        res[0:self.nbdof] = arg
        res[self.nbdof:] =   self.mat.dot( (rhs - op.dot(arg))[self.mattoglobal[:-1]] )
        return res

methodFactory["Lagrange".lower()] = Lagrange

class Projection(BOO):
    def __init__(self):
        super(Projection,self).__init__()
        self.slaves = None
        self.masters = None
        self.X = None
        self.D = None

    def __str__(self):

        res = "Projection method:\n"
        res += "Slaves " + str(self.slaves) + "\n"
        res += "Masters " + str(self.masters) + "\n"
        if self.X is not None:
            res += "X " + str(self.X.toarray()) + "\n"
        else:
            res += "X None\n"

        res += "D " + str(self.D) + "\n"
        return res

    def GetNumberOfDofs(self):
        return (self.nbMaster)

    def matvec(self,op,arg):
        return self.X.T.dot(op.dot(self.X.dot(arg) ) )

    def GetCOp(self,op):
        return self.X.T.dot(op.dot(self.X ) )   #  X^T*Op*X

    def GetCRhs(self, op,rhs):
        return self.X.T.dot( (rhs ) - op.dot(self.D))

    def RestoreSolution(self,op,rhs,arg):
        return  self.X.dot(arg) +self.D

    def RestrictSolution(self,op,rhs,arg):
        return arg[self.masters]

    def UpdateSystem(self,CH):
        self.PrintVerbose('Compute Factorisation')
        nbdof = CH.GetNumberOfDofsOnOriginalSystem()

        mat, self.mattoglobal, slaves = CH.GetCleanConstraintBase()
        slavesInMat = np.require(slaves,dtype=PBasicIndexType)
        slaveIds = self.mattoglobal[slavesInMat]

        mat = mat.tocsr()

        self.slaves = np.require(slaveIds,dtype=PBasicIndexType)

        masterMatMask = np.ones(mat.shape[1]-1,dtype=bool )
        masterMatMask[slavesInMat] = False
        masterInMat = np.where(masterMatMask)[0]

        masterMask = np.ones((nbdof),dtype=bool )
        masterMask[self.slaves] = False
        self.masters = np.where(masterMask)[0]

        self.nbMaster = len(self.masters)

        Ms = mat[:,slavesInMat]
        Mm = mat[:,masterInMat]

        self.Ms_1 = splinalg.splu(Ms)
        Ms_1Mm = sparse.coo_matrix(self.Ms_1.solve(Mm.toarray()))
        Xdata = np.hstack( (np.ones(self.nbMaster,dtype=PBasicFloatType), -Ms_1Mm.data  ) )
        Xrow  = np.hstack( ( self.masters, slaveIds[Ms_1Mm.row] ) )
        Xcol  = np.hstack( ( np.arange(self.nbMaster,dtype=PBasicIndexType), Ms_1Mm.col ) )

        self.X = sparse.coo_matrix( (Xdata,(Xrow,Xcol )), shape=(nbdof,self.nbMaster )   )

        self.D = np.zeros( (nbdof), dtype=PBasicFloatType )
        v = mat[:,-1].toarray()
        self.D[self.slaves] =   self.Ms_1.solve(v).ravel()

methodFactory["Projection".lower()] = Projection

def TestQR():

    CH = ConstraintsHolder()
    CH.SetGlobalDebugMode()
    CH.SetNumberOfDofs(2)
    CH.AddEquation([1,-1],2)
    CH.AddEquation([1.,0],3)
    CH.Compact()
    Ab = CH.ToSparse()[0].toarray()
    print(Ab)
    Cg= CH.GetCleanConstraintBase()[0].toarray()
    C = Cg[:,:-1]
    g = Cg[:,-1]
    print("C",C)
    print("g",g)

    print(C.T.dot(C))
    print(C.dot(C.T))
#TestQR()
#exit()

def CheckIntegrity(GUI=False):
    typeToCheck = list(methodFactory.keys())
    for ttc in typeToCheck:
        res = CheckIntegrityTTC(ttc,GUI=GUI)
        if res.lower() != 'ok':
               return res
    return 'ok'

def CheckIntegrityTTC(ttc,GUI=False):
    print('-------------------------------  '+ttc+'  ------------------------------')
    CH = ConstraintsHolder()
    CH.SetConstraintsMethod(ttc)

    if GUI:
        CH.SetGlobalDebugMode()
    CH.SetNumberOfDofs(4)

    # Reference solution
    refsol = np.array([0,1,2,3],dtype=PBasicFloatType )


    sys = 100*np.array([[ 1,-1, 0, 0],
                         [-1, 1, 0, 0],
                         [ 0, 0, 1,-1],
                         [ 0, 0,-1, 1]],dtype=PBasicFloatType )

    rhs = np.array([0,0,0,0],dtype=PBasicFloatType)
    # Dirichlet zero left
    CH.AddEquation([1,0,0,0],0)
    CH.AddEquationsFromIJV([0],[0],[1])
    #kinematic relation
    CH.AddEquation([0.,-1,1.,0.],1)
    #to test the elimination of redundant equations
    CH.AddEquation([0.,-2,2.,0.],2)
    #Dirichlet right side
    CH.AddEquation([0,0,0,1.],3.)
    CH.AddEquationSparse([3],[1.],3)
    CH.Compact()

    np.set_printoptions(threshold=np.inf,linewidth=np.inf)
    print("System to Solve ")
    print("Number Of constraints: ", CH.GetNumberOfConstraints())
    print("K ")
    print(sys)
    print("F ")
    print(rhs)
    print("Under Constraints ")
    Ab = CH.ToSparse()[0].toarray()
    A = Ab[:,:-1]
    b = Ab[:,-1]
    print("A ")
    print(A)
    print("b ")
    print(b)
    print("Reference solution ")
    print(refsol)

    CH.UpdateCSystem()


    Mv =  CH.GetCleanConstraintBase()[0].toarray()
    MvSlow = CH.GetCleanConstraintBase(purePython=True)[0].toarray()
    M = Mv[:,:-1]
    v = Mv[:,-1]
    print("Clean Contraints Base ")
    print("M ")
    print(M)
    print("v ")
    print(v)
    if len(v) != 3:
        raise Exception("Error in CheckIntegrity")

    CH.SetOp(sys,rhs)

    Kc = CH.GetCOp().tocsc()
    Fc = CH.GetCRhs()
    print("\nConstrained operator  : \n",Kc.toarray() )
    conditionNumber = np.linalg.cond(Kc.toarray(), p='fro')
    print("\nCondition number of Constrained operator  : \n", conditionNumber )
    print("\nConstrained rhs : \n",Fc )


    Uc = CH.RestrictSolution(refsol)
    print("\nRestricted solution (from refsol)  : \n", Uc)

    contrainedSolutionError =  np.linalg.norm(Kc.dot(Uc) - Fc)/np.linalg.norm(Fc)
    print("\nContraines solution error (|Kc(Uc)-Fc|/|Fc|)  : \n", contrainedSolutionError)
    if contrainedSolutionError > 1e-8:
        raise Exception("Error in CheckIntegrity")

    # try to solve the constrained system Using the matvec
    # we use the gmres to correctrly treat non defined matrices
    from scipy.sparse import linalg

    contrainedSysNbDofs = CH.GetNumberOfDofsOnCSystem()
    op = linalg.LinearOperator((contrainedSysNbDofs,contrainedSysNbDofs),matvec=CH.matvec,rmatvec=CH.rmatvec,dtype=PBasicFloatType)

    # the penalisation is not a good method with a iterative solver so we need
    # a simple precoditioning to help the solver
    diag = Kc.diagonal()
    diag[diag == 0] = 1.0
    M = sparse.dia_matrix((1./diag,0), shape=Kc.shape)

    sol = linalg.gmres(op,Fc,atol=0.,tol=0., M=M)[0]
    print("---------- Solution using gmres --------------  \n", sol )
    errorA = np.linalg.norm(Kc.dot(sol)- Fc)/np.linalg.norm(Fc)
    print('Error on the computed solution (using matricexs) |Kc(Uc)-Fc|/|Fc| \n', errorA)
    if errorA > 10e-8:
        raise Exception("Error in CheckIntegrity")
    errorB = np.linalg.norm(CH.matvec(sol)- Fc)/np.linalg.norm(Fc)
    print('Error on the computed solution (using matvec) | Op(Uc)-Fc |/|Fc| \n', errorB)
    if errorB > 10e-8:
        raise Exception("Error in CheckIntegrity")

    computed_solution = CH.RestoreSolution(sol)
    print("--------------------------")
    print("Reference solution \n", refsol )
    print("Computed solution \n", computed_solution )
    errorC = np.linalg.norm(computed_solution - refsol )/np.linalg.norm(refsol)
    print("relative error |ref_sol - computed_sol|/|ref_sol|",errorC)
    if errorC > 1e-8:
        raise Exception("Error in CheckIntegrity")

    print("---------- Solution using spsolve --------------  \n", sol )
    sol = sparse.linalg.spsolve(Kc,Fc )
    #sol = np.linalg.inv(K).dot(CH.GetCRhs())
    dfsol = CH.RestoreSolution( sol )
    print("Final direct solution", dfsol )
    errorD = np.linalg.norm(dfsol - refsol )/np.linalg.norm(refsol)
    print("Error on the solution |ref_sol - computed_sol|/|ref_sol| \n ", errorD)
    if errorD > 1e-8:
        raise Exception("Error in CheckIntegrity")
    print(CH)
    print('-------------------  '+ttc+' DONE  -----------------')

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))# pragma: no cover

