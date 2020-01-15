# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

class ConstraintsHolder(BOO):
    def __init__(self,nbdof=0):
        super(ConstraintsHolder,self).__init__()
        self.SetNumberOfDofs(nbdof)
        self.constraints = []
        self.Reset()

        self.SetConstraintsMethod("Projection")

#----------------------- for population of the class --------------------------

    def Reset(self):
        self.numberOfEquations=0
        self.rows = []
        self.cols = []
        self.vals = []

    def AddConstraint(self,constraints):
        self.constraints.append(constraints)

    def ComputeConstraintsEquations(self, mesh, unkownFields):
        self.Reset()
        for cons in self.constraints:
            cons.GenerateEquations(mesh,unkownFields,self)
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
        for j,val in zip(range(len(vals)),vals):
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
        ej = np.array(ej)
        ev = np.array(ev)
        s = np.unique(ei)

        for i in s:
            mask = np.equal(ei,i)
            self.AddEquationSparse(ej[mask],ev[mask],0)

    def NextEquation(self):
        self.numberOfEquations += 1

    def Compact(self):
        r = max(self.rows)+1
        c = max(self.cols)+1

        res = coo_matrix((self.vals, (self.rows, self.cols)), shape=((r, c)), copy=False )

        self.rows = res.row
        self.cols = res.col
        self.vals = res.data

    def SetOp(self,op,rhs=None):
        self.op = sparse.csr_matrix(op)
        self.rhs = rhs
        self.nbdof = self.op.shape[1]

#----------------------- To recover the system of equations -------------------

    def ToSparse(self):
        from scipy.sparse import coo_matrix
        mask = np.zeros(self.nbdof+1,dtype=float)

        usedDofs = np.sort(np.unique(self.cols))
        if len(usedDofs) == 0:
            return coo_matrix(([], ([], [])), shape=((0, 0 )), copy= True ), usedDofs

        if usedDofs[-1] != self.nbdof:
        	#we add the independent term
            usedDofs= np.append(usedDofs, [self.nbdof])

        nbUsedDofs = len(usedDofs)

        mask[usedDofs] = range(nbUsedDofs)
        cols = mask[self.cols]
        res = coo_matrix((self.vals, (self.rows, cols)), shape=((self.numberOfEquations, nbUsedDofs )), copy= True )

        return res, usedDofs

    def GetNumberOfConstraints(self):
        return self.numberOfEquations

    def GetNumberOfDofsOnOriginalSystem(self):
        return self.nbdof

    def _Factorise(self):
        #algorithm form https://rosettacode.org/wiki/Reduced_row_echelon_form#Python
        # with some modification to treat sparse matrices in dense mode ( the ifs)

        M, usedDofs = self.ToSparse()
        M = M.toarray()
        def expand(M,r):
           res = sparse.coo_matrix(M[0:r,:])
           return (res, usedDofs[0:-1])

        rowCount = self.GetNumberOfConstraints()
        columnCount = len(usedDofs)

        lead = 0
        for r in range(rowCount):
            if lead >= columnCount:
                return expand(M,r)
            i = r
            rollcpt = columnCount - 2
            while M[i,lead] == 0:
                i += 1
                # if no non zero in the col
                # we dont do the permutation of the last col
                if i == rowCount and rollcpt > lead  :
                    i = r
                    M[:,[lead,rollcpt]] = M[:,[rollcpt,lead]]
                    usedDofs[[lead,rollcpt]] = usedDofs[[rollcpt,lead]]
                    rollcpt -= 1
                    continue

                if i == rowCount:
                    i = r
                    lead += 1
                    if columnCount == lead:
                        return expand(M,r)
            if i != r:
                M[[i,r],:] = M[[r,i],:]

            lv = M[r,lead]
            if lv != 1:
                M[r,:] = M[r,:]/lv

            for i in range(rowCount):
                if i != r:
                    lv = M[i,lead]
                    if lv != 0:
                        M[i,:] -=  M[r,:]*lv
            lead += 1
        return expand(M,r+1)
#-----------------------  External API ------------------

    def SetConstraintsMethod(self,method):

        if method.lower() == "lagrangemultipliers":
            self.method = Lagrange()
        elif method.lower() == "projection":
            self.method = Projection()
        elif method.lower() == "penalisation":
            self.method = Penalisation()
        else:
            raise(Exception("Error"))

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
        return self.method.RestoreSolution(arg)

    def RestrictSolution(self,arg):
        return self.method.RestrictSolution(arg)

    def __str__(self):
        res = "Constraints Holder: \n"
        res = "   Number of contraints: "+str(self.GetNumberOfConstraints()) +"\n"
        res += str(self.method)
        return res

#--------------------------------------------------------
class Penalisation(BOO):
    def __init__(self):
        super(Penalisation, self).__init__()
        self.factor = 1e8
        self.maxdiag =1.

    def UpdateSystem(self,CH):
        self.nbdof = CH.GetNumberOfDofsOnOriginalSystem()

        mat, self.mattoglobal = CH._Factorise()

        self.mat = mat.tocsr()[:,0:-1].tocoo()
        nbpot = len(self.mat.data)

        op =  self.mat.T.dot(self.mat).tocoo()

        data = op.data
        rows = self.mattoglobal[op.row]
        cols = self.mattoglobal[op.col]
        self.penalOp =  sparse.coo_matrix((data,(rows,cols)),shape=(self.nbdof,self.nbdof)).tocsr()

        self.penalRhs = np.zeros(self.nbdof)
        rhs = self.mat.T.dot(mat.tocsr()[:,-1]).toarray().ravel()
        self.penalRhs[self.mattoglobal] = rhs

    def GetCOp(self,op):
        self.maxdiag = max(op.diagonal())
        return op + self.maxdiag*self.factor * self.penalOp

    def GetCRhs(self, op,rhs):
        return rhs + self.maxdiag*self.factor*self.penalRhs

    def GetNumberOfDofs(self):
    	return self.nbdof

    def matvec(self,op,arg):
        return op.dot(arg) +self.maxdiag*self.penalOp.dot(self.factor*arg)

    def RestoreSolution(self,arg):
        return arg

    def RestrictSolution(self,arg):
        return arg

class Lagrange(BOO):
    def __init__(self):
        super(Lagrange,self).__init__()

    def __str__(self):
        res = "Lagrange method:\n"
        return res

    def UpdateSystem(self,CH):
        mat, self.mattoglobal = CH._Factorise()
        self.mat = mat.tocsr()[:,0:-1].tocoo()
        self.rhs = mat.tocsr()[:,-1].toarray()
        self.nbdof = CH.GetNumberOfDofsOnOriginalSystem()

    def GetCOp(self,op):
    	nbdof = op.shape[0]
    	nbcons = self.mat.shape[0]
    	op= op.tocoo()
    	data = np.copy(op.data)
    	rows =np.copy(op.row)
    	cols =np.copy(op.col)

    	data = np.hstack((data,self.mat.data,self.mat.data))
    	rows = np.hstack((rows, self.mat.row+nbdof,self.mattoglobal[self.mat.col]))
    	cols = np.hstack((cols, self.mattoglobal[self.mat.col], self.mat.row+nbdof))

    	return sparse.coo_matrix((data,(rows,cols)),shape=(nbdof+nbcons,nbdof+nbcons)).tocsr()

    def GetCRhs(self, op,rhs):
        return np.hstack((rhs.ravel(),self.rhs.ravel()))

    def GetNumberOfDofs(self):
    	nbcons = self.mat.shape[0]
    	return self.nbdof+nbcons

    def matvec(self,op,arg):
        u = np.zeros(self.GetNumberOfDofs())
        u[0:self.nbdof] += op.dot(arg[0:self.nbdof])
        u[self.nbdof:]  += self.mat.dot(arg[self.mattoglobal])
        u[self.mattoglobal] += self.mat.T.dot(arg[self.nbdof:])
        return u

    def RestoreSolution(self,arg):
        return arg[0:self.nbdof]

    def RestrictSolution(self,arg):
        res = np.zeros(self.GetNumberOfDofs())
        res[0:self.nbdof] = arg
        return arg

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
        return len(self.masters)

    def rmatvec(self,op,arg):
        return np.dot(arg,self.X.T.dot(op.dot(self.X ) ) )

    def matvec(self,op,arg):
        return self.X.T.dot(op.dot(self.X.dot(arg) ) )

    def GetCOp(self,op):
        return self.X.T.dot(op.dot(self.X ) )

    def GetCRhs(self, op,rhs):
        return self.X.T.dot( (rhs - op.dot(self.D) )   )

    def RestoreSolution(self,arg):
        return self.X.dot(arg) + self.D

    def RestrictSolution(self,arg):
        return self.X.T.dot(arg)

    def UpdateSystem(self,CH):
        self.PrintVerbose('Compute Factorisation')

        mat,mattoglobal = CH._Factorise()
        mat = mat.tocsr()

        self.PrintVerbose('Compute Factorisation done')

        PP = []

        for i in range(mat.shape[0]):
            for j in range(i,len(mattoglobal)):
                if mat[i,j] != 0:
                    PP.append(j)
                    break

        slavesInMat = np.array(PP)

        np.set_printoptions(threshold=np.inf,linewidth=np.inf)

        slaves = mattoglobal[PP]
        self.PrintVerbose('Creation of the projector')
        nbdof = CH.GetNumberOfDofsOnOriginalSystem()
        masterMask = np.ones((nbdof+1),dtype=bool )

        masterMatMask = np.ones(mat.shape[1],dtype=bool )
        self.slaves = np.array(slaves,dtype=int)

        masterMask[self.slaves] = False
        masterMask[-1] = False
        masterMatMask[PP] = False

        self.masters = np.where(masterMask)[0]

        nbMasterDofs = nbdof-len(slaves)
        self.PrintVerbose(nbMasterDofs)
        self.PrintDebug(mat.shape)

        submat = mat[:,masterMatMask[:-1]].tocoo()

        subTT = mattoglobal[masterMatMask[:-1]]

        data = np.hstack((np.ones(nbMasterDofs), -submat.data ))
        row = np.hstack((self.masters, mattoglobal[slavesInMat[submat.row]]   ))

        TT = np.empty(nbdof,dtype=int)
        TT[self.masters] = range(len(self.masters))
        col = np.hstack((np.arange(nbMasterDofs), TT[subTT[submat.col]] ))

        P = np.hstack( (self.masters,self.slaves) )
        self.X = sparse.coo_matrix( (data,(row,col )), shape=(nbdof,nbMasterDofs)   )

        self.D = np.zeros( (nbdof ), dtype=float )
        slaveValues = np.squeeze(mat[:,-1].toarray())
        self.D[self.slaves] = slaveValues

def CheckIntegrity(GUI=False):
    typeToCheck = ["lagrangemultipliers", "penalisation","projection",]
    res = True

    for ttc in typeToCheck:
        res = CheckIntegrityTTC(ttc,GUI=GUI)
        if res.lower() != 'ok':
               return res
    return 'ok'

def CheckIntegrityTTC(ttc,GUI=False):
    print('-------------------  '+ttc+'  -----------------')
    CH = ConstraintsHolder()
    CH.SetConstraintsMethod(ttc)

    CH.SetGlobalDebugMode()
    CH.SetNumberOfDofs(4)

    refsol = np.array([0,1,2,3],dtype=float )

    sys = 1000*np.array([[ 1,-1, 0, 0],
                         [-1, 1, 0, 0],
                         [ 0, 0, 1,-1],
                         [ 0, 0,-1, 1]],dtype=float )

    rhs = np.array([0,0,0,0],dtype=float)

    # Dirichlet zero left
    CH.AddEquation([1,0,0,0],0)
    CH.AddEquation([1.,0,0,0],0)
    #kinematic relation
    CH.AddEquation([0.,-1,1.,0.],1)
    #to test the elimination of redundant equations
    CH.AddEquation([0.,-2,2.,0.],2)


    #Dirichlet right side
    CH.AddEquation([0,0,0,1.],3.)
    CH.AddEquation([0,0,0,1.],3.)
    CH.Compact()
    print(CH)
    CH.UpdateCSystem()
    print(CH)

    CH.SetOp( sys , rhs )
    print("op \n", CH.op.toarray())
    print("\nROp  : \n",CH.GetCOp().toarray() )
    print("\nRRhs : \n",CH.GetCRhs() )

    from scipy.sparse import linalg

    op = linalg.LinearOperator((CH.GetNumberOfDofsOnCSystem(),CH.GetNumberOfDofsOnCSystem()),matvec=CH.matvec,rmatvec=CH.rmatvec,dtype=float)
    sol = linalg.gmres(op,CH.GetCRhs())[0]

    print("sol ", sol )
    refA = np.linalg.norm(CH.GetCOp().dot(sol)- CH.GetCRhs())/np.linalg.norm(CH.GetCRhs())
    print('verification reduced A', refA)
    refB = np.linalg.norm(CH.matvec(sol)- CH.GetCRhs())/np.linalg.norm(CH.GetCRhs())

    print('verification reduced B', refB)

    ifsol = CH.RestoreSolution(sol)
    print("refsol solution", refsol )
    print("--------------------------")
    print("Final iterative solution", ifsol )
    print("error",np.linalg.norm(ifsol - refsol ))
    sol = sparse.linalg.spsolve(CH.GetCOp(),CH.GetCRhs() )

    K = CH.GetCOp().toarray()
    sol = np.linalg.inv(K).dot(CH.GetCRhs())
    dfsol = CH.RestoreSolution( sol )
    print("Final direct solution", dfsol )
    error = np.linalg.norm(dfsol - refsol )
    print("error", error)
    if error/np.linalg.norm(refsol) > 1e-6:
        return 'not ok'
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))# pragma: no cover

