#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

class ConstraintsHolder(BOO):
    def __init__(self,nbdof=0):
        super(ConstraintsHolder,self).__init__()
        self.SetNumberOfDofs(nbdof)
        self.constraints = []
        self.Reset()

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

    def NextEquation(self):
        self.numberOfEquations += 1

    def ToSparse(self):
        from scipy.sparse import coo_matrix
        mask = np.zeros(self.nbdof+1,dtype=float)
        unique = np.sort(np.unique(self.cols))
        if unique[-1] != self.nbdof:
            unique = np.append(unique, [self.nbdof])
        mask[unique] = range(len(unique))
        cols = mask[self.cols]
        res = coo_matrix((self.vals, (self.rows, cols)), shape=((self.numberOfEquations, len(unique))), copy=False )
        return res,unique

    def _Factorise(self):
        #algorithm form https://rosettacode.org/wiki/Reduced_row_echelon_form#Python
        # with some modification to treat sparse matrices in dense mode ( the ifs)

        M,unique = self.ToSparse()
        M = M.toarray()

        def expand(M,r):
           res = sparse.coo_matrix(M[0:r,:])
           return (res, unique)

        rowCount = self.numberOfEquations
        columnCount = len(unique)

        lead = 0
        for r in range(rowCount):
            #print(r,lead)
            if lead >= columnCount:
                return expand(M,r)
            i = r
            while M[i,lead] == 0:
                i += 1
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
        return expand(M,r)
    def ComputeProjector(self):
        #mat = self.ToArray()
        #for l in mat.tolist():
        #    print("\n ", l)
        #print("\n-----------------------------")
        self.PrintVerbose('Compute Factorisation')
        mat,unique = self._Factorise()
        mat = mat.tocsr()
        self.PrintVerbose('Compute Factorisation done')

        self.PrintVerbose(self.numberOfEquations)
        self.PrintVerbose(mat.shape)
        self.PrintVerbose(len(unique))

        line = mat.shape[0]
        slaves = unique[range(line)]

        self.PrintVerbose('Creation of the projector')

        mask = np.ones((self.nbdof+1),dtype=bool )

        self.slaves = np.array(slaves,dtype=int)
        mask[self.slaves] = False
        mask[-1] = False

        self.masters = np.where(mask)[0]
        nbMasterDofs = self.nbdof-len(slaves)
        self.PrintVerbose(nbMasterDofs)

        self.PrintDebug(mat.shape)
        self.PrintDebug(line)

        submat = mat[:,line:-1].tocoo()

        data = np.hstack((np.ones(nbMasterDofs), submat.data ))
        row = np.hstack((np.arange(nbMasterDofs), submat.row +nbMasterDofs))
        col = np.hstack((np.arange(nbMasterDofs), unique[submat.col] ))

        P = np.hstack( (self.masters,self.slaves) )

        self.X = sparse.coo_matrix( (data,(P[row],col )), shape=(self.nbdof,nbMasterDofs)   )

        self.D = np.zeros( (self.nbdof ), dtype=float )
        self.slaveValues = np.squeeze(mat[0:line,-1].toarray())
        self.D[self.slaves] = self.slaveValues

    def GetNumberOfMasterDegrees(self):
        return len(self.masters)

    def SetOp(self,op,rhs=None):
        self.op = sparse.csr_matrix(op)
        self.rhs = rhs

    def rmatvec(self,arg):
        return np.dot(arg,self.X.T.dot(self.op.dot(self.X ) ) )

    def matvec(self,arg):
        return self.X.T.dot(self.op.dot(self.X.dot(arg) ) )

    def GetReducedOp(self):
        return self.X.T.dot(self.op.dot(self.X ) )

    def GetReducedRhs(self,rhs=None):
        if rhs is not None:
            self.rhs = rhs
        return self.X.T.dot( (self.rhs - self.op.dot(self.D) )   )

    def ExpandSolution(self,arg):
        return self.X.dot(arg) + self.D

def CheckIntegrity(GUI=False):

    # traction bar with 3 elements

    CH = ConstraintsHolder()
    #CH.SetGlobalDebugMode(True)

    CH.SetNumberOfDofs(4)
    refsol = np.array([0, 3, 6 , 9],dtype=float )+1

    sys = np.array([[ 1, -1.,  0, 0],
                    [-1,  2, -1, 0],
                    [ 0, -1,  2,-1],
                    [ 0,  0, -1, 1]],dtype=float )

    rhs = np.array([0,0,0,0],dtype=float)



    CH.AddEquation([1,0,0,0],1)
    CH.AddEquationSparse([3],[1],10)
    CH.AddEquation([0,0,0,1],10)

    # this class can handle redundant equations

    print("\nAumented Constaint matrix :")
    print(CH.ToSparse())

    CH.ComputeProjector()
    print("\nslaves ",CH.slaves)
    print("\nMasters ",CH.masters)
    print("\nX ", CH.X.toarray())
    print("\nD ",CH.D)

    print('verification original', sys.dot(refsol) - rhs )
    print('verification', sys.dot(CH.X.dot(refsol[CH.masters]) + CH.D)  -rhs )
    print("\nMasters ",CH.masters)

    CH.SetOp( sys , rhs )

    print("op \n", CH.op.toarray())
    print("\nROp  : ",CH.GetReducedOp().toarray() )
    print("\nRRhs : ",CH.GetReducedRhs() )

    from scipy.sparse import linalg


    op = linalg.LinearOperator((CH.GetNumberOfMasterDegrees(),CH.GetNumberOfMasterDegrees()),matvec=CH.matvec,rmatvec=CH.rmatvec,dtype=float)
    sol = linalg.gmres(op,CH.GetReducedRhs())[0]
    print("sol ", sol )
    print('verification reduced A', CH.GetReducedOp().dot(sol)- CH.GetReducedRhs())
    print('verification reduced B', CH.matvec(sol)- CH.GetReducedRhs())

    ifsol = CH.ExpandSolution(sol)
    print("Final iterative solution", ifsol )
    print("error",ifsol - refsol )
    print("--------------------------")
    dfsol = CH.ExpandSolution( sparse.linalg.spsolve(CH.GetReducedOp(),CH.GetReducedRhs() )  )
    print("Final direct solution", dfsol )
    print("error",dfsol - refsol )

    return "Ok"

if __name__ == '__main__':

    print(CheckIntegrity(GUI=True))# pragma: no cover
