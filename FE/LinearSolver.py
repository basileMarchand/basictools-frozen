# -*- coding: utf-8 -*-

import scipy.sparse as sps
import scipy.sparse.linalg as spslin
import numpy as np
from OTTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

class LinearProblem(BOO):
    def __init__(self):
        super(LinearProblem,self).__init__()
        self.op = None   # the operator to solve
        self.solver = None
        self.type = "CG"
        self.u = None
        
    def SetOp(self, op):
        self.op = op
        self.PrintDebug('In SetOp')
        if self.type == "Direct":
            self.PrintDebug('Starting Factorisaton')
            self.solver = sps.linalg.factorized(op)
        elif self.type == "CG":
            pass
        self.PrintDebug('In SetOp Done')
    def Solve(self, rhs):
        if self.type == "Direct":
            self.u =  self.solver(rhs)
        elif self.type == "CG":
            M = sps.dia_matrix((1./self.op.diagonal(),0), shape=self.op.shape)
            if self.u is None:
                res = spslin.cg(self.op, rhs, M = M)#,tol=  1e-4) 
            else:
                res = spslin.cg(self.op, rhs, M = M, x0 = self.u)#,tol=  1e-4) 
            self.u = res[0]
        return self.u
        
def CheckIntegrity():
    
    LS = LinearProblem ()
    LS.SetOp(np.array([[0.5,0],[0,1]]))
    sol = LS.Solve(np.array([[1],[2]]))
    sol = LS.Solve(np.array([[1],[2]]))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()    
    
    LS.type = "Direct"
    LS.SetOp(np.array([[0.5,0],[0,1]]))
    sol = LS.Solve(np.array([[1],[2]]))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()    

        
    
    
    return "OK"     

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover