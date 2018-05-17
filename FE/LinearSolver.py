# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spslin

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from BasicTools.Helpers.TextFormatHelper import TFormat as TF

# you must set SetAlgo before setting the Op

class LinearProblem(BOO):
    def __init__(self):
        super(LinearProblem,self).__init__()
        self.op = None   # the operator to solve
        self.solver = None
        self.type = "CG"
        self.u = None
        self.tol = 1.e-6

    def SetOp(self, op):
        self.PrintVerbose('In SetOp')
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
        elif self.type == "cholesky":
            from sksparse.cholmod import cholesky
            self.solver = cholesky(op)


        self.PrintDebug('In SetOp Done')

    def SetAlgo(self, algoType):
        if algoType in  ["Direct" ,"CG", "lstsq","cholesky" ] :
            self.type = algoType
        else:
            self.Print(TF.InRed("Error : ") + "Type not allowed ("+algoType+")"  ) #pragma: no cover

    def Solve(self, rhs):
        self.PrintDebug('In LinearProblem Solve')
        if self.type == "Direct":
            self.u =  self.solver(rhs)
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

            if res[1] > 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Convergence to tolerance not achieved"))) #pragma: no cover
            if res[1] < 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Illegal input or breakdown"))) #pragma: no cover
        elif self.type == "lstsq":
            self.u = np.linalg.lstsq(self.op, rhs)[0]
        elif self.type == "cholesky":
            self.u = self.solver(rhs)

        self.PrintDebug("Done Linear solver")
            #norm_rhs = np.linalg.norm(rhs)
            #norm = np.linalg.norm(self.op.dot(self.u)-rhs)/norm_rhs
            #self.Print(TF.InYellowBackGround(TF.InRed("norm(error)/norm(rhs) : " + str(norm) )))
            #self.Print(TF.InYellowBackGround(TF.InRed("norm(rhs)  : " + str(norm_rhs))))

        return self.u

def CheckIntegrity():

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


    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
