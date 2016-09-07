# -*- coding: utf-8 -*-

import scipy.sparse as sps
import scipy.sparse.linalg as spslin
import numpy as np
from OTTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from OTTools.Helpers.TextFormatHelper import TFormat as TF

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
        self.PrintDebug('In SetOp')
        self.op = op
        if self.type == "Direct":
            self.PrintDebug('Starting Factorisaton')
            self.solver = sps.linalg.factorized(op)
        elif self.type == "CG":
            self.solver = None
            pass
        self.PrintDebug('In SetOp Done')
    def SetAlgo(self, algoType):
        if algoType == "Direct" or algoType == "CG":
            self.type = algoType
        else:
            self.Print(TF.InRed("Error : ") + "Type not allowed ("+algoType+")"  )

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
                self.Print(TF.InYellowBackGround(TF.InRed("Convergence to tolerance not achieved")))
            if res[1] < 0 :
                self.Print(TF.InYellowBackGround(TF.InRed("Illegal input or breakdown")))

        if self.InDebugMode():
            self.Print("Done Linear solver")
            #norm_rhs = np.linalg.norm(rhs)
            #norm = np.linalg.norm(self.op.dot(self.u)-rhs)/norm_rhs
            #self.Print(TF.InYellowBackGround(TF.InRed("norm(error)/norm(rhs) : " + str(norm) )))
            #self.Print(TF.InYellowBackGround(TF.InRed("norm(rhs)  : " + str(norm_rhs))))

        return self.u

def CheckIntegrity():

    LS = LinearProblem ()
    LS.SetGlobalDebugMode()

    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]])) )
    sol = LS.Solve(np.array([[1],[2]]))
    if abs(sol[0] - 2.) >1e-15  :
        print(sol[0]-2)
        raise Exception()
    if abs(sol[1] - 2.) > 1e-15 : raise Exception()

    LS.type = "Direct"

    LS.SetOp(sps.csc_matrix(np.array([[0.5,0],[0,1]])))
    sol = LS.Solve(np.array([[1],[2]]))
    if sol[0] != 2. : raise Exception()
    if sol[1] != 2. : raise Exception()

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
