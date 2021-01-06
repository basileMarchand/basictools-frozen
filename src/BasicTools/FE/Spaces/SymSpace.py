# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from sympy import Symbol, DiracDelta, Matrix
from sympy.utilities.lambdify import lambdify

from BasicTools.FE.Spaces.SpaceBase import SpaceBase
from BasicTools.FE.Spaces.SpaceBase import SpaceAtIntegrationPoints

class SymSpaceBase(SpaceBase):
    def __init__(self):
        super(SymSpaceBase,self).__init__()
        #----symbolic part use Create() to pass to the next step ------------
        self.xi  = Symbol("xi")
        self.eta = Symbol("eta")
        self.phi = Symbol("phi")
        self.symN = None
        self.symdNdxi = None
        self.__created__ = False

    def __getstate__(self):
        return None

    def __setstat__(self,state):
        self.__init__()


    def GetNumberOfShapeFunctions(self):
        return len(self.symN)

    def Create(self,force=False):
        if self.__created__ and not force:
            return
        self.__created__ = True
        def DiractDeltaNumeric(data,der=None):
            if abs(data)>1e-15 :
                return 0
            else:
                return 1

        allcoords = (self.xi,self.eta,self.phi)
        self.lcoords = tuple(  (self.xi,self.eta,self.phi)[x] for x in range(self.GetDimensionality())  )
        nbSF = self.GetNumberOfShapeFunctions()
        nbDim = self.GetDimensionality()


        subsList = [ (DiracDelta(0),1.), (DiracDelta(0,1),1.), (DiracDelta(0,2),1.) ]
        lambdifyList =  [ {"DiracDelta":DiractDeltaNumeric}, "numpy"]

        ############# shape function treatement ########################

        clean_N = self.symN.subs(subsList)
        self.fct_N_Matrix =  lambdify(allcoords,[ clean_N[i] for i in  range(nbSF)  ], lambdifyList )

        ############# shape functions first derivative #################
        self.symdNdxi = [[None]*nbSF for i in range(nbDim)]

        for i in range(nbDim ) :
            for j in range(nbSF) :
                self.symdNdxi[i][j] = self.symN[j].diff(self.lcoords[i])

        self.symdNdxi = Matrix(self.symdNdxi)
        if self.symdNdxi.shape == (0,0):
           self.fct_dNdxi_Matrix = lambda xi,chi,phi: np.empty((0,0))
        else:
           self.fct_dNdxi_Matrix =  lambdify(allcoords,self.symdNdxi.subs(subsList) , lambdifyList )
        ############ shape functions second derivative ################

        self.symdNdxidxi = [ None ]*nbSF
        self.fct_dNdxidxi_Matrix = [ None ]*nbSF

        for i in range(nbSF) :
            self.symdNdxidxi[i] = [[0]*nbDim for j in range(nbDim ) ]
            for j in range(nbDim ) :
                for k in range(self.GetDimensionality() ) :
                    func = self.symN[i].diff(self.lcoords[j]).diff(self.lcoords[k])
                    self.symdNdxidxi[i][j][k] = func
            self.symdNdxidxi[i] = Matrix(self.symdNdxidxi[i])

            self.fct_dNdxidxi_Matrix[i] = lambdify(allcoords,self.symdNdxidxi[i].subs(subsList) , lambdifyList )

    def SetIntegrationRule(self, points, weights):
        self.Create()

        res = SpaceAtIntegrationPoints()
        res.SetIntegrationRule(self,points,weights)
        return res

    def GetPosOfShapeFunction(self,i,Xi):
        valN = self.GetShapeFunc(self.posN[i,:])
        return np.dot(valN,Xi).T

    def GetShapeFunc_default(self,xi=0,chi=0,phi=0):
        return self.fct_N_Matrix(xi,chi,phi)

    def GetShapeFunc(self,qcoor):
        return np.asarray(self.GetShapeFunc_default(*qcoor), dtype=float)

    def GetShapeFuncDer_default(self,xi=0,chi=0,phi=0):
        return self.fct_dNdxi_Matrix(xi,chi,phi)

    def GetShapeFuncDer(self,qcoor):
        return np.asarray(self.GetShapeFuncDer_default(*qcoor), dtype=float)

    def GetShapeFuncDerDer_default(self,xi=0,chi=0,phi=0):
        return [ np.asarray(x(xi,chi,phi),dtype=float) for x in self.fct_dNdxidxi_Matrix ]

    def GetShapeFuncDerDer(self,qcoor):
        return self.GetShapeFuncDerDer_default(*qcoor)

def CheckIntegrity():
    return "ok"
