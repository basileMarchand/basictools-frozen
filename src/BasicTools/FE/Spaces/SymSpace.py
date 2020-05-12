# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from sympy import Symbol, DiracDelta, Matrix
from sympy.utilities.lambdify import lambdify

from BasicTools.FE.FElement import FElement

class SpaceClasification(object):
    def __init__(self,family = "Lagrange", degree = 1, discontinuous = False):
        self.family = family
        self.degree = degree
        self.discontinuous = discontinuous
    def GetName(self):
        return self.family[0] + str(self.degree)+ ("D" if self.discontinuous else "")

class SpaceBase(FElement):
    def __init__(self):
        self.classification = SpaceClasification()
    def GetName(self) :
        return self.geoSupport.name +"_" + self.classification.GetName()
    def GetDimensionality(self):
        return self.geoSupport.dimensionality

    def ClampParamCoorninates(self,xietaphi):
        res = xietaphi.copy()
        for cpt, d in enumerate(xietaphi):
            if cpt < self. GetDimensionality():
                res[cpt] = max(0.,d)
                res[cpt] = min(1.,res[cpt])
            else:
                res[cpt] = 0
        return res

class SymSpaceBase(SpaceBase):
    def __init__(self):
        super(SymSpaceBase,self).__init__()
        #----symbolic part use Create() to pass to the next step ------------
        self.xi  = Symbol("xi")
        self.eta = Symbol("eta")
        self.phi = Symbol("phi")
        self.symN = None
        self.symdNdxi = None

    def GetNumberOfShapeFunctions(self):
        return len(self.symN)

    def Create(self):

        def DiractDeltaNumeric(data,der=None):
            if data :
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
       self.int_Weights = weights
       self.int_Points = points
       self.int_nbPoints = len(weights)

       self.valN = [None]*self.int_nbPoints
       self.valdphidxi = [None]*self.int_nbPoints

       for pp in range(self.int_nbPoints):
           point = points[pp]
           self.valN[pp] = self.GetShapeFunc(point)
           self.valdphidxi[pp] = self.GetShapeFuncDer(point)

    def GetPosOfShapeFunction(self,i,Xi):
        valN = self.GetShapeFunc(self.posN[i,:])
        return np.dot(valN,Xi).T

    def GetShapeFunc_default(self,xi=0,chi=0,phi=0):
        return self.fct_N_Matrix(xi,chi,phi)

    def GetShapeFunc(self,qcoor):
        return np.array(self.GetShapeFunc_default(*qcoor), dtype=float)

    def GetShapeFuncDer_default(self,xi=0,chi=0,phi=0):
        return self.fct_dNdxi_Matrix(xi,chi,phi)

    def GetShapeFuncDer(self,qcoor):
        return np.array(self.GetShapeFuncDer_default(*qcoor), dtype=float)


    def GetShapeFuncDerDer_default(self,xi=0,chi=0,phi=0):
        nsf = self.GetNumberOfShapeFunctions()
        dim = self.GetDimensionality()
        return [ np.array(x(xi,chi,phi),dtype=float) for x in self.fct_dNdxidxi_Matrix ]

    def GetShapeFuncDerDer(self,qcoor):
        return self.GetShapeFuncDerDer_default(*qcoor)

    def GetNormal(self,Jack):
        # Edge in 2D
        if Jack.shape[0] == 1 and Jack.shape[1] == 2 :
            res = np.array([Jack[0,1],-Jack[0,0]],dtype =np.float)
            #res = np.array([Jack[1,:] -Jack[0,:]],dtype =np.float) #ANCIENNE VERSION
        # surface in 3D
        elif Jack.shape[0] == 2 and Jack.shape[1] == 3 :
            res =  np.cross(Jack[0,:],Jack[1,:])
        else:
            print("Shape of Jacobian not coherent. Possible error: an elset has the same name of the considered faset")

        #normalisation
        res /= np.linalg.norm(res)
        return res

    def GetJackAndDetI(self, pp, xcoor):
       return self.GetJackAndDet(self.valdphidxi[pp], xcoor)

    def GetJackAndDet(self, Nfder, xcoor):

       Jack = np.dot(Nfder,xcoor)

       dim = self.GetDimensionality()

       s = xcoor.shape[1]

       if dim == s:
           Jdet = np.linalg.det(Jack)

           if dim == 3:
               def jinv(vec,jack=Jack):
                   m1, m2, m3, m4, m5, m6, m7, m8, m9 = jack.flatten()
                   determinant = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9
                   return np.dot(np.array([[m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5],
                       [m6*m7-m4*m9, m1*m9-m3*m7, m3*m4-m1*m6],
                       [m4*m8-m5*m7, m2*m7-m1*m8, m1*m5-m2*m4]]),vec) /determinant
               return Jack,Jdet,jinv

       elif dim == 0:
           Jdet = 1
       elif dim == 1:
           Jdet = np.linalg.norm(Jack)
       elif dim == 2:
           Jdet = np.linalg.norm(np.cross (Jack[0,:],Jack[1,:]))

       q,r = np.linalg.qr(Jack)
       qt = q.T

       jinv = lambda vec,qt=qt,r=r: np.linalg.lstsq(r, np.dot(qt,vec), rcond = None)[0]

       return Jack, Jdet, jinv

    def Eval_FieldI(self,I,Xi,J,Jinv,der=-1):

        if der ==-1:
            #print (self.valN[I])
            #print (Xi)
            res = np.dot(self.valN[I],Xi).T
        else:
            res = np.dot(Jinv(self.valdphidxi[I])[der,:],Xi)
        return res

    def ComputeNfder(self):
       Nfer = []
       for i in range(self.NumIntegRule['nbGaussPoints']):
         Nfer.append(self.GetShapeFuncDer(self.NumIntegRule['p'][i,:]))
       return Nfer

    def GetBMeca(self,BxByBzI):
        nbsf = self.GetNumberOfShapeFunctions()
        B = np.zeros((6,nbsf*3), dtype=np.float)
        for i in range(nbsf):
             B[0,i] =   BxByBzI[0,i]
             B[1,i+nbsf] = BxByBzI[1,i]
             B[2,i+2*nbsf] = BxByBzI[2,i]

             B[3,i] =   BxByBzI[1,i]
             B[3,i+nbsf] =   BxByBzI[0,i]

             B[4,i] =   BxByBzI[2,i]
             B[4,i+2*nbsf] =   BxByBzI[0,i]

             B[5,i+nbsf] =   BxByBzI[2,i]
             B[5,i+2*nbsf] =   BxByBzI[1,i]

def CheckIntegrity():
    return "ok"
