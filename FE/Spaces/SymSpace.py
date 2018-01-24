# -*- coding: utf-8 -*-
from BasicTools.FE.FElement import FElement
from sympy import Symbol
import numpy as np


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


class SymSpaceBase(SpaceBase):
    def __init__(self):
        super(SymSpaceBase,self).__init__()
        #----symbolic part use Create() to pass to the next step ------------
        self.xi  = Symbol("xi")
        self.eta = Symbol("eta")
        self.phi = Symbol("phi")
        self.symN = None
        self.symdNdxi = None
        self.symdNdeta = None
        self.symdNdphi = None

        #---- Generic part use SetIntegrationRule(...) rule to generate shape funct values and derivatives ----------
        #self.GetShapeFunc = None (self,qcoor)
        #---- Geneal element part, use xcoor with the coordinate of the elemetn to get values

    def GetNumberOfShapeFunctions(self):
        return len(self.symN)

    def Create(self):

        self.lcoords = (self.xi,self.eta,self.phi)
        self.lcoords = self.lcoords[0:self.GetDimensionality()]

#        if self.geoSupport.dimensionality >= 1:
#            self.symdNdxi = self.symN.diff(self.xi)
#        if self.geoSupport.dimensionality >= 2:
#            self.symdNdeta = self.symN.diff(self.eta)
#        if self.geoSupport.dimensionality >= 3:
#            self.symdNdphi = self.symN.diff(self.phi)

        self.fct_N = [None]*self.GetNumberOfShapeFunctions()
        self.fct_dNdxi = [[None]*self.GetDimensionality() for i in range(self.GetNumberOfShapeFunctions() ) ]
#        print self.fct_dNdxi
#        if self.GetDimensionality() == 1:
#            self.coords = [ self.xi]
#        if self.GetDimensionality() == 1:
#            self.coords = [ self.xi, self.eta]
#        if self.GetDimensionality() == 1:
#            self.coords = [ self.xi, self.eta, self.phi]

        for i in range(self.GetNumberOfShapeFunctions()) :
            self.fct_N[i] = lambda p,i=i: self.symN[i].subs(zip(self.lcoords,p)).evalf()
            for j in range(self.GetDimensionality() ) :
                func = self.symN[i].diff(self.lcoords[j])
                self.fct_dNdxi[i][j] =  lambda p,func=func: func.subs(zip(self.lcoords,p)).evalf()
        #print(self.fct_dNdxi)

    def SetIntegrationRule(self, points, weights):
       self.int_Weights = weights
       self.int_Points = points
       self.int_nbPoints = len(weights)

       nsf = self.GetNumberOfShapeFunctions()
       dim = self.GetDimensionality()

       self.valN = [None]*self.int_nbPoints
       self.valdphidxi = [None]*self.int_nbPoints

       for pp in range(self.int_nbPoints):
           point = points[pp]
           self.valN[pp] = self.GetShapeFunc(point)
           #self.valdphidxi[pp] = np.empty((dim,nsf),dtype=np.float)
           self.valdphidxi[pp] = self.GetShapeFuncDer(point)

    def GetShapeFunc(self,qcoor):
        res = np.array([self.fct_N[i](qcoor) for i in range(self.GetNumberOfShapeFunctions()) ], dtype=np.float)
        return res

    def GetShapeFuncDer(self,qcoor):
        nsf = self.GetNumberOfShapeFunctions()
        dim = self.GetDimensionality()
        res = np.empty((dim,nsf),dtype=np.float)
        for fct in range(nsf):
            for direc in range(dim):
                #print self.fct_dNdxi[fct]
                res[direc,fct]   = self.fct_dNdxi[fct][direc](qcoor)
        return res

    def GetNormal(self,Jack):
        # Edge in 2D
        if Jack.shape[0] == 1 and Jack.shape[1] == 2 :
            res = np.array([Jack[1,:] -Jack[0,:]],dtype =np.float)
        # surface in 3D
        elif Jack.shape[0] == 2 and Jack.shape[1] == 3 :
            res =  np.cross (Jack[0,:],Jack[1,:])

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
       elif dim == 0:
           Jdet = 1
       elif dim == 1:
           Jdet = np.linalg.norm(Jack)
       elif dim == 2:
           Jdet = np.linalg.norm(np.cross (Jack[0,:],Jack[1,:]))

       q,r = np.linalg.qr(Jack)
       qt = q.T

       jinv = lambda vec,qt=qt,r=r: np.linalg.lstsq(r, np.dot(qt,vec))[0]

       return Jack, Jdet, jinv

#    def GetBxByBzI(self,I,xcoor):
#        Nfder = [self.valdNdxi[I], self.valdNdeta[I], self.valdNdphi[I]]
#        return self.GetBxByBz(Nfder,xcoor)


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
