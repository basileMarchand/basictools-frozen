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

       #tetra:
       """self.valdphidxi = np.array([[[-1,-0.3333333333,0,1.333333333,0.6666666667,-0.6666666667,-0.6666666667,0.6666666667,0,0],[-1,0,-0.3333333333,-0.6666666667,0.6666666667,1.333333333,-0.6666666667,0,0.6666666667,0],[-1,0,0,-0.6666666667,0,-0.6666666667,1.333333333,0.6666666667,0.6666666667,-0.3333333333]],
[[0.3333333333,1,0,-1.333333333,0.6666666667,-0.6666666667,-0.6666666667,0.6666666667,0,0],[0.3333333333,0,-0.3333333333,-2,2,2.220446049e-16,-0.6666666667,0,0.6666666667,0],[0.3333333333,0,0,-2,0,-0.6666666667,2.220446049e-16,2,0.6666666667,-0.3333333333]],
[[0.3333333333,-0.3333333333,0,2.220446049e-16,2,-2,-0.6666666667,0.6666666667,0,0],[0.3333333333,0,1,-0.6666666667,0.6666666667,-1.333333333,-0.6666666667,0,0.6666666667,0],[0.3333333333,0,0,-0.6666666667,0,-2,2.220446049e-16,0.6666666667,2,-0.3333333333]],
[[0.3333333333,-0.3333333333,0,3.330669074e-16,0.6666666667,-0.6666666667,-2,2,0,0],[0.3333333333,0,-0.3333333333,-0.6666666667,0.6666666667,3.330669074e-16,-2,0,2,0],[0.3333333333,0,0,-0.6666666667,0,-0.6666666667,-1.333333333,0.6666666667,0.6666666667,1]],
[[0,0,0,0,1,-1,-1,1,0,0],[0,0,0,-1,1,0,-1,0,1,0],[0,0,0,-1,0,-1,0,1,1,0]]])"""
       #tri:
       """self.valdphidxi = np.array([[[0.5675879273,-1.351381891,0.7837939637,1.783793964,0,-1.783793964],[0.5675879273,-1.783793964,0,1.783793964,0.7837939637,-1.351381891]],
[[-0.7837939637,1.351381891,-0.5675879273,1.783793964,0,-1.783793964],[-0.7837939637,-0.4324120727,0,0.4324120727,0.7837939637,0]],
[[-0.7837939637,0,0.7837939637,0.4324120727,0,-0.4324120727],[-0.7837939637,-1.783793964,0,1.783793964,-0.5675879273,1.351381891]],
[[-2.267390292,2.901085438,-0.633695146,0.366304854,0,-0.366304854],[-2.267390292,-0.366304854,0,0.366304854,-0.633695146,2.901085438]],
[[0.633695146,-2.901085438,2.267390292,0.366304854,0,-0.366304854],[0.633695146,-3.267390292,0,3.267390292,-0.633695146,1.110223025e-16]],
[[0.633695146,2.775557562e-16,-0.633695146,3.267390292,0,-3.267390292],[0.633695146,-0.366304854,0,0.366304854,2.267390292,-2.901085438]]])"""


       """print("self.valdphidxi2 =", self.valdphidxi2)
       print("self.valdphidxi =", self.valdphidxi)
       print("np.array(self.valdphidxi).shape, np.linalg.norm(np.array(self.valdphidxi)) =", np.array(self.valdphidxi).shape, np.linalg.norm(np.array(self.valdphidxi)))
       print("np.array(self.valdphidxi2).shape, np.linalg.norm(np.array(self.valdphidxi2)) =", np.array(self.valdphidxi2).shape, np.linalg.norm(np.array(self.valdphidxi2)))
       print("====")
       print(np.array(self.valN2).shape, np.linalg.norm(np.array(self.valN2)))"""

       #tetra
       """self.valN = np.array([[1.110223025e-16,-0.1111111111,-0.1111111111,0.3333333333,0.1111111111,0.3333333333,0.3333333333,0.1111111111,0.1111111111,-0.1111111111],
[-0.1111111111,-0,-0.1111111111,0.3333333333,0.3333333333,0.1111111111,0.1111111111,0.3333333333,0.1111111111,-0.1111111111],
[-0.1111111111,-0.1111111111,-0,0.1111111111,0.3333333333,0.3333333333,0.1111111111,0.1111111111,0.3333333333,-0.1111111111],
[-0.1111111111,-0.1111111111,-0.1111111111,0.1111111111,0.1111111111,0.1111111111,0.3333333333,0.3333333333,0.3333333333,-0],
[-0.125,-0.125,-0.125,0.25,0.25,0.25,0.25,0.25,0.25,-0.125]])"""
       #tri:
       """self.valN = np.array([[-0.08473049309,0.1928335113,-0.04820837782,0.7954802262,-0.04820837782,0.1928335113],
[-0.04820837782,0.1928335113,-0.08473049309,0.1928335113,-0.04820837782,0.7954802262],
[-0.04820837782,0.7954802262,-0.04820837782,0.1928335113,-0.08473049309,0.1928335113],
[0.517632342,0.299215231,-0.07480380775,0.03354481152,-0.07480380775,0.299215231],
[-0.07480380775,0.299215231,0.517632342,0.299215231,-0.07480380775,0.03354481152],
[-0.07480380775,0.03354481152,-0.07480380775,0.299215231,0.517632342,0.299215231]])"""


       """print("self.valN2 =", self.valN2)
       print("self.valN =", self.valN)
       print("np.array(self.valN).shape, np.linalg.norm(np.array(self.valN)) =", np.array(self.valN).shape, np.linalg.norm(np.array(self.valN)))
       print("np.array(self.valN2).shape, np.linalg.norm(np.array(self.valN2)) =", np.array(self.valN2).shape, np.linalg.norm(np.array(self.valN2)))
       print("====")
       print("====")"""
       

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
