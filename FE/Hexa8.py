# -*- coding: utf-8 -*-
__author__ = "Falipe Bordeu, Fabien Casenave"

import numpy as np

import BasicTools.FE.MaterialHelp as MH
from BasicTools.FE.FElement import FElement
from BasicTools.FE.FemHelp import Integral
import BasicTools.FE.ElementNames as ElementsNames

class Hexa8(FElement):
    def __init__(self):
        super(Hexa8,self).__init__()
        self.nnodes = 8
        self.dimensionality = 3
        self.name = ElementsNames.Hexaedron_8
        self.NumIntegRule = {}

    def GetShapeFunc(self, qcoor):
        xi = float(qcoor[0]);   eta = float(qcoor[1]); phi  = float(qcoor[2]);

        Nf = (1./8.)*np.array([(1-xi)*(1-eta)*(1-phi), (1+xi)*(1-eta)*(1-phi), (1+xi)*(1+eta)*(1-phi), (1-xi)*(1+eta)*(1-phi), (1-xi)*(1-eta)*(1+phi), (1+xi)*(1-eta)*(1+phi), (1+xi)*(1+eta)*(1+phi), (1-xi)*(1+eta)*(1+phi)])
        return Nf

    def ShapeFuncDerBaryCoord(self, qcoor):
       xi = float(qcoor[0]);   eta = float(qcoor[1]); phi  = float(qcoor[2]);
       Nfxi  = (1./8.)*np.array([(-1)  *(1-eta)*(1-phi),    (1)*(1-eta)*(1-phi),    (1)*(1+eta)*(1-phi),  (-1) *(1+eta)*(1-phi), (-1)  *(1-eta)*(1+phi),   (1) *(1-eta)*(1+phi), (1)   *(1+eta)*(1+phi), (-1)  *(1+eta)*(1+phi)])
       Nfeta = (1./8.)*np.array([(1-xi)*   (-1)*(1-phi), (1+xi)*   (-1)*(1-phi), (1+xi)*(1)    *(1-phi), (1-xi)*(1)    *(1-phi), (1-xi)*(-1)   *(1+phi), (1+xi)*(-1)   *(1+phi), (1+xi)*(1)    *(1+phi), (1-xi)*(1)    *(1+phi)])
       Nfphi = (1./8.)*np.array([(1-xi)*(1-eta)*(-1)   , (1+xi)*(1-eta)*(-1)   , (1+xi)*(1+eta)*(-1)   , (1-xi)*(1+eta)*(-1)   , (1-xi)*(1-eta)*(1)    , (1+xi)*(1-eta)*(1)    , (1+xi)*(1+eta)*(1)    , (1-xi)*(1+eta)*(1)])
       return np.array([Nfxi, Nfeta, Nfphi])

#---------------------------------------------------------------

    def GetBxByBz(self, Nfder, xcoor):
       xcoorx = xcoor[:,0]; xcoory = xcoor[:,1]; xcoorz = xcoor[:,2]

       J11 = np.vdot(Nfder[0],xcoorx); J21 = np.vdot(Nfder[1],xcoorx); J31 = np.vdot(Nfder[2],xcoorx)
       J12 = np.vdot(Nfder[0],xcoory); J22 = np.vdot(Nfder[1],xcoory); J32 = np.vdot(Nfder[2],xcoory)
       J13 = np.vdot(Nfder[0],xcoorz); J23 = np.vdot(Nfder[1],xcoorz); J33 = np.vdot(Nfder[2],xcoorz)

       Jdet = J11*J22*J33+J21*J32*J13+J31*J12*J23-J31*J22*J13-J11*J32*J23-J21*J12*J33

       Jinv = np.array([[J22*J33-J32*J23,J32*J13-J12*J33,J12*J23-J22*J13],
                        [J31*J23-J21*J33,J11*J33-J31*J13,J21*J13-J11*J23],
                        [J21*J32-J31*J22,J31*J12-J11*J32,J11*J22-J21*J12]])

       B = np.dot(Jinv,Nfder)/Jdet
       return B


#---------------------------------------------------------------

    def SetGaussQuad222IntRule(self):
       a = np.sqrt(1./3.)
       self.NumIntegRule['p'] = np.array([[-a,-a,-a],
                                          [ a,-a,-a],
                                          [-a, a,-a],
                                          [ a, a,-a],
                                          [-a,-a, a],
                                          [ a,-a, a],
                                          [-a, a, a],
                                          [ a, a, a]]);
       self.NumIntegRule['w'] = np.ones(8)/8.;
       self.NumIntegRule['nbGaussPoints'] = 8;

    def ComputeNfder(self):
       Nfer = []
       for i in range(self.NumIntegRule['nbGaussPoints']):
         Nfer.append(self.ShapeFuncDerBaryCoord(self.NumIntegRule['p'][i,:]))
       return Nfer



def CheckIntegrity():
    H = Hexa8()
    H.SetGaussQuad222IntRule()
    print(H.GetShapeFunc([0,0,0]))
    print(H.GetShapeFunc([-1,-1,-1]))
    print(H.ShapeFuncDerBaryCoord([0,0,0])[0])
    print(H.ComputeNfder()[0][0])
    return 'OK'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
