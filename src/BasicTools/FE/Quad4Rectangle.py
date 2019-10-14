# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       
import numpy as np

import BasicTools.Containers.ElementNames as ElementsNames

from BasicTools.FE.FElement import FElement
from BasicTools.FE.FemHelp import Integral
import BasicTools.FE.MaterialHelp as MH


class Quad4Rectangle(FElement):
    def __init__(self):
        super(Quad4Rectangle,self).__init__()
        self.delta = np.array([1,1]);
        self.nnodes = 4
        self.dimensionality = 2
        self.name = ElementsNames.Quadrangle_4+"Rectangle"

    def GetDetJack(self,qcoor):
        dx = float(self.delta[0]);   dy = float(self.delta[1]);
        return dx*dy/4.;

    def GetShapeFunc(self, qcoor):
        xi = float(qcoor[0]);   eta = float(qcoor[1]);

        Nf = (1./4.)*np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
        return Nf

    def ShapeFuncDer(self, qcoor):
       dx = float(self.delta[0]);   dy = float(self.delta[1]);
       dxidx = 2./dx;
       detady = 2./dy;
       xi = float(qcoor[0]);   eta = float(qcoor[1]);

       Nfx = (1./4.)*dxidx *np.array([(-1)  *(1-eta),    (1)*(1-eta),    (1)*(1+eta),  (-1) *(1+eta)])
       Nfy = (1./4.)*detady*np.array([(1-xi)*   (-1), (1+xi)*   (-1), (1+xi)*(1)    , (1-xi)*(1)    ])
       return np.array([Nfx, Nfy])

#---------------- Derivative matrices ----------

    def IsoDispB(self,qcoor,pos=None):
        [Nfx, Nfy] = self.ShapeFuncDer(qcoor)
        B = np.array([[Nfx[0], 0     , Nfx[1], 0     , Nfx[2], 0     , Nfx[3], 0     ],
                      [0     , Nfy[0], 0     , Nfy[1], 0     , Nfy[2], 0     , Nfy[3]],
                      [Nfy[0], Nfx[0], Nfy[1], Nfx[1], Nfy[2], Nfx[2], Nfy[3], Nfx[3]]]);
        return B, self.GetDetJack(qcoor)

    def IsoLaplaceB(self,qcoor,pos=None):
         [Nfx, Nfy] = self.ShapeFuncDer(qcoor)
         B = np.array([[Nfx[0], Nfx[1], Nfx[2], Nfx[3]],
                       [Nfy[0], Nfy[1], Nfy[2], Nfy[3]]]);
         return B, self.GetDetJack(qcoor)


#----------------mass matrices ----------
    def IsoDispM(self,qcoor,pos=None):
        N = self.GetShapeFunc(qcoor)
        B = np.array([[N[0], 0   , N[1], 0   , N[2], 0   , N[3], 0   ],
                      [0   , N[0], 0   , N[1], 0   , N[2], 0   , N[3]]])
        return B, self.GetDetJack(qcoor)

    def IsoLaplaceM(self,qcoor,pos=None):
         N = self.GetShapeFunc(qcoor)
         B = np.array([[N[0], N[1], N[2], N[3]]])
         return B, self.GetDetJack(qcoor)

#-------------------------

    def GetIsotropLaplaceK(self,k):
        #IsoHexaCubeKLaplace(k,delta):
        K = MH.LaplaceOrtho(k,k,dim =2)
        return Integral(K,self.IsoLaplaceB,self,self.nnodes)

    def GetOrthoLaplaceK(self,k):
        K = MH.LaplaceOrtho(k[0],k[1],dim =2)
        return Integral(K,self.IsoLaplaceB,self,self.nnodes)


    def GetIsotropLaplaceM(self,rho):
        K = np.identity(1)*rho
        return Integral(K,self.IsoLaplaceM,self,self.nnodes)

#-------------------------

    def GetIsotropDispK(self,E,nu):
        #IsoHexaCubeK
        k = MH.HookeIso(E,nu, dim=2)
        return Integral(k,self.IsoDispB,self,self.nnodes*2)

    def GetIsotropDispM(self,rho):
         K = np.identity(2)*rho
         return Integral(K,self.IsoDispM,self,self.nnodes*2)

#-------------------------

def CheckIntegrity():

    import numpy.linalg as lin
    myElement = Quad4Rectangle()
    myElement.delta = [1,1]
    KK = myElement.GetIsotropDispK(1,0.3)
    VV =lin.eig(KK)
    if np.sum(VV[0].real < 1e-10) != 3:
        return "not ok"# pragma: no cover

    K2 = myElement.GetIsotropLaplaceK(1)
    VV =lin.eig(K2)
    if np.sum(VV[0].real < 1e-10) != 1:
        return "not ok "# pragma: no cover

    #print(VV[0])
    c = myElement.GetOrthoLaplaceK([1, 2])
    print(c)

    b,_ = myElement.IsoLaplaceB([0,0])
    #print(myElement.GetDetJack())
    #print(b)
    u = np.matrix([0, 1 ,1+2 ,0+2  ]).T;
    res1 =(b*u).T
    res2 = [1./myElement.delta[0],2./myElement.delta[1] ]
    if np.all(res1==res2):
        return 'OK'
    return 'not OK  '# pragma: no cover

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover