# -*- coding: utf-8 -*-

import numpy as np
import BasicTools.FE.MaterialHelp as MH
from BasicTools.FE.FElement import FElement
from BasicTools.FE.FemHelp import Integral
import BasicTools.FE.ElementNames as ElementsNames


class Hexa8Cuboid(FElement):
    def __init__(self):
        super(Hexa8Cuboid,self).__init__()
        self.delta = np.array([1,1,1]);
        self.nnodes = 8
        self.dimensionality = 3
        self.name = ElementsNames.Hexaedron_8+"Cuboid"

    def GetDetJack(self,qcoor):
        dx = float(self.delta[0]);   dy = float(self.delta[1]);  dz = float(self.delta[2]);
        return dx*dy*dz/8.;

    def GetShapeFunc(self, qcoor):
        xi = float(qcoor[0]);   eta = float(qcoor[1]); phi  = float(qcoor[2]);

        Nf = (1./8.)*np.array([(1-xi)*(1-eta)*(1-phi), (1+xi)*(1-eta)*(1-phi), (1+xi)*(1+eta)*(1-phi), (1-xi)*(1+eta)*(1-phi), (1-xi)*(1-eta)*(1+phi), (1+xi)*(1-eta)*(1+phi), (1+xi)*(1+eta)*(1+phi), (1-xi)*(1+eta)*(1+phi)])
        return Nf

    def ShapeFuncDer(self, qcoor):
       dx = float(self.delta[0]);   dy = float(self.delta[1]);  dz = float(self.delta[2]);
       dxidx = 2./dx;
       detady = 2./dy;
       dphidz = 2./dz;
       xi = float(qcoor[0]);   eta = float(qcoor[1]); phi  = float(qcoor[2]);

       Nfx = (1./8.)*dxidx *np.array([(-1)  *(1-eta)*(1-phi),    (1)*(1-eta)*(1-phi),    (1)*(1+eta)*(1-phi),  (-1) *(1+eta)*(1-phi), (-1)  *(1-eta)*(1+phi),   (1) *(1-eta)*(1+phi), (1)   *(1+eta)*(1+phi), (-1)  *(1+eta)*(1+phi)])
       Nfy = (1./8.)*detady*np.array([(1-xi)*   (-1)*(1-phi), (1+xi)*   (-1)*(1-phi), (1+xi)*(1)    *(1-phi), (1-xi)*(1)    *(1-phi), (1-xi)*(-1)   *(1+phi), (1+xi)*(-1)   *(1+phi), (1+xi)*(1)    *(1+phi), (1-xi)*(1)    *(1+phi)])
       Nfz = (1./8.)*dphidz*np.array([(1-xi)*(1-eta)*(-1)   , (1+xi)*(1-eta)*(-1)   , (1+xi)*(1+eta)*(-1)   , (1-xi)*(1+eta)*(-1)   , (1-xi)*(1-eta)*(1)    , (1+xi)*(1-eta)*(1)    , (1+xi)*(1+eta)*(1)    , (1-xi)*(1+eta)*(1)])
       return np.array([Nfx, Nfy, Nfz])

#---------------- Derivative matrices mass matrices ----------

    def IsoDispB(self,qcoor,pos):
        [Nfx, Nfy, Nfz] = self.ShapeFuncDer(qcoor)
        B = np.array([[Nfx[0], 0     , 0     , Nfx[1], 0     , 0     , Nfx[2], 0     , 0     , Nfx[3], 0     , 0     , Nfx[4], 0     , 0     , Nfx[5], 0     , 0     , Nfx[6], 0     , 0     , Nfx[7], 0     , 0     ],
                      [0     , Nfy[0], 0     , 0     , Nfy[1], 0     , 0     , Nfy[2], 0     , 0     , Nfy[3], 0     , 0     , Nfy[4], 0     , 0     , Nfy[5], 0     , 0     , Nfy[6], 0     , 0     , Nfy[7], 0     ],
                      [0     , 0     , Nfz[0], 0     , 0     , Nfz[1], 0     , 0     , Nfz[2], 0     , 0     , Nfz[3], 0     , 0     , Nfz[4], 0     , 0     , Nfz[5], 0     , 0     , Nfz[6], 0     , 0     , Nfz[7]],
                      [Nfy[0], Nfx[0], 0     , Nfy[1], Nfx[1], 0     , Nfy[2], Nfx[2], 0     , Nfy[3], Nfx[3], 0     , Nfy[4], Nfx[4], 0     , Nfy[5], Nfx[5], 0     , Nfy[6], Nfx[6], 0     , Nfy[7], Nfx[7], 0     ],
                      [0     , Nfz[0], Nfy[0], 0     , Nfz[1], Nfy[1], 0     , Nfz[2], Nfy[2], 0     , Nfz[3], Nfy[3], 0     , Nfz[4], Nfy[4], 0     , Nfz[5], Nfy[5], 0     , Nfz[6], Nfy[6], 0     , Nfz[7], Nfy[7]],
                      [Nfz[0], 0     , Nfx[0], Nfz[1], 0     , Nfx[1], Nfz[2], 0     , Nfx[2], Nfz[3], 0     , Nfx[3], Nfz[4], 0     , Nfx[4], Nfz[5], 0     , Nfx[5], Nfz[6], 0     , Nfx[6], Nfz[7], 0     , Nfx[7]]]);
        return B, self.GetDetJack(qcoor)

    def IsoLaplaceB(self,qcoor,pos):
         [Nfx, Nfy, Nfz] = self.ShapeFuncDer(qcoor)
         B = np.array([[Nfx[0], Nfx[1], Nfx[2], Nfx[3], Nfx[4], Nfx[5], Nfx[6], Nfx[7]],
                   [Nfy[0], Nfy[1], Nfy[2], Nfy[3], Nfy[4], Nfy[5], Nfy[6], Nfy[7]],
                   [Nfz[0], Nfz[1], Nfz[2], Nfz[3], Nfz[4], Nfz[5], Nfz[6], Nfz[7]]]);
         return B, self.GetDetJack(qcoor)


#----------------mass matrices ----------
    def IsoDispM(self,qcoor,pos):
        N = self.GetShapeFunc(qcoor)
        B = np.array([[N[0], 0   , 0   , N[1], 0   , 0   , N[2], 0   , 0   , N[3], 0   , 0   , N[4], 0   , 0   , N[5], 0   , 0   , N[6], 0   , 0   , N[7], 0   , 0     ],
                   [0   , N[0], 0   , 0   , N[1], 0   , 0   , N[2], 0   , 0   , N[3], 0   , 0   , N[4], 0   , 0   , N[5], 0   , 0   , N[6], 0   , 0   , N[7], 0     ],
                   [0   , 0   , N[0], 0   , 0   , N[1], 0   , 0   , N[2], 0   , 0   , N[3], 0   , 0   , N[4], 0   , 0   , N[5], 0   , 0   , N[6], 0   , 0   , N[7]]])
        return B, self.GetDetJack(qcoor)

    def IsoLaplaceM(self,qcoor,pos):
         N = self.GetShapeFunc(qcoor)
         B = np.array([[N[0], N[1], N[2], N[3], N[4], N[5], N[6], N[7]]])
         return B,  self.GetDetJack(qcoor)

#-------------------------

    def GetIsotropLaplaceK(self,k):
        #IsoHexaCubeKLaplace(k,delta):
        K = MH.LaplaceOrtho(k,k,k)
        return Integral(K,self.IsoLaplaceB,self,self.nnodes)

    def GetOrthoLaplaceK(self,k):
        K = MH.LaplaceOrtho(k[0],k[1],k[2])
        return Integral(K,self.IsoLaplaceB,self,self.nnodes)


    def GetIsotropLaplaceM(self,rho):
        K = np.identity(1)*rho
        return Integral(K,self.IsoLaplaceM,self,self.nnodes)

#-------------------------

    def GetIsotropDispK(self,E,nu):
        #IsoHexaCubeK
        k = MH.HookeIso(E,nu)
        return Integral(k,self.IsoDispB,self,self.nnodes*3)

    def GetIsotropDispM(self,rho):
         K = np.identity(3)*rho
         return Integral(K,self.IsoDispM,self,self.nnodes*3)

#-------------------------

def CheckIntegrity():

    import numpy.linalg as lin
    myElement = Hexa8Cuboid()
    myElement.delta = [1,1,1]
    KK = myElement.GetIsotropDispK(1,0.3)
    VV =lin.eig(KK)
    if np.sum(VV[0].real < 1e-10) != 6:
        return "not ok "# pragma: no cover

    K2 = myElement.GetIsotropLaplaceK(1)
    VV =lin.eig(K2)
    if np.sum(VV[0].real < 1e-10) != 1:
        return "not ok "# pragma: no cover

    #print(VV[0])
    c = myElement.GetOrthoLaplaceK([1, 2, 3])
    print(c)

    b,_ = myElement.IsoLaplaceB([0,0,0],None)
    #print(myElement.GetDetJack())
    #print(b)
    u = np.matrix([0, 0+2 ,0+2+3 ,0+3 ,1 ,1+2 ,1+2+3, 1+3 ]).T;
    res1 =(b*u).T
    res2 = [2./myElement.delta[0],3./myElement.delta[1],1./myElement.delta[2] ]
    if np.all(res1==res2):
        return 'OK'
    return 'notOK'# pragma: no cover

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
