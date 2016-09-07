# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:58:44 2016

@author: d584808
"""

# -*- coding: utf-8 -*-

import numpy as np
from OTTools.FE.FElement import FElement
from OTTools.FE.FemHelp import Integral
import OTTools.FE.MaterialHelp as MH
import numpy.linalg as linalg
import OTTools.FE.ElementNames as ElementsNames

class Tri3(FElement):
    def __init__(self):
        super(Tri3,self).__init__()
        self.delta = np.array([1,1]);
        self.nnodes = 3
        self.dimensionality = 2
        self.name = ElementsNames.Triangle_3

    def GetDetJack(self,qcoor, pos=None):
        if pos == None:
            raise(Exception('Need position to calculate gradient'))

        dx = float(self.delta[0]);   dy = float(self.delta[1]);
        return dx*dy/4.;

    def GetShapeFunc(self, qcoor):
        xi = float(qcoor[0]);   eta = float(qcoor[1]);

        Nf = np.array([(1-xi-eta), xi,eta])
        return Nf

    #def compute_gradient(pg,e,dim):
        #[B det_jac]=
    def compute_gradient(self, qcoor, pos=None):
        #[B det_jac]=
       if pos is  None:
            raise(Exception('Need position to calculate gradient'))

       dN_xi = np.array([-1, 1, 0])
       dN_eta = np.array([-1, 0, 1])

       jac = [[ dN_xi.dot(pos[:,0]) , dN_eta.dot(pos[:,0])],
              [ dN_xi.dot(pos[:,1]) , dN_eta.dot(pos[:,1])]   ];

       #print(pos)
       #print(jac)
       jaci = linalg.inv(jac);
       dN = jaci.T.dot(np.vstack((dN_xi,dN_eta)));
       dN_x = dN[0,:];
       dN_y = dN[1,:];
       return (dN_x,dN_y, linalg.det(jac) )

    def IsoDispB(self,qcoor,pos=None):
        [dN_x, dN_y, det_jac] = self.compute_gradient(qcoor,pos)
        B = np.array([[dN_x[0], 0      , dN_x[1], 0      , dN_x[2], 0      ],
                      [0      , dN_y[0], 0      , dN_y[1], 0      , dN_y[2]],
                      [dN_y[0], dN_x[0], dN_y[1], dN_x[1], dN_y[2], dN_x[2]]]);
        return B,det_jac

    def IsoLaplaceB(self,qcoor,pos=None):
         [dN_x, dN_y, det_jac] = self.compute_gradient(qcoor,pos)
         B = np.array([[dN_x[0], dN_x[1], dN_x[2]],
                       [dN_y[0], dN_y[1], dN_y[2]]]);
         return B,det_jac
#-------------------------

    def GetIsotropLaplaceK(self,k,pos):
        #IsoHexaCubeKLaplace(k,delta):
        K = MH.LaplaceOrtho(k,k,dim =2)
        return Integral(K,self.IsoLaplaceB,self,self.nnodes,pos)

    def GetOrthoLaplaceK(self,k,pos):
        K = MH.LaplaceOrtho(k[0],k[1],dim =2)
        return Integral(K,self.IsoLaplaceB,self,self.nnodes,pos)


    def GetIsotropLaplaceM(self,rho,pos):
        K = np.identity(1)*rho
        return Integral(K,self.IsoLaplaceM,self,self.nnodes,pos)
    #-------------------------

    def GetIsotropDispK(self,E,nu,pos):
        #IsoHexaCubeK
        k = MH.HookeIso(E,nu, dim=2)
        return Integral(k,self.IsoDispB,self,self.nnodes*2,pos)

    def GetIsotropDispM(self,rho,pos):
         K = np.identity(2)*rho
         return Integral(K,self.IsoDispM,self,self.nnodes*2,pos)

def CheckIntegrity():

    import numpy.linalg as lin
    myElement = Tri3()
    pos = np.array([ [0,0], [1,0] , [0,1] ])
    KK = myElement.GetIsotropDispK(1,0.3,pos)
    print("--------------------------")
    print(KK)
    VV =lin.eig(KK)
    #print(VV)
    if np.sum(VV[0].real < 1e-10) != 3:
        return "not ok"# pragma: no cover

    K2 = myElement.GetIsotropLaplaceK(1,pos)
    VV =lin.eig(K2)
    if np.sum(VV[0].real < 1e-10) != 1:
        VV[0].real
        return "not ok "# pragma: no cover

    #print(VV[0])
    c = myElement.GetOrthoLaplaceK([1, 2],pos)
    print(c)

    b,_ = myElement.IsoLaplaceB([0,0],pos)
    #print(myElement.GetDetJack())
    print(b)
    u = np.matrix([0,1,2],dtype=float).T;
    res1 =(b*u).T
    res2 = [1./myElement.delta[0],2./myElement.delta[1] ]
    if np.all(res1==res2):
        return 'OK'
    return 'not OK  '# pragma: no cover

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover