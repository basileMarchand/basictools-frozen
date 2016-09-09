# -*- coding: utf-8 -*-
import numpy as np

def HookeIso(E,nu, dim = 3, PlaneStress = True):
    if dim == 2:
        if PlaneStress:
          return ((float(E)/(1.-nu**2))*
          np.array([[1 , nu, 0     ],
                    [nu, 1 , 0     ],
                    [0 , 0 , 1-nu]]));
        else:
          return  ((float(E)/((1+nu)*(1-2*nu)))*
          np.array([[1-nu, nu  , 0 ],
                    [nu  , 1-nu, 0 ],
                    [0   , 0   , 0.5-nu]]));

    res= ((float(E)/((1+nu)*(1-2*nu)))*
        np.array([[1-nu, nu  ,  nu  , 0      ,0      ,0 ],
                  [nu  , 1-nu,  nu  , 0      ,0      ,0 ],
                  [nu  , nu  ,  1-nu, 0      ,0      ,0 ],
                  [0   , 0   ,  0   , 0.5-nu ,0      ,0 ],
                  [0   , 0   ,  0   , 0      ,0.5-nu ,0 ],
                  [0   , 0   ,  0   , 0      ,0      ,0.5-nu]]));
    return res


def LaplaceOrtho(k1,k2,k3=1, dim = 3):
    if dim == 2:
        return np.array([[k1 , 0 ],
                         [0  , k2]]);

    res= np.array([[k1 , 0 ,  0  ],
                   [0  , k2,  0  ],
                   [0  ,0  ,  k3 ]]);
    return res

def CheckIntegrity():
    HI3D = HookeIso(1,0.3)
    HI2D = HookeIso(1,0.3,dim=2)
    HI2D = HookeIso(1,0.3,dim=2,PlaneStress= False)

    LO3D = LaplaceOrtho(1,2,3)
    LO2D= LaplaceOrtho(1,2,dim=2)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover