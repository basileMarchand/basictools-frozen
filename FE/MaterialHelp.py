# -*- coding: utf-8 -*-
import numpy as np

def HookeIso(E,nu):
    res= ((float(E)/((1+nu)*(1-2*nu)))*
        np.array([[1-nu, nu  ,  nu  , 0      ,0      ,0 ],
                  [nu  , 1-nu,  nu  , 0      ,0      ,0 ],
                  [nu  , nu  ,  1-nu, 0      ,0      ,0 ],
                  [0   , 0   ,  0   , 0.5-nu ,0      ,0 ],
                  [0   , 0   ,  0   , 0      ,0.5-nu ,0 ],
                  [0   , 0   ,  0   , 0      ,0      ,0.5-nu]]));
    return res
    
    
def LaplaceOrtho(k1,k2,k3):
    res= np.array([[k1 , 0 ,  0  ],
                   [0  , k2,  0  ],
                   [0  ,0  ,  k3 ]]);
    return res