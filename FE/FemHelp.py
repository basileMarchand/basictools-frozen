# -*- coding: utf-8 -*-

import math
import numpy as np
    
def Integral(E,Bop,elem,ndofs):
  
    res = np.zeros((ndofs,ndofs),float)    
    
    ## from http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
    ## 1 gauss point
    #p = [0. ];
    #w = [2.]
    ## 2 gauss points
    p = [-math.sqrt(1./3), math.sqrt(1./3) ];
    w = [1., 1.]
    ## 3 gauss points
    #p = [-math.sqrt(3./5.),0, math.sqrt(3./5.) ];
    #w = [5./9., 8./9.,5./9.]
    ## 4 gauss points
    #w1 =(1./2.)-math.sqrt(5./6.)/6.
    #w2 = (1./2.)+math.sqrt(5./6.)/6.
    #w = [w1,w2 , w2, w1];
    #p1 = math.sqrt((3+2*math.sqrt(6./5.))/7.);
    #p2 = math.sqrt((3-2*math.sqrt(6./5.))/7.);
    #p=[-p1,-p2 , p2, p1]
 
    for i in range(len(p)):
        for j in range(len(p)):
            for k in range(len(p)):
                qcoor = [p[i], p[j], p[k]]
                wp =  w[i]*w[j]*w[k]
                B = Bop(qcoor);
                Jdet = elem.GetDetJack(qcoor);
                res = res + wp*Jdet*B.T.dot(E.dot(B))
    return res
    
def CheckIntegrity():
    
    return 'ok'
    
    
if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover

