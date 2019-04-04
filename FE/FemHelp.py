# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import math
import numpy as np

import BasicTools.Containers.ElementNames as ElementsNames


def Integral(E,Bop,elem,ndofs,pos=None):
    """
    Function to calculate the integral of E(u)*Bop*E(u) over an element:

    inputes:
        linear operators E : matrix line
        Bop operator : Bop(integration point coordinates in element space, the coordinates of the nodes of the element )
            must return the B operator and the determinant of the jacobian.
       elem: the element
       ndofs: the number of dofs (E.shape[0])
       pos, the positions of the nodes of the element (if needed by Bop).

    # for the moment this is not a nice implementation
    """


    res = np.zeros((ndofs,ndofs),float)

    ## from http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
    ## 1 gauss point
    #p = [0. ];
    #w = [2.]
    ## 2 gauss points

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
    if elem.name == ElementsNames.Triangle_3:
        p=1./6.*np.array([[1, 1] ,[4, 1],[ 1 ,4] ]);
        w=1./6.*np.array([1 , 1 , 1]);
        for i in range(len(w)):
           B,Jdet = Bop(p[i,:],pos);
           res = res + w[i]*Jdet*B.T.dot(E.dot(B))
        return res

    if elem.name == ElementsNames.Tetrahedron_4:
        p0 = 1./4.
        p1 = 1./6.
        p2 = 1./2.
        p=np.array([[p0,p0,p0],
                    [p1,p1,p1],
                    [p2,p0,p0],
                    [p0,p2,p0],
                    [p0,p0,p2]]);
        w0 = -2./15
        w1 = 3./40
        w= np.array([ w0, w1, w1, w1, w1]);
        for i in range(len(w)):
           B,Jdet = Bop(p[i,:],pos);
           res = res + w[i]*Jdet*B.T.dot(E.dot(B))
        return res

    p = [-math.sqrt(1./3), math.sqrt(1./3) ];
    w = [1., 1.]
    if elem.dimensionality == 2:
        for i in range(len(p)):
            for j in range(len(p)):
                    qcoor = [p[i], p[j]]
                    wp =  w[i]*w[j]
                    B,Jdet = Bop(qcoor,pos);
                    #Jdet = elem.GetDetJack(qcoor);
                    res = res + wp*Jdet*B.T.dot(E.dot(B))
        return res
    else:
        for i in range(len(p)):
            for j in range(len(p)):
                for k in range(len(p)):
                    qcoor = [p[i], p[j], p[k]]
                    wp =  w[i]*w[j]*w[k]
                    B,Jdet = Bop(qcoor,pos);
                    Jdet = elem.GetDetJack(qcoor);
                    res = res + wp*Jdet*B.T.dot(E.dot(B))
        return res

def CheckIntegrity():
    from BasicTools.FE.Hexa8Cuboid import Hexa8Cuboid
    from BasicTools.FE.Quad4Rectangle import  Quad4Rectangle
    from BasicTools.FE.Tri3 import  Tri3

    Hexa8Cuboid().GetIsotropDispK(1.,0.3);
    Quad4Rectangle().GetIsotropDispK(1.,0.3);
    Tri3().GetIsotropDispK(1.,0.3,np.array([ [0,0], [1,0] , [0,1] ]));
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
