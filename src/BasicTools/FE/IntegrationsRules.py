# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

import BasicTools.Containers.ElementNames as EN

def TensorProductPoints(dim,npoints=2):
    import math
    if npoints == 2:
        p = [-math.sqrt(1./3)/2+.5, math.sqrt(1./3)/2+.5 ];
        w = [1./2., 1./2.]
    elif npoints == 3:
        #https://fr.wikipedia.org/wiki/M%C3%A9thodes_de_quadrature_de_Gauss
        p = [-math.sqrt(3./5.)/2+.5, 0.5 ,math.sqrt(3./5.)/2+.5 ];
        w = [5./18., 8./18., 5./18.]
    else:
        raise
    Pres = []
    Wres = []
    if dim ==1:
        return (np.array([p,]).T,np.array(w))
    if dim == 2:
        for i in range(len(p)):
            for j in range(len(p)):
                    Pres.append([p[i], p[j]])
                    Wres.append( w[i]*w[j])
        return (np.array(Pres),np.array(Wres))
    elif dim ==3:
        for i in range(len(p)):
            for j in range(len(p)):
                for k in range(len(p)):
                    Pres.append([p[k], p[j],p[i] ])
                    Wres.append( w[i]*w[j]*w[k])
        return (np.array(Pres),np.array(Wres))
    else :
        raise


IntegrationRulesAlmanac = {}

### integration Point in the center of the elemetn #####
ElementCenter = {}
IntegrationRulesAlmanac["ElementCenterEval"] = ElementCenter

# 1D elements
ElementCenter[EN.Bar_2] = ( 1./2.*np.array([[1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Bar_3] = ElementCenter[EN.Bar_2]

# 2D elements
ElementCenter[EN.Triangle_3] = ( 1./3.*np.array([[1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Triangle_6] = ElementCenter[EN.Triangle_3]

ElementCenter[EN.Quadrangle_4] = ( 1./2.*np.array([[1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Quadrangle_8] = ElementCenter[EN.Quadrangle_4]
ElementCenter[EN.Quadrangle_9] = ElementCenter[EN.Quadrangle_4]

# 3D elements
ElementCenter[EN.Tetrahedron_4] = ( 1./4.*np.array([[1., 1.,1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Tetrahedron_10] = ElementCenter[EN.Tetrahedron_4]

ElementCenter[EN.Hexaedron_8] = ( 1./2.*np.array([[1., 1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Hexaedron_20] = ElementCenter[EN.Hexaedron_8]
ElementCenter[EN.Hexaedron_27] = ElementCenter[EN.Hexaedron_8]


################# Lagrange P1 ##################################
LagrangeP1 = {}
IntegrationRulesAlmanac["LagrangeP1"] = LagrangeP1

LagrangeP1[EN.Point_1] = ( np.array([[0] ]),
                          np.array([1.]))

LagrangeP1[EN.Triangle_3] = ( 1./6.*np.array([[1., 1.] ,[4., 1.],[ 1. ,4.] ]),
                                  1./6.*np.array([1. , 1. , 1]))
LagrangeP1[EN.Triangle_6] = LagrangeP1[EN.Triangle_3]

#from https://www.code-aster.org/V2/doc/v13/en/man_r/r3/r3.01.01.pdf
p0 = 1./4.
p1 = 1./6.
p2 = 1./2.
p=np.array([[p0,p0,p0],
            [p1,p1,p1],
            [p2,p1,p1],
            [p1,p2,p1],
            [p1,p1,p2]]);
w0 = -2./15
w1 = 3./40
w= np.array([ w0, w1, w1, w1, w1]);
LagrangeP1[EN.Tetrahedron_4] = ( p,w)
LagrangeP1[EN.Tetrahedron_10] = LagrangeP1[EN.Tetrahedron_4]
LagrangeP1[EN.Bar_2] = TensorProductPoints(dim=1,npoints=2)
LagrangeP1[EN.Bar_3] = LagrangeP1[EN.Bar_2]
LagrangeP1[EN.Quadrangle_4] = TensorProductPoints(dim=2,npoints=2)
LagrangeP1[EN.Quadrangle_8] = LagrangeP1[EN.Quadrangle_4]
LagrangeP1[EN.Quadrangle_9] = LagrangeP1[EN.Quadrangle_4]
LagrangeP1[EN.Hexaedron_8]  = TensorProductPoints(dim=3,npoints=2)
LagrangeP1[EN.Hexaedron_20] = LagrangeP1[EN.Hexaedron_8]
LagrangeP1[EN.Hexaedron_27] = LagrangeP1[EN.Hexaedron_8]

################# Lagrange P2 ##################################
LagrangeP2 = {}
IntegrationRulesAlmanac["LagrangeP2"] = LagrangeP2

LagrangeP2[EN.Point_1] = ( np.array([[0] ]), np.array([1.]))
LagrangeP2[EN.Bar_2]  = TensorProductPoints(dim=1,npoints=3)
LagrangeP2[EN.Bar_3] = LagrangeP2[EN.Bar_2]
LagrangeP2[EN.Quadrangle_4] = TensorProductPoints(dim=2,npoints=3)
LagrangeP2[EN.Quadrangle_8] = LagrangeP2[EN.Quadrangle_4]
LagrangeP2[EN.Quadrangle_9] = LagrangeP2[EN.Quadrangle_4]
LagrangeP2[EN.Hexaedron_8]  = TensorProductPoints(dim=3,npoints=3)
LagrangeP2[EN.Hexaedron_20] = LagrangeP2[EN.Hexaedron_8]
LagrangeP2[EN.Hexaedron_27] = LagrangeP2[EN.Hexaedron_8]

p = np.array([[0.1666666667,0.1666666667,0.1666666667],
[0.5,0.1666666667,0.1666666667],
[0.1666666667,0.5,0.1666666667],
[0.1666666667,0.1666666667,0.5],
[0.25,0.25,0.25]]);
w = np.array([0.075 ,0.075 ,0.075 ,0.075 ,-0.1333333333]);

LagrangeP2[EN.Tetrahedron_4] = (p,w)
LagrangeP2[EN.Tetrahedron_10] = LagrangeP2[EN.Tetrahedron_4]

p = np.array([[0.4459484909,0.4459484909,0],
[0.1081030182,0.4459484909,0],
[0.4459484909,0.1081030182,0],
[0.09157621351,0.09157621351,0],
[0.816847573,0.09157621351,0],
[0.09157621351,0.816847573,0]])
w = np.array([0.1116907948 ,0.1116907948 ,0.1116907948 ,0.05497587183 ,0.05497587183 ,0.05497587183])

LagrangeP2[EN.Triangle_3] = (p,w)
LagrangeP2[EN.Triangle_6] = LagrangeP2[EN.Triangle_3]





# vvvvvv this is not needed any more  vvvvv
#Define dictionnary for LagrangeP1 and LangrangeP2, based on "linear" dictionnary from ElementNames
linearLagrange = {}
linearLagrange[True]  = LagrangeP1
linearLagrange[False] = LagrangeP2

def Lagrange(elementName):
    return  LagrangeIsoParam[elementName]
    #return linearLagrange[EN.linear[elementName]][EN.geoSupport[elementName]]
# ^^^^^^ this is not needed any more ^^^^^^

LagrangeIsoParam  = {}
for name in LagrangeP1:
    if EN.linear[name]:
        LagrangeIsoParam[name] = LagrangeP1[name]
    else:
        LagrangeIsoParam[name] = LagrangeP2[name]


IntegrationRulesAlmanac["LagrangeIsoParam"] = LagrangeIsoParam

##### Nodal P1 Itegration points for the evaluation of post quantities at nodes  ######
NodalEvaluationP1 = {}
IntegrationRulesAlmanac["NodalEvalP1"] = NodalEvaluationP1

import BasicTools.FE.Spaces.TriSpaces as TrS
tri = TrS.Tri_P1_Lagrange()
NodalEvaluationP1[EN.Triangle_3] = (tri.posN , np.ones(tri.posN.shape[0]) )
NodalEvaluationP1[EN.Triangle_6] = NodalEvaluationP1[EN.Triangle_3]

import BasicTools.FE.Spaces.TetSpaces as TS
tet = TS.Tet_P1_Lagrange()
NodalEvaluationP1[EN.Tetrahedron_4] = (tet.posN , np.ones(tet.posN.shape[0]) )
NodalEvaluationP1[EN.Tetrahedron_10] = NodalEvaluationP1[EN.Tetrahedron_4]

import BasicTools.FE.Spaces.HexaSpaces as HS
hexa = HS.Hexa_P1_Lagrange()
NodalEvaluationP1[EN.Hexaedron_8] = (hexa.posN , np.ones(hexa.posN.shape[0]) )
#NodalEvaluationP1[EN.Hexaedron_20] = NodalEvaluationP1[EN.Hexaedron_8]
NodalEvaluationP1[EN.Hexaedron_27] = NodalEvaluationP1[EN.Hexaedron_8]

##### Nodal P2 Itegration points for the evaluation of post quantities at nodes  ######
NodalEvaluationP2 = {}
IntegrationRulesAlmanac["NodalEvalP2"] = NodalEvaluationP2

import BasicTools.FE.Spaces.TriSpaces as TrS
tri = TrS.Tri_P2_Lagrange()
NodalEvaluationP2[EN.Triangle_3] = (tri.posN , np.ones(tri.posN.shape[0]) )
NodalEvaluationP2[EN.Triangle_6] = NodalEvaluationP2[EN.Triangle_3]

import BasicTools.FE.Spaces.TetSpaces as TS
tet = TS.Tet_P2_Lagrange()
NodalEvaluationP2[EN.Tetrahedron_4] = (tet.posN , np.ones(tet.posN.shape[0]) )
NodalEvaluationP2[EN.Tetrahedron_10] = NodalEvaluationP2[EN.Tetrahedron_4]

import BasicTools.FE.Spaces.HexaSpaces as HS
hexa = HS.Hexa_P2_Lagrange()
NodalEvaluationP2[EN.Hexaedron_8] = (hexa.posN , np.ones(hexa.posN.shape[0]) )
#NodalEvaluationP2[EN.Hexaedron_20] = NodalEvaluationP2[EN.Hexaedron_8]
NodalEvaluationP2[EN.Hexaedron_27] = NodalEvaluationP2[EN.Hexaedron_8]


NodalEvalIsoGeo  = {}
IntegrationRulesAlmanac["NodalEvalGeo"] = NodalEvalIsoGeo
for name in NodalEvaluationP1:
    if EN.linear[name]:
        NodalEvalIsoGeo[name] = NodalEvaluationP1[name]
    else:
        NodalEvalIsoGeo[name] = NodalEvaluationP2[name]


def CheckIntegrity(GUI=False):
    for rulename,rule  in IntegrationRulesAlmanac.items():
        print(rulename)
        for elemName,data in rule.items():
            lw = len(data[0] )
            print(rulename + " " + elemName + " has " + str(lw) + " integration points")
            if lw != len(data[1]):
                raise( Exception("incompatible rule") )
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover

