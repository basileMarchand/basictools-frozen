# -*- coding: utf-8 -*-
import numpy as np

import BasicTools.Containers.ElementNames as EN

# Triangle GaussLobato 3 poins

IntegrationRulesAlmanac = {}

ElementCenter = {}

ElementCenter[EN.GeoHex] = ( 1./2.*np.array([[1., 1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.GeoQuad] = ( 1./2.*np.array([[1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.GeoTet] = ( 1./4.*np.array([[1., 1.,1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.GeoTri] = ( 1./3.*np.array([[1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.GeoBar] = ( 1./2.*np.array([[1.] ]),
                                  np.array([1. ]))



LagrangeP1 = {}

LagrangeP1[EN.GeoPoint] = ( np.array([[0] ]),
                          np.array([1.]))

LagrangeP1[EN.GeoTri] = ( 1./6.*np.array([[1., 1.] ,[4., 1.],[ 1. ,4.] ]),
                                  1./6.*np.array([1. , 1. , 1]))


#LagrangeP1[EN.GeoTri] = ( 1./3.*np.array([[1, 1]  ]),
#                                        np.array([ 0.5 ]))

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
LagrangeP1[EN.GeoTet] = ( p,w)


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

        """a = np.sqrt(1./3.)
        Pres = np.array([[-a,-a,-a],
                                          [ a,-a,-a],
                                          [-a, a,-a],
                                          [ a, a,-a],
                                          [-a,-a, a],
                                          [ a,-a, a],
                                          [-a, a, a],
                                          [ a, a, a]]);
        Wres = np.ones(8)/8.;"""
        return (np.array(Pres),np.array(Wres))
    else :
        raise


LagrangeP1[EN.GeoBar]  = TensorProductPoints(dim=1,npoints=2)
LagrangeP1[EN.GeoQuad] = TensorProductPoints(dim=2,npoints=2)
LagrangeP1[EN.GeoHex]  = TensorProductPoints(dim=3,npoints=2)


LagrangeP2 = {}

LagrangeP2[EN.GeoPoint] = ( np.array([[0] ]), np.array([1.]))
LagrangeP2[EN.GeoBar]  = TensorProductPoints(dim=1,npoints=3)
LagrangeP2[EN.GeoQuad] = TensorProductPoints(dim=2,npoints=3)
LagrangeP2[EN.GeoHex]  = TensorProductPoints(dim=3,npoints=3)


p = np.array([[0.1666666667,0.1666666667,0.1666666667],
[0.5,0.1666666667,0.1666666667],
[0.1666666667,0.5,0.1666666667],
[0.1666666667,0.1666666667,0.5],
[0.25,0.25,0.25]]);
w = np.array([0.075 ,0.075 ,0.075 ,0.075 ,-0.1333333333]);

LagrangeP2[EN.GeoTet] = (p,w)



p = np.array([[0.4459484909,0.4459484909,0],
[0.1081030182,0.4459484909,0],
[0.4459484909,0.1081030182,0],
[0.09157621351,0.09157621351,0],
[0.816847573,0.09157621351,0],
[0.09157621351,0.816847573,0]])
w = np.array([0.1116907948 ,0.1116907948 ,0.1116907948 ,0.05497587183 ,0.05497587183 ,0.05497587183])

LagrangeP2[EN.GeoTri] = (p,w)

IntegrationRulesAlmanac["LagrangeP1"] = LagrangeP1
IntegrationRulesAlmanac["LagrangeP2"] = LagrangeP2
IntegrationRulesAlmanac["IsoParam"] = LagrangeP1
IntegrationRulesAlmanac["ElementEvalGeo"] = ElementCenter
#LagrangeIsoParam = {} # Lagrange iso parametrique
#
#for name, linear in EN.linear.items():
#    if linear:
#        LagrangeIsoParam[name] = LagrangeP1[name]
#    else:
#        LagrangeIsoParam[name] = LagrangeP2[name]
#


#Nodal Itegration for the evaluation of post quantities at nodes
NodalEvaluation = {}
import BasicTools.FE.Spaces.HexaSpaces as HS
hexa = HS.Hexa_P1_Lagrange()
NodalEvaluation[EN.GeoHex] = (hexa.posN , np.ones(hexa.posN.shape[0]) )
import BasicTools.FE.Spaces.TetSpaces as TS
tet = TS.Tet_P1_Lagrange()
NodalEvaluation[EN.GeoTet] = (tet.posN , np.ones(tet.posN.shape[0]) )
import BasicTools.FE.Spaces.TriSpaces as TrS
tri = TrS.Tri_P1_Lagrange()
NodalEvaluation[EN.GeoTri] = (tri.posN , np.ones(tri.posN.shape[0]) )

IntegrationRulesAlmanac["NodalEvalGeo"] = NodalEvaluation
IntegrationRulesAlmanac["NodalEvalP1"] = NodalEvaluation
#Nodal Itegration for the evaluation of post quantities at nodes
NodalEvaluation = {}
import BasicTools.FE.Spaces.HexaSpaces as HS
hexa = HS.Hexa_P2_Lagrange()
NodalEvaluation[EN.GeoHex] = (hexa.posN , np.ones(hexa.posN.shape[0]) )
import BasicTools.FE.Spaces.TetSpaces as TS
tet = TS.Tet_P2_Lagrange()
NodalEvaluation[EN.GeoTet] = (tet.posN , np.ones(tet.posN.shape[0]) )
import BasicTools.FE.Spaces.TriSpaces as TrS
tri = TrS.Tri_P2_Lagrange()
NodalEvaluation[EN.GeoTri] = (tri.posN , np.ones(tri.posN.shape[0]) )

IntegrationRulesAlmanac["NodalEvalP2"] =NodalEvaluation

#Define dictionnary for LagrangeP1 and LangrangeP2, based on "linear" dictionnary from ElementNames
linearLagrange = {}
linearLagrange[True]  = LagrangeP1
linearLagrange[False] = LagrangeP2

def Lagrange(elementName):
  return linearLagrange[EN.linear[elementName]][EN.geoSupport[elementName]]

def CheckIntegrity(GUI=False):
    for rulename,rule  in IntegrationRulesAlmanac.items():
        print(rulename)
        for geoname,data in rule.items():
            lw = len(data[0] )
            print(rulename + " " + geoname.name + " has " + str(lw) + " integration points")
            if lw != len(data[1]):
                raise( Exception("incompatible rule") )
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover

