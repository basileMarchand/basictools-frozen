# -*- coding: utf-8 -*-
import numpy as np

import BasicTools.FE.ElementNames as EN
# Triangle GaussLobato 3 poins



IntegrationRulesAlmanac = {}


LagrangeP1 = {}
LagrangeP1[EN.GeoTri] = ( 1./6.*np.array([[1, 1] ,[4, 1],[ 1 ,4] ]),
                                  1./6.*np.array([1 , 1 , 1]))

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
    else:
        raise
    Pres = []
    Wres = []
    if dim ==1:
        return (np.array([p,]).T,w)
    if dim == 2:
        for i in range(len(p)):
            for j in range(len(p)):
                    Pres.append([p[i], p[j]])
                    Wres.append( w[i]*w[j])
        return (np.array(Pres),Wres)
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
        return (np.array(Pres),Wres)
    else :
        raise


LagrangeP1[EN.GeoBar]  = TensorProductPoints(dim=1,npoints=2)
LagrangeP1[EN.GeoQuad] = TensorProductPoints(dim=2,npoints=2)
LagrangeP1[EN.GeoHex]  = TensorProductPoints(dim=3,npoints=2)

IntegrationRulesAlmanac["IsoGeo"] = LagrangeP1


#Nodal Itegration for the evaluation of post quantities at nodes
NodalEvaluation = {}
import BasicTools.FE.Spaces.HexaSpaces as HS
hexa = HS.Hexa_P1_Lagrange()
NodalEvaluation[EN.GeoHex] = (hexa.posN , list(np.ones(hexa.posN.shape[0])) )
import BasicTools.FE.Spaces.TetSpaces as TS
tet = TS.Tet_P1_Lagrange()
NodalEvaluation[EN.GeoTet] = (tet.posN , list(np.ones(tet.posN.shape[0])) )
IntegrationRulesAlmanac["NodalEvalGeo"] = NodalEvaluation

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

