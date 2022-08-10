# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import NewType, Dict, Tuple, List
import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType, ArrayLike
import BasicTools.Containers.ElementNames as EN

ElementIntegrationRuleType = NewType('ElementIntegrationRuleType', Tuple[np.ndarray, np.ndarray])
IntegrationRulesType = NewType('IntegrationRulesType', Dict[str,ElementIntegrationRuleType] )

IntegrationRulesAlmanac = {} # type: Dict[str,IntegrationRulesType]

def GetRule(ruleName:str=None, rule:ElementIntegrationRuleType=None) -> ElementIntegrationRuleType:
    """Get the integration rule using a name or a integration rule. In the case where
    both are none then LagrangeIsoParam is returned

    Parameters
    ----------
    ruleName : str, optional
        name of the integration rule, by default None
    rule : IntegrationRuleType, optional
        integration rule in the case ruleName is None, by default None

    Returns
    -------
    IntegrationRuleType
        the integration rule
    """
    if ruleName is None :
        if rule is None:
            return IntegrationRulesAlmanac["LagrangeIsoParam"]
        else:
            return rule
    else:
        if rule is None:
            return IntegrationRulesAlmanac[ruleName]
        else:
            raise(Exception("must give ruleName or rule not both"))

def TensorProductGauss(dim:int, npoints:int=2)-> ElementIntegrationRuleType:
    """Tensor product of gauss point in dim dimension using npoints per dimension

    Parameters
    ----------
    dim : int
        the dimensionality of the tensor product
    npoints : int, optional
        Number of points per dimension , by default 2

    Returns
    -------
    IntegrationRuleType
        Integration Rule (points and weights)
    """
    import math
    if npoints == 2:
        p = [-math.sqrt(1./3)/2+.5, math.sqrt(1./3)/2+.5 ]
        w = [1./2., 1./2.]
    elif npoints == 3:
        #https://fr.wikipedia.org/wiki/M%C3%A9thodes_de_quadrature_de_Gauss
        p = [-math.sqrt(3./5.)/2+.5, 0.5 ,math.sqrt(3./5.)/2+.5 ]
        w = [5./18., 8./18., 5./18.]
    else:
        raise # pragma: no cover

    if dim == 1:
        return ElementIntegrationRuleType( (np.array([p,]).T,np.array(w)) )

    return TensorProdHomogeneous(dim,np.array([p]).T,np.array(w))

def TensorProdHomogeneous(dim:int,p:ArrayLike,w:ArrayLike) ->ElementIntegrationRuleType:
    #pp,ww = TensorProd(p,w,p,w)
    #for i in range(dim-2):
    #    pp,ww = TensorProd(p,w,pp,ww)
    #return pp,ww
    if dim == 2:
        return TensorProd(p,w,p,w)
    elif dim ==3:
        return TensorProd(p,w,p,w,p,w)
    else :
        raise# pragma: no cover

def TensorProd(p1:ArrayLike, w1:ArrayLike,
               p2:ArrayLike, w2:ArrayLike,
               p3:ArrayLike=None, w3:ArrayLike=None) -> ElementIntegrationRuleType:
    Pres = []
    Wres = []
    if p3 is None:
        for i in range(len(p1)):
            for j in range(len(p2)):
                    res: List[PBasicFloatType] = []
                    res.extend(p1[i,:])
                    res.extend(p2[j,:])
                    Pres.append(res)

                    #Pres.append([p1[i,0], p2[j,0]])
                    Wres.append( w1[i]*w2[j])
        return ElementIntegrationRuleType((np.array(Pres),np.array(Wres)))
    else:
        for k in range(len(p3)):
            for j in range(len(p2)):
                for i in range(len(p1)):
                    res = []
                    res.extend(p1[i,:])
                    res.extend(p2[j,:])
                    res.extend(p3[k,:])
                    Pres.append(res)

                    #Pres.append([p1[i,0], p2[j,0]])
                    Wres.append( w1[i]*w2[j]*w3[k])
        return ElementIntegrationRuleType( (np.array(Pres),np.array(Wres)) )

        #p23,w23 = TensorProd(p2,w2,p3,w3)
        #return TensorProd(p1,w1,p23,w23)

### integration Point in the center of the element #####
ElementCenter = IntegrationRulesType({})

IntegrationRulesAlmanac["ElementCenterEval"] = ElementCenter

# 1D elements
ElementCenter[EN.Bar_2] = ( 1./2.*np.array([[1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Bar_3] = ElementCenter[EN.Bar_2]

# 2D elements
ElementCenter[EN.Triangle_3] = ( 1./3.*np.array([[1., 1.] ]),
                                  np.array([0.5 ]))
ElementCenter[EN.Triangle_6] = ElementCenter[EN.Triangle_3]

ElementCenter[EN.Quadrangle_4] = ( 1./2.*np.array([[1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Quadrangle_8] = ElementCenter[EN.Quadrangle_4]
ElementCenter[EN.Quadrangle_9] = ElementCenter[EN.Quadrangle_4]

# 3D elements
ElementCenter[EN.Tetrahedron_4] = ( 1./4.*np.array([[1., 1.,1.] ]),
                                  np.array([ 1./6. ]))
ElementCenter[EN.Tetrahedron_10] = ElementCenter[EN.Tetrahedron_4]

ElementCenter[EN.Hexaedron_8] = ( 1./2.*np.array([[1., 1., 1.] ]),
                                  np.array([1. ]))
ElementCenter[EN.Hexaedron_20] = ElementCenter[EN.Hexaedron_8]
ElementCenter[EN.Hexaedron_27] = ElementCenter[EN.Hexaedron_8]
ElementCenter[EN.Wedge_6] = ( np.array([[1./3, 1./3, 0.5] ]),
                                  np.array([0.5 ]))

################# Lagrange P1 ##################################
LagrangeP1 = {}
IntegrationRulesAlmanac["LagrangeP1"] = LagrangeP1

LagrangeP1[EN.Point_1] = ( np.array([[0] ]),
                          np.array([1.]))

LagrangeP1[EN.Triangle_3] = ( 1./6.*np.array([[1., 1.] ,[4., 1.],[ 1. ,4.] ]),
                                  1./6.*np.array([1. , 1. , 1]))
LagrangeP1[EN.Triangle_6] = LagrangeP1[EN.Triangle_3]

LagrangeP1[EN.Tetrahedron_4] = ( np.array([[1/4,1/4,1/4]]), np.array([ 1./6.]))
LagrangeP1[EN.Tetrahedron_10] = LagrangeP1[EN.Tetrahedron_4]
LagrangeP1[EN.Bar_2] = TensorProductGauss(dim=1,npoints=2)
LagrangeP1[EN.Bar_3] = LagrangeP1[EN.Bar_2]
LagrangeP1[EN.Quadrangle_4] = TensorProductGauss(dim=2,npoints=2)
LagrangeP1[EN.Quadrangle_8] = LagrangeP1[EN.Quadrangle_4]
LagrangeP1[EN.Quadrangle_9] = LagrangeP1[EN.Quadrangle_4]
LagrangeP1[EN.Hexaedron_8]  = TensorProductGauss(dim=3,npoints=2)
LagrangeP1[EN.Hexaedron_20] = LagrangeP1[EN.Hexaedron_8]
LagrangeP1[EN.Hexaedron_27] = LagrangeP1[EN.Hexaedron_8]
LagrangeP1[EN.Wedge_6] = TensorProd(LagrangeP1[EN.Triangle_3][0][:,0:2],LagrangeP1[EN.Triangle_3][1],
                                       LagrangeP1[EN.Bar_2][0],LagrangeP1[EN.Bar_2][1])

################# Lagrange P2 ##################################
LagrangeP2 = {}
IntegrationRulesAlmanac["LagrangeP2"] = LagrangeP2

LagrangeP2[EN.Point_1] = ( np.array([[0] ]), np.array([1.]))
LagrangeP2[EN.Bar_2]  = TensorProductGauss(dim=1,npoints=3)
LagrangeP2[EN.Bar_3] = LagrangeP2[EN.Bar_2]
LagrangeP2[EN.Quadrangle_4] = TensorProductGauss(dim=2,npoints=3)
LagrangeP2[EN.Quadrangle_8] = LagrangeP2[EN.Quadrangle_4]
LagrangeP2[EN.Quadrangle_9] = LagrangeP2[EN.Quadrangle_4]
LagrangeP2[EN.Hexaedron_8]  = TensorProductGauss(dim=3,npoints=3)
LagrangeP2[EN.Hexaedron_20] = LagrangeP2[EN.Hexaedron_8]
LagrangeP2[EN.Hexaedron_27] = LagrangeP2[EN.Hexaedron_8]

#from https://www.code-aster.org/V2/doc/v13/en/man_r/r3/r3.01.01.pdf
p0 = 1./4.
p1 = 1./6.
p2 = 1./2.
p=np.array([[p1,p1,p1],
            [p2,p1,p1],
            [p1,p2,p1],
            [p1,p1,p2],
            [p0,p0,p0]])
w0 = -2./15
w1 = 3./40
w= np.array([w1, w1, w1, w1, w0])
LagrangeP2[EN.Tetrahedron_4] = (p,w)
LagrangeP2[EN.Tetrahedron_10] = LagrangeP2[EN.Tetrahedron_4]

p = np.array([[0.445948490915,0.445948490915,0],
[0.108103018168,0.445948490915,0],
[0.445948490915,0.108103018168,0],
[0.091576213509,0.091576213509,0],
[0.816847572980,0.091576213509,0],
[0.091576213509,0.816847572980,0]])
w = np.array([0.1116907948390 ,0.1116907948390 ,0.1116907948390 ,0.0549758718276 ,0.0549758718276 ,0.0549758718276])

LagrangeP2[EN.Triangle_3] = (p,w)
LagrangeP2[EN.Triangle_6] = LagrangeP2[EN.Triangle_3]
LagrangeP2[EN.Wedge_6] = TensorProd(LagrangeP2[EN.Triangle_3][0][:,0:2],LagrangeP2[EN.Triangle_3][1],
                                       LagrangeP2[EN.Bar_2][0],LagrangeP2[EN.Bar_2][1])

LagrangeIsoParam  = {}
for name in LagrangeP1:
    if EN.degree[name] == 1:
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

import BasicTools.FE.Spaces.WedgeSpaces as WS
wedge = WS.Wedge_P1_Lagrange()
NodalEvaluationP1[EN.Wedge_6] = (wedge.posN , np.ones(wedge.posN.shape[0]) )

##### Nodal P2 Itegration points for the evaluation of post quantities at nodes  ######
NodalEvaluationP2 = {}
IntegrationRulesAlmanac["NodalEvalP2"] = NodalEvaluationP2

import BasicTools.FE.Spaces.BarSpaces as BaS
bar = BaS.Bar_P2_Lagrange()
NodalEvaluationP2[EN.Bar_2] = (bar.posN , np.ones(tri.posN.shape[0]) )
NodalEvaluationP2[EN.Bar_3] = NodalEvaluationP2[EN.Bar_2]

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

import BasicTools.FE.Spaces.WedgeSpaces as WS
wedgep2 = WS.Wedge_P2_Lagrange()
NodalEvaluationP2[EN.Wedge_6] = (wedgep2.posN , np.ones(wedgep2.posN.shape[0]) )



import BasicTools.FE.Spaces.QuadSpaces as QS
quad = QS.Quad_P2_Lagrange()
NodalEvaluationP2[EN.Quadrangle_4] = (quad.posN , np.ones(quad.posN.shape[0]) )
NodalEvaluationP2[EN.Quadrangle_8] = NodalEvaluationP2[EN.Quadrangle_4]
NodalEvaluationP2[EN.Quadrangle_9] = NodalEvaluationP2[EN.Quadrangle_4]


NodalEvalIsoGeo  = {}
IntegrationRulesAlmanac["NodalEvalGeo"] = NodalEvalIsoGeo
for name in NodalEvaluationP1:
    if EN.degree[name] == 1:
        NodalEvalIsoGeo[name] = NodalEvaluationP1[name]
    else:
        NodalEvalIsoGeo[name] = NodalEvaluationP2[name]

trapezoidalOrderGenerated = 3
# trapezoidal Rules
for i in range(1,trapezoidalOrderGenerated):
    traprule = {}
    IntegrationRulesAlmanac["TrapezoidalP"+str(i)] = traprule
    if i == 1:
        w = np.ones(1)
        p = np.array([[0.5]])
    elif i == 2:
        w = np.ones(2)/2
        p = np.array([[0, 1]]).T
    else:
        w = np.ones(i)/(2*i-2)
        w[1:-1] *=2
        p = np.arange(i)[:,np.newaxis]/(i-1)
    traprule[EN.Bar_2] = ( p, w)
    traprule[EN.Bar_3] = ( p , w)
    traprule[EN.Quadrangle_4] = TensorProdHomogeneous(2,p,w)
    traprule[EN.Quadrangle_8] = TensorProdHomogeneous(2,p,w)
    traprule[EN.Quadrangle_9] = TensorProdHomogeneous(2,p,w)
    traprule[EN.Hexaedron_8]  = TensorProdHomogeneous(3,p,w)
    traprule[EN.Hexaedron_20] = TensorProdHomogeneous(3,p,w)
    traprule[EN.Hexaedron_27] = TensorProdHomogeneous(3,p,w)

for i in range(1,trapezoidalOrderGenerated):
    traprule = IntegrationRulesAlmanac["TrapezoidalP"+str(i)]
    if i == 1:
        p = 1./2.*np.array([[1,1.],])
        w = 1./2.*np.array([1. ])
    else:
        p = []
        w = np.ones(i*(i+1)//2)
        cpt = 0
        for j in range(i):
            phi0 = 1./(i-1)*j
            for k in range(i-j):
                phi1 = 1./(i-1)*(k)
                p.append([phi0,phi1])
                #halp of the contribution of the border points
                if j == 0  or  j == i-1 or  k == 0 or  k == (i-j-1):
                    w[cpt] *= 0.5
                # 1/6 for the corner points
                if (j == 0 and k == 0) or (j == 0 and k == i-j-1) or ( j == i-1 and k == 0):
                    w[cpt] = 1./6.
                cpt += 1

        w = np.array(w)
        # normalization of the weight to sum = 1/2
        w *= 1/(np.sum(w)*2)
        p = np.array(p)
    traprule[EN.Triangle_3] = (p,w)
    traprule[EN.Triangle_6] = (p,w)
    traprule[EN.Wedge_6] = TensorProd(p[:,0:2],w,traprule[EN.Bar_2][0],traprule[EN.Bar_2][1])

def CheckIntegrity(GUI=False):
    GetRule()
    GetRule(ruleName="TrapezoidalP1")
    GetRule(rule=NodalEvalIsoGeo)
    error = False
    #must fail
    try:
        GetRule(ruleName="TrapezoidalP1",rule=NodalEvalIsoGeo)
        error = True# pragma: no cover
    except:
        pass
    if error: raise

    VolumeOne = [EN.Point_1,
                 EN.Bar_2, EN.Bar_3,
                 EN.Quadrangle_4, EN.Quadrangle_8,EN.Quadrangle_9,
                 EN.Hexaedron_8, EN.Hexaedron_20, EN.Hexaedron_27]
    VolumeHalf = [EN.Triangle_3,EN.Triangle_6, EN.Wedge_6, EN.Wedge_15, EN.Wedge_18]
    VolumeSixth = [EN.Tetrahedron_4, EN.Tetrahedron_10]

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

    for rulename,rule  in IntegrationRulesAlmanac.items():
        print("---" + rulename + "----------------------------------------------")
        for elemName,data in rule.items():
            #if elemName is not EN.Triangle_3 :
            #    continue
            p = data[0]
            w = data[1]
            lp = len(data[0] )
            lw = len(data[1] )
            print(rulename + " " + elemName + " has " + str(lw) + " integration points")
            if lp != lw:# pragma: no cover
                print(data)
                raise( Exception("incompatible rule") )

            #Nodal... not designed for integration
            if rulename in [ "NodalEvalP1" ,"NodalEvalP2", "NodalEvalGeo"] :
                continue

            if elemName in VolumeOne:
                volref = 1.
            elif elemName in VolumeHalf:
                volref = 0.5
            else:
                volref = 1./6.

            if np.abs(np.sum(data[1])- volref)  >  1e-10:# pragma: no cover
                print(data)
                print("Mesure : ",np.sum(data[1]))
                print("Mesure Error: ",np.abs(np.sum(data[1])- volref))
                raise(Exception(rulename + " " + elemName + " does not itegrate constast funciton"))

            #if rulename in ["ElementCenterEval"]:
            #    continue

            # End check of constant functiton

            # Checking linear function

            #cant integrate over a point
            if elemName in [EN.Point_1] :
                continue

            # TrapezoidalP1 cant linear funtions
            if rulename in ['TrapezoidalP1']:
               continue

            def f1(x):
                return x[0]-0.5

            if elemName in VolumeOne:
                volref = 0.
            elif elemName in VolumeHalf:
                volref = (1/3-0.5)/2
            else:
                volref = -0.25/6

            space_ipvalues = LagrangeSpaceGeo[elemName].SetIntegrationRule(p,w)
            integral = 0
            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,LagrangeSpaceGeo[elemName].posN)
                integral +=  Jdet*w[ip]*f1(p[ip])

            #integral = data[1].dot([f1(x) for x in data[0]])
            if np.abs(integral - volref)  >  1e-10:# pragma: no cover
                print("function :  f(x) = x-0.5")
                print("int S f(x)dx =  {}".format(integral) )
                print("Integral Exact: ",volref)
                print("Integral Error: ",np.abs(integral- volref))
                #print(data)
                raise(Exception(rulename + " " + elemName + " does not itegrate f(x) = x-0.5"))

            # ElementCenterEval and LagrangeP1 cant integrate quadratic funtions
            if rulename in ['ElementCenterEval', "LagrangeP1"] or rulename[0:len("Trapezoidal")] == "Trapezoidal" :
                continue

            # this is a degre 1 integration
            if rulename == "LagrangeIsoParam" and elemName == 'tet4':
                continue

            # Checking quadratic integration
            def f2(x):
                return (x[0]-0.5)**2

            if elemName in VolumeOne:
                volref = 1./12.
            elif elemName in VolumeHalf:
                volref = 1/24.
            elif elemName in VolumeSixth:
                volref = 1./60.
            else:
                raise # pragma: no cover

            space_ipvalues = LagrangeSpaceGeo[elemName].SetIntegrationRule(p,w)
            integral = 0
            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,LagrangeSpaceGeo[elemName].posN)
                integral +=  Jdet*w[ip]*f2(p[ip])

            #integral = data[1].dot([f1(x) for x in data[0]])
            if np.abs(integral - volref)  >  1e-10:# pragma: no cover
                print("function :  f(x) = (x-0.5)**2")
                print("int S f(x)dx =  {}".format(integral) )
                print("Integral Exact: ",volref)
                print("Integral Error: ",np.abs(integral- volref))
                #print(data)
                raise(Exception(rulename + " " + elemName + " does not itegrate f(x) = (x-0.5)**2"))




    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover

