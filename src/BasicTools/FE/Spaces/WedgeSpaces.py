# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np
from sympy.matrices import Matrix

import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase


class WedgeSpaceBase(SymSpaceBase):
    def __init__(self):
        super(WedgeSpaceBase,self).__init__()
        self.dimensionality = 3
        self.geoSupport = EN.GeoWed
    def ClampParamCoorninates(self,xietaphi):
        res = xietaphi.copy()
        if res[0]+res[1] > 1:
            dif = res[0]-res[1]
            res[0] = (1+dif)/2.
            res[1] = 1.-res[0]
        res[0] = max(min(res[0],1),0)
        res[1] = max(min(res[1],1),0)
        res[2] = max(min(res[2],1),0)
        return res

class Wedge_P0_Global(WedgeSpaceBase):
    def __init__(self):
        super(Wedge_P0_Global,self).__init__()

        self.symN = Matrix([1])
        self.posN = np.array([ [ None, None, None] ])
        self.dofAttachments = [("G",None,None)]

class Wedge_P0_Lagrange(WedgeSpaceBase):
    def __init__(self):
        super(Wedge_P0_Lagrange,self).__init__()

        self.symN = Matrix([1])
        self.posN = np.array([ [ 1./3, 1./3, 0.5] ])
        self.dofAttachments = [("C",0,None)]

class Wedge_P1_Lagrange(WedgeSpaceBase):
    def __init__(self):
        super(Wedge_P1_Lagrange,self).__init__()
        xi = self.xi
        eta = self.eta
        phi = self.phi
        self.symN = Matrix([(1-xi-eta)*(1-phi),
                            xi*(1-phi),
                            eta*(1-phi),
                            (1-xi-eta)*phi,
                            xi*phi,
                            eta*phi
                            ])
        self.posN = np.array([[0,0,0],
                              [1,0,0],
                              [0,1,0],
                              [0,0,1],
                              [1,0,1],
                              [0,1,1]])

        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ("P",4,None),
                               ("P",5,None) ]


class Wedge_P2_Lagrange(WedgeSpaceBase):
    def __init__(self):
        super(Wedge_P2_Lagrange,self).__init__()
        xi = self.xi
        eta = self.eta
        phi = self.phi
        T = (1-xi-eta)

        trisf = Matrix([T*(1-2*xi-2*eta),
                            xi *(2*xi-1),
                            eta*(2*eta-1),
                            4*xi*T,
                            4*xi*eta,
                            4*eta*T])



        L1 = 1-phi
        L2 = phi
        barsymN = Matrix([L1*(2*L1-1),L2*(2*L2-1),4*L1*L2 ])


        self.posN = np.array([[0,0,0],
                             [1,0,0],
                             [0,1,0],
                             [0,0,1],
                             [1,0,1],
                             [0,1,1],

                            [0.5,  0,0],
                            [0.5,0.5,0],
                            [0  ,0.5,0],

                            [0.5,  0,1],
                            [0.5,0.5,1],
                            [0  ,0.5,1],

                            [0,0,0.5],
                            [1,0,0.5],
                            [0,1,0.5],

                            [0.5,  0,0.5],
                            [0.5,0.5,0.5],
                            [0  ,0.5,0.5],

                             ])

        self.symN = Matrix([trisf[0]*barsymN[0], #P 0
                            trisf[1]*barsymN[0], #P 1
                            trisf[2]*barsymN[0], #P 2
                            trisf[0]*barsymN[1], #P 3
                            trisf[1]*barsymN[1], #P 4
                            trisf[2]*barsymN[1], #P 5

                            trisf[3]*barsymN[0], #F2 0
                            trisf[4]*barsymN[0], #F2 1
                            trisf[5]*barsymN[0], #F2 2
                            trisf[3]*barsymN[1], #F2 6
                            trisf[4]*barsymN[1], #F2 7
                            trisf[5]*barsymN[1], #F2 8

                            trisf[0]*barsymN[2], #F2 3
                            trisf[1]*barsymN[2], #F2 4
                            trisf[2]*barsymN[2], #F2 5
                            trisf[3]*barsymN[2], #F 0
                            trisf[4]*barsymN[2], #F 1
                            trisf[5]*barsymN[2]] #F 2
                            )



        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ("P",4,None),
                               ("P",5,None),

                               ("F2",0,None),
                               ("F2",1,None),
                               ("F2",2,None),
                               ("F2",6,None),
                               ("F2",7,None),
                               ("F2",8,None),

                               ("F2",3,None),
                               ("F2",4,None),
                               ("F2",5,None),
                               ("F",2,None),
                               ("F",3,None),
                               ("F",4,None),

                                ]
