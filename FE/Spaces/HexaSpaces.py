# -*- coding: utf-8 -*-

import numpy as np
from sympy.matrices import Matrix

import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase


class Hexa_P0_Global(SymSpaceBase):
    def __init__(self):
        super(Hexa_P0_Global,self).__init__()
        self.geoSupport = EN.GeoQuad


        self.symN = Matrix([1])
        self.posN = np.array([ [ None, None, None] ])
        self.dofAttachments = [("G",None,None)]
        self.Create()


class Hexa_P0_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Hexa_P0_Lagrange,self).__init__()
        self.geoSupport = EN.GeoQuad


        self.symN = Matrix([1])
        self.posN = np.array([ [ 0.5, 0.5, 0.5] ])
        self.dofAttachments = [("C",0,None)]
        self.Create()

class Hexa_P1_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Hexa_P1_Lagrange,self).__init__()
        self.geoSupport = EN.GeoHex


        xi = self.xi
        eta = self.eta
        phi = self.phi

        self.symN = Matrix([(1-xi)*(1-eta)*(1-phi),
                            ( +xi)*(1-eta)*(1-phi),
                            ( +xi)*( +eta)*(1-phi),
                            (1-xi)*( +eta)*(1-phi),
                            (1-xi)*(1-eta)*( +phi),
                            ( +xi)*(1-eta)*( +phi),
                            ( +xi)*( +eta)*( +phi),
                            (1-xi)*( +eta)*( +phi)])

        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ("P",4,None),
                               ("P",5,None),
                               ("P",6,None),
                               ("P",7,None),
                               ]

        self.posN = np.array([[ 0, 0, 0],
                              [ 1, 0, 0],
                              [ 1, 1, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 1],
                              [ 1, 0, 1],
                              [ 1, 1, 1],
                              [ 0, 1, 1]])
        self.Create()

class Hexa_P2_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Hexa_P2_Lagrange,self).__init__()
        self.geoSupport = EN.GeoHex


        x1 = 1-self.xi
        x2 = self.xi

        e1 = 1-self.eta
        e2 = self.eta

        p1 = 1-self.phi
        p2 = self.phi

        xA = x1*(2*x1-1)
        xB = 4*x1*x2
        xC = x2*(2*x2-1)

        eA = e1*(2*e1-1)
        eB = 4*e1*e2
        eC = e2*(2*e2-1)

        pA = p1*(2*p1-1)
        pB = 4*p1*p2
        pC = p2*(2*p2-1)

        self.symN = Matrix([xA*eA*pA, # linear part
                            xC*eA*pA,
                            xC*eC*pA,
                            xA*eC*pA,
                            xA*eA*pC,
                            xC*eA*pC,
                            xC*eC*pC,
                            xA*eC*pC,

                            xB*eA*pA, #edges of base
                            xC*eB*pA,
                            xB*eC*pA,
                            xA*eB*pA,

                            xB*eA*pC, #edges of  top
                            xC*eB*pC,
                            xB*eC*pC,
                            xA*eB*pC,

                            xA*eA*pB, # vertical edges
                            xC*eA*pB,
                            xC*eC*pB,
                            xA*eC*pB,

                            xA*eB*pB,      # central faces
                            xC*eB*pB,
                            xB*eA*pB,
                            xB*eC*pB,
                            xB*eB*pA,
                            xB*eB*pC,

                            xB*eB*pB # central element

                            ])

        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ("P",4,None),
                               ("P",5,None),
                               ("P",6,None),
                               ("P",7,None),
                               ("F2",0,None),
                               ("F2",1,None),
                               ("F2",2,None),
                               ("F2",3,None),
                               ("F2",4,None),
                               ("F2",5,None),
                               ("F2",6,None),
                               ("F2",7,None),
                               ("F2",8,None),
                               ("F2",9,None),
                               ("F2",10,None),
                               ("F2",11,None),
                               ("F",0,None), #20
                               ("F",1,None), #21
                               ("F",2,None), #22
                               ("F",3,None), #23
                               ("F",4,None), #24
                               ("F",5,None), #25
                               ("C",0,None),

                               ]

        self.posN = np.array([[ 0, 0, 0],
                              [ 1, 0, 0],
                              [ 1, 1, 0],
                              [ 0, 1, 0],
                              [ 0, 0, 1],
                              [ 1, 0, 1],
                              [ 1, 1, 1],
                              [ 0., 1, 1],

                            [0.5,0.0,0.0], #edges of base
                            [1.0,0.5,0.0],
                            [0.5,1.0,0.0],
                            [0.0,0.5,0.0],

                            [0.5,0.0,1.0], #edges of  top
                            [1.0,0.5,1.0],
                            [0.5,1.0,1.0],
                            [0.0,0.5,1.0],

                            [0.0,0.0,0.5], # vertical edges
                            [1.0,0.0,0.5],
                            [1.0,1.0,0.5],
                            [0.0,1.0,0.5],

                            [0.0,0.5,0.5],      # central faces
                            [1.0,0.5,0.5],
                            [0.5,0.0,0.5],
                            [0.5,1.0,0.5],
                            [0.5,0.5,0.0],
                            [0.5,0.5,1.0],

                            [0.5,0.5,0.5] # central element

                              ])
        self.Create()

def CheckIntegrity(GUI=False):
    p0G =Hexa_P0_Global()
    p0L =Hexa_P0_Lagrange()
    p1L =Hexa_P1_Lagrange()
    p2L =Hexa_P2_Lagrange()
    print(p2L)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
