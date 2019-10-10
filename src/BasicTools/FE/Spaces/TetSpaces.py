#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-
import numpy as np
from sympy.matrices import Matrix

import BasicTools.Containers.ElementNames as EN

from BasicTools.FE.Spaces.SymSpace import SymSpaceBase


class TetSpaceBase(SymSpaceBase):
    def __init__(self):
        super(TetSpaceBase,self).__init__()
        self.dimensionality = 3
        self.geoSupport = EN.GeoTet #tri


class Tet_P0_Global(TetSpaceBase):
    def __init__(self):
        super(Tet_P0_Global,self).__init__()
        self.symN = Matrix([1])
        self.posN = np.array([[None,None,None]])
        self.dofAttachments = [("G",None,None)]
        self.Create()


class Tet_P0_Lagrange(TetSpaceBase):
    def __init__(self):
        super(Tet_P0_Lagrange,self).__init__()
        self.symN = Matrix([1])
        self.posN = np.array([[0.25,0.25,0.25]])
        self.dofAttachments = [("C",0,None)]
        self.Create()

class Tet_P1_Lagrange(TetSpaceBase):
    def __init__(self):
        super(Tet_P1_Lagrange,self).__init__()
        xi = self.xi
        eta = self.eta
        phi = self.phi
        self.symN = Matrix([(1-xi-eta-phi), xi,eta, phi])
        self.posN = np.array([[0,0,0],
                              [1,0,0],
                              [0,1,0],
                              [0,0,1]])

        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None) ]
        self.Create()

class Tet_P2_Lagrange(TetSpaceBase):
    def __init__(self):
        super(Tet_P2_Lagrange,self).__init__()
        xi = self.xi
        eta = self.eta
        phi = self.phi

        T = (1-xi-eta-phi)

        self.symN = Matrix([eta*(2*eta-1),
                            T*(2*T-1),
                            xi*(2*xi-1),
                            phi*(2*phi-1),
                            4*eta*T,
                            4*T*xi,
                            4*xi*eta,
                            4*eta*phi,
                            4*T*phi,
                            4*xi*phi])

        self.posN = np.array([[0,1,0],
                              [0,0,0],
                              [1,0,0],
                              [0,0,1],
                              [0,0.5,0],
                              [0.5,0,0],
                              [0.5,0.5,0],
                              [0,0.5,0.5],
                              [0,0,0.5],
                              [0.5,0,0.5]])

        self.dofAttachments = [("P",2,None),
                               ("P",0,None),
                               ("P",1,None),
                               ("P",3,None),
                               ("F2",2,None),
                               ("F2",0,None),
                               ("F2",1,None),
                               ("F2",5,None),
                               ("F2",3,None),
                               ("F2",4,None)
                               ]


        self.Create()

def CheckIntegrity(GUI=False):
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
