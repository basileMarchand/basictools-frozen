# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

import numpy as np
from sympy.matrices import Matrix
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase


class BarSpaceBase(SymSpaceBase):
    def __init__(self):
        super(BarSpaceBase,self).__init__()
        self.dimensionality = 1
        self.geoSupport = EN.GeoBar

class Bar_P0_Global(BarSpaceBase):
    def __init__(self):
        super(Bar_P0_Global,self).__init__()
        self.symN = Matrix([1])
        self.posN = np.array([[None]])
        self.dofAttachments = [("G",None,None)]
        self.Create()

class Bar_P0_Lagrange(BarSpaceBase):
    def __init__(self):
        super(Bar_P0_Lagrange,self).__init__()
        self.symN = Matrix([1])

        self.posN = np.array([[0.5]])
        self.dofAttachments = [("C",0,None) ]
        self.Create()

class Bar_P0_Legendre(Bar_P0_Lagrange):
    def __init__(self):
        super(Bar_P0_Legendre,self).__init__()

class Bar_P1_Lagrange(BarSpaceBase):
    def __init__(self):
        super(Bar_P1_Lagrange,self).__init__()
        xi = self.xi
        self.symN = Matrix([(1-xi), xi])
        self.posN = np.array([[0],
                              [1]])

        self.dofAttachments = [("P",0,None),
                               ("P",1,None)]
        self.Create()

# work in progress
class Bar_P2_Lagrange(BarSpaceBase):
    def __init__(self):
        super(Bar_P2_Lagrange,self).__init__()
        xi = self.xi
        L1 = 1-xi
        L2 = xi
        self.symN = Matrix([L1*(2*L1-1),L2*(2*L2-1),4*L1*L2 ])
        self.posN = np.array([[0],
                              [1],
                              [0.5],])

        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("C",0,None),]
        self.Create()


def CheckIntegrity(GUI=False):
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
