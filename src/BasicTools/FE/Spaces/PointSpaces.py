#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import numpy as np
from sympy.matrices import Matrix
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase


class PointSpaceBase(SymSpaceBase):
    def __init__(self):
        super(PointSpaceBase,self).__init__()
        self.dimensionality = 1
        self.geoSupport = EN.GeoBar

class Point_P0_Global(PointSpaceBase):
    def __init__(self):
        super(Point_P0_Global,self).__init__()
        self.symN = Matrix([1])
        self.posN = np.array([[None]])
        self.dofAttachments = [("G",None,None)]
        self.Create()

class Point_P0_Lagrange(PointSpaceBase):
    def __init__(self):
        super(Point_P0_Lagrange,self).__init__()
        self.symN = Matrix([1])

        self.posN = np.array([[0]])
        self.dofAttachments = [("C",0,None) ]
        self.Create()

def CheckIntegrity(GUI=False):
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
