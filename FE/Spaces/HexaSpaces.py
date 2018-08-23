# -*- coding: utf-8 -*-
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase
import BasicTools.FE.ElementNames as EN
from sympy.matrices import Matrix
import numpy as np


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

def CheckIntegrity():
    return "ok"
