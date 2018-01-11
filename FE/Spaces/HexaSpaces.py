# -*- coding: utf-8 -*-
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase
import BasicTools.FE.ElementNames as EN
from sympy.matrices import Matrix
import numpy as np

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
        self.Create()

def CheckIntegrity():
    return "ok"
