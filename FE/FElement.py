# -*- coding: utf-8 -*-

from OTTools.Helpers.BaseOutputObject import BaseOutputObject

class FElement(BaseOutputObject):
    def __init__(self):
        super(FElement,self).__init__()
        self.nnodes = -1
        self.dim = -1
        self.name = None
        self.dimensionality = -1

    def GetDetJack(self,qcoor):
        raise# pragma: no cover

    def GetDimensionality(self):
        return self.dimensionality; #pragma: no cover

def CheckIntegrity():
    print(FElement().GetDimensionality())

    return "OK"
