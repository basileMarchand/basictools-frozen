# -*- coding: utf-8 -*-

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class FElement(BaseOutputObject):
    def __init__(self):
        super(FElement,self).__init__()
        self.nnodes = -1
        self.name = None
        self.dimensionality = -1

    def GetDetJack(self,qcoor):
        raise# pragma: no cover

    def GetDimensionality(self):
        return self.dimensionality; #pragma: no cover

def CheckIntegrity():
    print(FElement().GetDimensionality())

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
