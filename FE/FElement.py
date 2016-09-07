# -*- coding: utf-8 -*-

from OTTools.Helpers.BaseOutputObject import BaseOutputObject

class FElement(BaseOutputObject):
    def __init__(self):
        super(FElement,self).__init__()
        self.nnodes = 0
        self.dim = 0
        self.name = None
        
    def GetDetJack(self,qcoor):
        raise# pragma: no cover 
        
    def GetDimensionality(self):
        return self.dimensionality;
        
def CheckIntegrity():
    return "OK"
