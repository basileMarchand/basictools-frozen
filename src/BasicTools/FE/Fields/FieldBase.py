# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

unaryOps = {"__neg__":np.negative,
            "__abs__":np.abs}

binaryOps = {"__add__":np.add,
             "__mul__":np.multiply,
             "__pow__":np.power,
             "__sub__":np.subtract,
             "__rmul__":np.multiply,
             "__truediv__":np.divide,

             "__gt__":np.ndarray.__gt__,
             "__ge__":np.ndarray.__ge__,
             "__lt__":np.ndarray.__lt__,
             "__le__":np.ndarray.__le__,

             }


class FieldBase(BaseOutputObject):
    def __init__(self,name=None,mesh=None):
        super(FieldBase,self).__init__()
        if name is None:
            self.name = ""
        else:
            self.name = name
        self.mesh = mesh

    def GetName(self):
        return self.name

    def SetName(self,name):
        self.name = name

    def ConvertDataForNativeTreatment(self):
        self.data =  self.unaryOp(np.ascontiguousarray).data

    def __neg__(self):
        return self.unaryOp(np.negative)
    def __abs__(self):
        return self.unaryOp(np.abs)

    def __add__(self,other):
        return self.binaryOp(other,binaryOps["__add__"])
    def __radd__(self,other):
        return self.binaryOp(other,binaryOps["__add__"])
    def __mul__(self,other):
        return self.binaryOp(other,binaryOps["__mul__"])
    def __rmul__(self,other):
        return self.binaryOp(other,binaryOps["__mul__"])
    def __pow__(self,other):
        return self.binaryOp(other,binaryOps["__pow__"])
    def __gt__(self,other):
        return self.binaryOp(other,binaryOps["__gt__"])
    def __ge__(self,other):
        return self.binaryOp(other,binaryOps["__ge__"])
    def __lt__(self,other):
        return self.binaryOp(other,binaryOps["__lt__"])
    def __le__(self,other):
        return self.binaryOp(other,binaryOps["__le__"])

    def __sub__(self,other):
        return self.binaryOp(other,binaryOps["__sub__"])

    def __rsub__(self,other):
        res = -self
        res += other
        return res

    def __truediv__(self,other):
        return self.binaryOp(other,binaryOps["__truediv__"])

    def __getattr__(self,name):
        op = getattr(np,name,None)
        if op is None:
            raise(AttributeError(str(type(self)) + " does not have the '"+str(name)+"' attribute."))
            return
        def newfunc():
           res = self.unaryOp(op)
           return res
        return newfunc

def CheckIntegrity(GUI=False):
    obj = FieldBase("temp")
    print(obj)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
