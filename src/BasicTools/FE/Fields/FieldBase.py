# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

unaryOps = {"__neg__":np.negative}

binaryOps = {"__add__":np.add,
             "__mul__":np.multiply,
             "__pow__":np.power,
             "__sub__":np.subtract,
             "__rmul__":np.multiply,
             "__truediv__":np.divide
             }


class FieldBase(BaseOutputObject):
    def __init__(self,name=None,mesh=None):
        super(FieldBase,self).__init__()
        self.name = name
        self.mesh = mesh

    def __neg__(self):
        return self.unaryOp(np.negative)

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

    def __sub__(self,other):
        return self.binaryOp(other,binaryOps["__sub__"])

    def __truediv__(self,other):
        return self.binaryOp(other,binaryOps["__truediv__"])

    def __getattr__(self,name):
        op = getattr(np,name)
        if op is None:
            return None
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
