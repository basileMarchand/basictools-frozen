# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class SymExprBase(BaseOutputObject):
    def __init__(self, string=None, symbols=None):
        super(SymExprBase,self).__init__()
        self._expression = ""
        self.constants = {}
        if symbols is None:
            self.SetConstant("t",0.0)
        else:
            for s in symbols:
                self.SetConstant(s,0.0)

        if string is not None:
            self.SetExpression(string)

    def SetExpression(self,string,_symbols=None):
        from sympy import Symbol
        from sympy import symbols
        from sympy import lambdify
        from sympy.parsing.sympy_parser import parse_expr

        if _symbols is None:
            self.stringSymbols = list(self.constants.keys())
            _symbols = self.stringSymbols

        self._expression = string
        _expr = parse_expr(self._expression)
        self.func = lambdify(symbols(self.stringSymbols),_expr)
        self.dfuncd = dict()
        self.d2funcd2 = dict()
        for s in  self.stringSymbols:
            self.dfuncd[s] = lambdify(symbols(self.stringSymbols),_expr.diff(Symbol(s)))
            for s2 in  self.stringSymbols:
                self.d2funcd2[(s,s2)] = lambdify(symbols(self.stringSymbols),_expr.diff(Symbol(s)).diff(Symbol(s)))

    def SetConstant(self,name,value):
        self.constants[name]= value

    def GetValue(self,pos=None):
        return self.func(**self.constants)

    def GetValueDerivative(self,coor,pos=None):
        return self.dfuncd[coor](**self.constants)

    def GetValueSecondDerivative(self,coor1,coor2,pos=None):
        return self.d2funcd2[(coor1,coor2)](**self.constants)

    def __call__(self,pos=None):
        return self.GetValue(pos)

class SymExprWithPos(SymExprBase):
    def __init__(self, string=None, symbols=None):
        super(SymExprWithPos,self).__init__(string=string, symbols=symbols)

    def SetExpression(self,string):
        self.stringSymbols = list(self.constants.keys())
        self.stringSymbols.extend("xyz")
        super().SetExpression(string,self.stringSymbols)

    def GetValue(self,pos):
        res = self.func(x=pos[:,0],y=pos[:,1],z=pos[:,2], **self.constants)
        if res.size == pos.shape[0]:
            return res
        else:
            return np.full((pos.shape[0],),fill_value=res)

    def GetValueDerivative(self,coor,pos):
        res =self.dfuncd[coor](x=pos[:,0],y=pos[:,1],z=pos[:,2], **self.constants)
        if res.size == pos.shape[0]:
            return res
        else:
            return np.full((pos.shape[0],),fill_value=res)

    def GetValueSecondDerivative(self,coor1,coor2,pos):
        res = np.asarray(self.d2funcd2[(coor1,coor2)](x=pos[:,0],y=pos[:,1],z=pos[:,2], **self.constants))
        if res.size == pos.shape[0]:
            return res
        else:
            return np.full((pos.shape[0],),fill_value=res)

    def __str__(self):
        res = f"SymExprWithPos('{self._expression}') "
        return res

    def __repr__(self):
        return self.__str__()

def CreateSymExprWithPos(ops):

    sym = SymExprWithPos()
    sym.SetExpression(ops["val"])
    return sym


def CheckIntegrity(GUI=False):
    #minthreshold="0.00000"
    string = """<Pressure  eTag="ET2" val="sin(3*t)+x**2" />"""

    import xml.etree.ElementTree as ET
    root = ET.fromstring(string)
    data = root.attrib

    data.pop("id",None)

    obj = CreateSymExprWithPos(data)


    obj.SetConstant("t",3.14159/6.)
    print(obj)
    print("data : ")
    data = np.array([[100.0,0.1,0.2 ],[0,0.1,0.2 ] ])
    print(data)
    print("f = ", obj._expression)
    print(obj.GetValue(data))
    print("dfdx :")
    print(obj.GetValueDerivative("x",data))
    print("dfdt :")
    c = obj.GetValueDerivative("t",data)
    print(c)
    print("d2fdx2 :")
    print(obj.GetValueSecondDerivative("x","x",data))

    import BasicTools.Containers.ElementNames as EN


    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube(dimensions=[10,11,12],spacing=[1.,1.,1.],ofTetras=False)

    if np.any(mesh.nodes[:,0]**2+np.sin(3*3.14159/6.) - obj.GetValue(mesh.nodes)):
        raise (ValueError("vectors does not match"))

    for name,data in mesh.elements.items():
        if EN.dimension[name] == 3:
             data.tags.CreateTag("Inside3D",False).SetIds(np.arange(data.GetNumberOfElements()))
             data.tags.CreateTag("Outside3D",False)
        if EN.dimension[name] == 2:
             data.tags.CreateTag("InterSurf",False).SetIds(np.arange(data.GetNumberOfElements()))

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
