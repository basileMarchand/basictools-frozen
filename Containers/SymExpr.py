# -*- coding: utf-8 -*-
import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class SymExprBase(BaseOutputObject):
    def __init__(self,string=None):
        super(SymExprBase,self).__init__()
        self._expression = ""
        self.constants = {}
        self.SetConstants("t",0.0)
        if string is not None:
            self.SetExpression(string)

    def SetExpression(self,string):
        from sympy import symbols
        from sympy import lambdify
        from sympy.parsing.sympy_parser import parse_expr
        self.stringSymbols = list(self.constants.keys())
        self._expression = string
        self.func = lambdify(symbols(self.stringSymbols),parse_expr(self._expression))

    def SetConstants(self,name,value):
        self.constants[name]= value

    def GetValue(self,pos=None):
        return self.func(**self.constants)

    def __call__(self,pos=None):
        return self.GetValue(pos)

class SymExprWithPos(SymExprBase):
    def __init__(self):
        super(SymExprWithPos,self).__init__()
        self.eTag = None
        self.on = ""

    def SetExpression(self,string):
        from sympy import symbols
        from sympy import lambdify
        from sympy.parsing.sympy_parser import parse_expr
        self.stringSymbols = list(self.constants.keys())
        self.stringSymbols.extend("xyz")
        self._expression = string
        self.func = lambdify(symbols(self.stringSymbols),parse_expr(self._expression))

    def GetValue(self,pos):
        return self.func(x=pos[:,0],y=pos[:,1],z=pos[:,2], **self.constants)

def CreateSymExprWithPos(ops):

    sym = SymExprWithPos()
    sym.SetExpression(ops["val"])
    return sym


def CheckIntegrity(GUI=False):
    #minthreshold="0.00000"
    string = """<Pressure  eTag="ET2" val="sin(3*t)+x" />"""

    import xml.etree.ElementTree as ET
    root = ET.fromstring(string)

    from TopoTools.CLApp.MainApp import MainApp
    app = MainApp()

    data = app.XmlToDic(root)
    data.pop("id",None)

    obj = CreateSymExprWithPos(data)


    obj.SetConstants("t",3.14159/6.)
    print(obj)
    print(obj.GetValue(np.array([[100.0,0.1,0.2 ] ])))



    import BasicTools.Containers.ElementNames as EN
    import TopoTools.TopoZones as TZ
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    mesh = CreateCube(dimensions=[10,11,12],spacing=[1.,1.,1.],ofTetras=False)

    if np.any(mesh.nodes[:,0]+np.sin(3*3.14159/6.) - obj.GetValue(mesh.nodes)):
        raise (ValueError("vectors does not match"))

    for name,data in mesh.elements.items():
        if EN.dimension[name] == 3:
             data.tags.CreateTag(TZ.Inside3D,False).SetIds(np.arange(data.GetNumberOfElements()))
             data.tags.CreateTag(TZ.Outside3D,False)
        if EN.dimension[name] == 2:
             data.tags.CreateTag(TZ.InterSurf,False).SetIds(np.arange(data.GetNumberOfElements()))


    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
