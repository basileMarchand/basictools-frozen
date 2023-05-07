# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import List, Dict
import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class SymExprBase(BaseOutputObject):
    """Store read from string and store a symbolic expression.
    the expression read using sympy. the first and the second derivative
    are automatically computed.
    """
    def __init__(self, string=None, symbols=None):
        """create and set a symbolic expression. A list of symbols (list[str]) can be provided
         to determine the free variable of the expression. if no symbols are provided the symbol
         't' (for time) with value 0.0 is defined by default. """
        super(SymExprBase,self).__init__()

        self._expression = ""
        """The string representation of the expression. """
        self.constants = {}
        """Constants used to evaluate the expression. (free variables). """

        if symbols is None:
            self.SetConstant("t",0.0)
        else:
            for s in symbols:
                self.SetConstant(s,0.0)

        if string is not None:
            self.SetExpression(string)

    def SetExpression(self, string:str, _symbols: List[str]=None):
        """Set the expression to be used. A list of symbols (list[str]) can be provided
         to determine the free variable of the expression. if no symbols are provided the symbol
         't' (for time) with value 0.0 is defined by default.

        Parameters
        ----------
        string : str
            the string representation of the expression to be parsed by sympy
        _symbols : list[str], optional
            List of symbol to be used to parse the expression, by default None
        """
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

    def SetConstant(self, name:str, value:np.number):
        """Add/Set the value of the free variable of the expression

        Parameters
        ----------
        name : str
            name of the variable
        value : np.number
            value
        """
        self.constants[name]= value

    def GetValue(self, pos=None)-> np.number:
        """Return the evaluated expression

        Parameters
        ----------
        pos : _type_, optional
            Not Used, by default None

        Returns
        -------
        np.number
            the evaluated expression
        """

        return self.func(**self.constants)

    def GetValueDerivative(self, coor:str, pos=None) -> np.number:
        """Return the first derivative of the expression with respect to coor

        Parameters
        ----------
        coor : str
            the name of the variable to be used by the derivative
        pos : _type_, optional
            Not Used, by default None

        Returns
        -------
        np.number
            the evaluation of the derivative of the expression with respect to coor.
        """
        return self.dfuncd[coor](**self.constants)

    def GetValueSecondDerivative(self, coor1: str, coor2: str, pos=None) -> np.number:
        """Return the second derivative of the expression with respect to coor1 and coor2

        Parameters
        ----------
        coor1 : str
            the name of the variable to be used by the derivative
        coor2 : str
            the name of the variable to be used by the derivative
        pos : _type_, optional
            Not Used, by default None

        Returns
        -------
        np.number
            the evaluation of the derivative of the expression with respect to coor1 and coord2.
            d/dcoord1 * d/dcoord2 * expr
        """
        return self.d2funcd2[(coor1,coor2)](**self.constants)

    def __call__(self,pos=None) -> np.number:
        """Wrapper for the GetValue.
        Return the evaluated expression

        Parameters
        ----------
        pos : _type_, optional
            Not used, by default None

        Returns
        -------
        np.number
            the evaluated expression
        """
        return self.GetValue(pos)

class SymExprWithPos(SymExprBase):
    """Store read from string and store a symbolic expression depending implicitly on (x,y,z).
    the expression read using sympy. the first and the second derivative
    are automatically computed.
    """
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

def CreateSymExprWithPos(ops:Dict)->SymExprWithPos:
    """Simple wrapper to create a SymExprWithPos from a dict.
    ["val"] is used to extract the expression to be used for the contruction of
    the SymExprWithPos

    Parameters
    ----------
    ops : Dict
        _description_

    Returns
    -------
    SymExprWithPos
        the Symbolic expression dependent of the position
    """
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
