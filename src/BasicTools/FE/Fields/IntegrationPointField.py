# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       


import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

from BasicTools.Containers.UnstructuredMesh import AllElements

from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac as IntegrationRulesAlmanac
from BasicTools.Containers import ElementNames as EN


class IntegrationPointField(BaseOutputObject):
    def __init__(self,name=None,mesh=None,rules=None,data=None):
        super(IntegrationPointField,self).__init__()
        self.name = name
        
        if data is None:
            self.data = {}
        else:
            self.data = data
            
        #self.data["tetra"] = numpy( nuberOfElements, numberOfIntegrationPoints)
        self.mesh = mesh
        self.rules = rules # dict like  [EN.GeoTri] = tuple(points, weights)
        #self.filter = ElementFilter(dimensionality=3,tag="Toto")
        
    def SetRule(self,ruleName):
        self.rules = IntegrationRulesAlmanac[ruleName]
    
    def GetRuleFor(self,elemtype):
        return self.rules[EN.geoSupport[elemtype]]
        
    def GetFieldFor(self,elemtype):
            return self.data[elemtype]
    
    def Allocate(self):
        for name,data in self.mesh.elements.items():
            nbItegPoints = len(self.GetRuleFor(name)[1])
            nbElements = data.GetNumberOfElements()
            self.data[name] = np.zeros((nbElements,nbItegPoints), dtype=np.float)

    def __binaryOp__(self,other,op):
        self.__CheckCompatiblility(self,other)
        
        res = IntegrationPointField(name = "" )
        res.data = { key:op(A.data,B.data) for A,B in zip(self.data,other.data)}
        res.mesh = self.mesh
        res.rules = self.rules
        return res

    def __add__(self,other):
        return self.__binaryOp(other,np.add)
        
    def __CheckCompatiblility(self,A,B):
        if id(A.mesh) != id(B.mesh):
            raise (Exception("The support of the fields are not the same"))
        if id(A.rules) != id(B.rules):
            raise (Exception("The rules of the fields are not the same"))
            
    def __str__(self):
        res =  "IntegrationPointField\n"
        # ("+str(self.ncoomp)+")
        res += " name : "+ str(self.name) + "\n"
        res += str({ name:data.shape for name,data in self.data.items()} )  + "\n"
        return res


def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    mesh = CreateCube([10.,10.,10.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])

    eps = IntegrationPointField("Epsilon_0",mesh=mesh,rule=LagrangeP1)
    eps.Allocate()
    
    print(eps)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
