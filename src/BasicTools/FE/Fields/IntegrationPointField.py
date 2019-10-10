#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-


import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

from BasicTools.Containers.UnstructuredMesh import AllElements

from BasicTools.FE.IntegrationsRules import Lagrange as Lagrange


class IntegrationPointField(BaseOutputObject):
    def __init__(self,name):
        super(IntegrationPointField,self).__init__()
        self.storage = {}
        self.name = name
        #self.ncoomp = 1
        pass
    def Allocate(self,mesh,integrationrule, tag = AllElements):
        for name,data in mesh.elements.items():
            nbItegPoints = len(Lagrange(name)[1])
            if tag == AllElements:
                nbElements = data.GetNumberOfElements()
            else :
                if tag in data.tags:
                    nbElements = len(data.tags[tag])
                else:
                    nbElements = 0
            self.storage[name] = np.zeros((nbElements,nbItegPoints), dtype=np.float)

    def __str__(self):
        res =  "IntegrationPointField\n"
        # ("+str(self.ncoomp)+")
        res += " name : "+ str(self.name) + "\n"
        res += str({ name:data.shape for name,data in self.storage.items()} )  + "\n"
        return res

    def GetValueAtIP(self,name,el,ip):
        return self.storage[name][el,ip]

    def SetValueAtIP(self,name,el,ip,val):
        self.storage[name][el,ip] = val

    def IncrementValueAtIP(self,name,el,ip,val):
        self.storage[name][el,ip] += val

    def GetStorage(self,name):
        return self.storage[name]

def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    mesh = CreateCube([10.,10.,10.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])

    eps = IntegrationPointField("Epsilon_0")
    eps.ncoomp = 6
    eps.Allocate(mesh,LagrangeP1)
    print(eps)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
