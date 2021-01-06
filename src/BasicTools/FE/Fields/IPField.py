# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac as IntegrationRulesAlmanac
from BasicTools.Containers import ElementNames as EN
from BasicTools.FE.Fields.FieldBase import FieldBase
from BasicTools.Helpers.TextFormatHelper import TFormat


class IPField(FieldBase):
    def __init__(self,name=None,mesh=None,rule=None,ruleName=None,data=None):
        super(IPField,self).__init__(name=name, mesh=mesh)
        if data is None:
            self.data = {}
        else:
            self.data = data
        self.rule = None
        self.SetRule(ruleName=ruleName,rule=rule)

    def SetRule(self,ruleName=None,rule=None):
        if ruleName is not None and rule is not None:
            raise(Exception("must give ruleName or rule not both"))
        if ruleName is not None:
            self.rule = IntegrationRulesAlmanac[ruleName]
        else :
            self.rule = rule

    def GetRuleFor(self,elemtype):
        return self.rule[elemtype]

    def GetFieldFor(self,elemtype):
        return self.data[elemtype]

    def Allocate(self,val=0):
        self.data = dict()
        for name,data in self.mesh.elements.items():
            nbItegPoints = len(self.GetRuleFor(name)[1])
            nbElements = data.GetNumberOfElements()
            self.data[name] = np.zeros((nbElements,nbItegPoints), dtype=np.float)+val

    def CheckCompatiblility(self,B):
        if isinstance(B,type(self)):
            if id(self.mesh) != id(B.mesh):
                raise (Exception("The support of the fields are not the same"))
            if id(self.rule) != id(B.rule):
                raise (Exception("The rules of the fields are not the same"))

    def __str__(self):
        res =  TFormat.InBlue("IPField")+"\n"
        # ("+str(self.ncoomp)+")
        if self.name is not None:
            res += " name : "+ str(self.name) + "\n"
            res += str({ name:data.shape for name,data in self.data.items()} )  + "\n"
        return res

    def unaryOp(self, op, out=None):
        if out is None:
            res = type(self)(name = None, mesh=self.mesh,rule=self.rule  )
        else:
            res = out
        res = type(self)(name = None, mesh=self.mesh,rule=self.rule )
        res.data = { key:op(self.data[key]) for key in self.data.keys()}
        return res

    def binaryOp(self,other,op,out=None):
        self.CheckCompatiblility(other)
        if out is None:
            res = type(self)(name = None, mesh=self.mesh,rule=self.rule  )
        else:
            res = out

        if isinstance(other,type(self)):
            res.data = { key:op(self.data[key],other.data[key]) for key in set(self.data.keys()).union(other.data.keys())}
            return res
        elif type(other).__module__ == np.__name__ and np.ndim(other) != 0:
            res = np.empty(other.shape,dtype=object)
            for res_data,other_data in np.nditer([res,other],flags=["refs_ok"],op_flags=["readwrite"]):
                res_data[...] = op(self,other_data)
            return res

        res.data = { key:op(self.data[key],other) for key in self.data.keys()}

        return res



class RestrictedIPField(IPField):
    def __init__(self,name=None,mesh=None,rule=None,ruleName=None,data=None,efmask=None):
        super(RestrictedIPField,self).__init__(name=name, mesh=mesh,rule=rule,ruleName=ruleName,data=data)
        if efmask == None:
            from BasicTools.Containers.Filters import ElementFilter
            self.efmask = ElementFilter()
        else:
            self.efmask = efmask

    def Allocate(self,val=0):
        self.efmask.SetMesh(self.mesh)
        self.data = dict()
        for name,data,ids in self.efmask :
            nbItegPoints = len(self.GetRuleFor(name)[1])
            nbElements = len(ids)
            self.data[name] = np.zeros((nbElements,nbItegPoints), dtype=np.float)+val

    def AllocateFromIpField(self,ipField):
        self.name = ipField.name
        self.mesh = ipField.mesh
        self.rule = ipField.rule
        self.efmask.SetMesh(self.mesh)
        self.data = dict()
        for name,data,ids in self.efmask :
            self.data[name] = ipField.data[name][ids,:]

    def GetIpFieldRepr(self,fillvalue=0):
        res = IPField(name=self.name,mesh=self.mesh,rule=self.rule)
        res.Allocate(fillvalue)
        for name,data,ids in self.efmask :
            res.data[name][ids,:] = self.data[name]
        return res

    def CheckCompatiblility(self,B):
        super(RestrictedIPField,self).CheckCompatiblility(B)
        if id(self.efmask) != id(B.efmask):
           raise (Exception("The efmask of the fields are not the same"))

    def unaryOp(self,op,out=None):
        res = type(self)(name = None, mesh=self.mesh,rule=self.rule, efmask=self.efmask)
        return super(RestrictedIPField,self).unaryOp(op,out=res)

    def binaryOp(self,other,op,out=None):
        res = type(self)(name = None, mesh=self.mesh,rule=self.rule, efmask=self.efmask)
        return super(RestrictedIPField,self).binaryOp(other,op,out=res)


def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube([2.,3.,4.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])

    sig11 = IPField("Sig_11",mesh=mesh,rule=LagrangeP1)
    sig11.Allocate()
    print(sig11)

    sig22 = sig11+0.707107
    sig12 = 2*(-sig22)*5/sig22

    A = sig11**2
    B = sig11*sig22
    C = sig22**2
    D = 1.5*sig12*2
    E = A-B+C+(D)**2
    vonMises = np.sqrt(E)

    print(vonMises.data)
    print("545454")
    print(np.linalg.norm([sig22, sig22 ] ).data )



    dummyField = IPField()
    dummyField.data[None] = np.arange(3)+1
    print("dummyField")
    print(dummyField*dummyField-dummyField**2/dummyField)


    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
