# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.FE.IntegrationsRules import GetRule
from BasicTools.Containers import ElementNames as EN
from BasicTools.FE.Fields.FieldBase import FieldBase
from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.Containers.Filters import ElementFilter


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
        self.rule = GetRule(ruleName=ruleName,rule=rule)

    def GetRuleFor(self,elemtype):
        return self.rule[elemtype]

    def GetDataFor(self,elemtype):
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

    def GetCellRepresentation(self,fillvalue=0,method='mean'):
        """
        Function to push the data from the field into a vector homogeneous to
        the mesh (for visualisation for example).
        """
        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfElements(),dtype=float)
        else:
            res = np.ones(self.mesh.GetNumberOfElements(),dtype=float)*fillvalue

        cpt =0
        for name,data in self.mesh.elements.items():
            nbelems = data.GetNumberOfElements()

            if name not in self.data:
                cpt += nbelems
                continue

            if method == 'mean':
                data = np.mean(self.data[name],axis=1)
            elif method == 'max':
                data = np.max(self.data[name],axis=1)
            elif method == 'min':
                data = np.min(self.data[name],axis=1)
            elif method == 'maxdiff' or method == "maxdifffraction":
                cols = self.data[name].shape[1]
                op = np.zeros( (cols,(cols*(cols-1))//2) )
                icpt = 0
                for i in range(0,cols-1):
                    for j in range(i+1,cols):
                        op[i,icpt] = 1
                        op[j,icpt] = -1
                        icpt += 1
                data = np.max(abs(self.data[name].dot(op)),axis=1)
                if method == "maxdifffraction":
                    data /= np.mean(self.data[name],axis=1)
            else:
                col = min(int(method),self.data[name].shape[1])
                data = self.data[name][:,col]

            res[cpt:cpt+nbelems] = data
            cpt += nbelems

        return res

    def GetPointRepresentation(self,fillvalue=0, method="mean",dim=None):
        cellData = self.GetCellRepresentation(fillvalue=fillvalue,method=method)
        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfNodes(),dtype=float)
        else:
            res = np.ones(self.mesh.GetNumberOfNodes(),dtype=float)*fillvalue

        pntcpt = np.zeros(self.mesh.GetNumberOfNodes(),dtype=int)
        cpt = 0
        for name,data in self.mesh.elements.items():
            if dim is not None:
                if EN.dimension[name] != dim:
                    cpt += data.GetNumberOfElements()
                    continue
            for i in range(data.GetNumberOfNodesPerElement()):
                res[data.connectivity[:,i]] += cellData[cpt:cpt+data.GetNumberOfElements()]
                pntcpt[data.connectivity[:,i]] += 1
            cpt += data.GetNumberOfElements()
        pntcpt[pntcpt==0] = 1
        res /= pntcpt
        return res

    def Flatten(self,dim=None):
        elements3D = ElementFilter(self.mesh,dimensionality=dim)
        nbvalues = 0
        for elemType,data,ids in elements3D:
            nbvalues += np.prod(self.data[elemType].shape)
        res = np.empty(nbvalues,dtype=float)
        cpt = 0
        for elemType,data,ids in elements3D:
            lsize= np.prod(self.data[elemType].shape)
            res[cpt:cpt+lsize] = self.data[elemType].flatten()
            cpt += lsize
        return res

    def SetDataFromNumpy(self, indata,dim=None):
        elements3D = ElementFilter(self.mesh,dimensionality=dim)
        nbvalues = 0
        for elemType,data,ids in elements3D:
            nbvalues += np.prod(self.data[elemType].shape)

        if indata.size != nbvalues:
            raise(Exception("incompatible sizes"))

        cpt = 0
        for elemType,data,ids in elements3D:
            lsize= np.prod(self.data[elemType].shape)

            self.data[elemType][:,:] = indata[cpt:cpt+lsize].reshape(self.data[elemType].shape)
            cpt += lsize


    def GetRestrictedIPField(self,efmask):
        res = RestrictedIPField(name=self.name,mesh=self.mesh,rule=self.rule,efmask=efmask)
        res.AllocateFromIpField(self)
        return res

class RestrictedIPField(IPField):
    def __init__(self,name=None,mesh=None,rule=None,ruleName=None,data=None,efmask=None):
        super(RestrictedIPField,self).__init__(name=name, mesh=mesh,rule=rule,ruleName=ruleName,data=data)
        if efmask == None:
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

    def GetRestrictedIPField(self,efmask):
        from BasicTools.Containers.Filters import IntersectionElementFilter
        res = RestrictedIPField(name=self.name,mesh=self.mesh,rule=self.rule,efmask= IntersectionElementFilter(filters=[efmask,self.efmask]) )
        print("self.efmask")
        print(self.efmask)
        print("efmask")
        print(efmask)
        print("res.efmask")
        print(res.efmask)
        res.AllocateFromIpField(self)
        return res

    def AllocateFromIpField(self,ipField):
        self.name = ipField.name
        self.mesh = ipField.mesh
        self.rule = ipField.rule
        self.efmask.SetMesh(self.mesh)
        self.data = dict()

        if isinstance(ipField,RestrictedIPField):
            for name,data,ids in self.efmask :
                print('toto')
                idsII =  ipField.efmask.GetIdsToTreat(data)
                idds = np.empty(data.GetNumberOfElements(),dtype=int)
                idds[idsII] = np.arange(len(idsII))
                self.data[name] = ipField.data[name][idds[ids],:]
        else:
            for name,data,ids in self.efmask :
                self.data[name] = ipField.data[name][ids,:]

    def GetIpFieldRepr(self,fillvalue=0):
        res = IPField(name=self.name,mesh=self.mesh,rule=self.rule)
        res.Allocate(fillvalue)
        for name,data,ids in self.efmask :
            if name not in self.data:
                continue
            res.data[name][ids,:] = self.data[name]
        return res

    def CheckCompatiblility(self,B):
        super(RestrictedIPField,self).CheckCompatiblility(B)
        if isinstance(B,type(self)):
            if not self.efmask.IsEquivalent(B.efmask):
                raise (Exception("The efmask of the fields are not the same"))

    def unaryOp(self,op,out=None):
        res = type(self)(name = None, mesh=self.mesh,rule=self.rule, efmask=self.efmask)
        return super(RestrictedIPField,self).unaryOp(op,out=res)

    def binaryOp(self,other,op,out=None):
        res = type(self)(name = None, mesh=self.mesh,rule=self.rule, efmask=self.efmask)
        return super(RestrictedIPField,self).binaryOp(other,op,out=res)

    def GetCellRepresentation(self,fillvalue=0,method='mean'):
        """
        Function to push the data from the field into a vector homogeneous to
        the mesh (for visualisation for example).
        """
        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfElements(),dtype=float)
        else:
            res = np.ones(self.mesh.GetNumberOfElements(),dtype=float)*fillvalue

        cpt =0
        self.efmask.SetMesh(self.mesh)
        for name, eldata in self.mesh.elements.items():
            ids = self.efmask.GetIdsToTreat(eldata)
            nbelems = eldata.GetNumberOfElements()

            if name not in self.data:
                cpt += nbelems
                continue

            if method == 'mean':
                data = np.mean(self.data[name],axis=1)
            elif method == 'max':
                data = np.max(self.data[name],axis=1)
            elif method == 'min':
                data = np.min(self.data[name],axis=1)
            elif method == 'maxdiff' or method == "maxdifffraction":
                cols = self.data[name].shape[1]
                op = np.zeros( (cols,(cols*(cols-1))//2) )
                icpt = 0
                for i in range(0,cols-1):
                    for j in range(i+1,cols):
                        op[i,icpt] = 1
                        op[j,icpt] = -1
                        icpt += 1
                data = np.max(abs(self.data[name].dot(op)),axis=1)
                if method == "maxdifffraction":
                    data /= np.mean(self.data[name],axis=1)
            else:
                col = min(int(method),self.data[name].shape[1])
                data = self.data[name][:,col]
            res[cpt+np.array(ids,dtype=int)]=data
            cpt += nbelems

        return res

def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube([2.,3.,4.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])
    print(mesh)
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

    data = sig22.GetCellRepresentation()
    sig22.SetDataFromNumpy(sig22.Flatten())
    data2 = sig22.GetCellRepresentation()

    if np.linalg.norm(data-data2) > 0  :
        raise() # pragma: no cover

    data2 = sig22.GetPointRepresentation(0,"max")
    data2 = sig22.GetPointRepresentation(0,"min")
    data2 = sig22.GetPointRepresentation(0,"maxdiff")
    data2 = sig22.GetPointRepresentation(0,"maxdifffraction")
    data2 = sig22.GetPointRepresentation(0)
    data2 = sig22.GetPointRepresentation(-1,0)

    dummyField = IPField()
    dummyField.data[None] = np.arange(3)+1
    print("dummyField")
    print(dummyField*dummyField-dummyField**2/dummyField)

    restrictedIPField = sig22.GetRestrictedIPField(ElementFilter(tag="Skin"))
    restrictedIPField = -(2*restrictedIPField)

    print(restrictedIPField )
    data2 = restrictedIPField.GetPointRepresentation(0,"max")
    data2 = restrictedIPField.GetPointRepresentation(0,"min")
    data2 = restrictedIPField.GetPointRepresentation(0,"maxdiff")
    data2 = restrictedIPField.GetPointRepresentation(0,"maxdifffraction")
    data2 = restrictedIPField.GetPointRepresentation(0)
    data2 = restrictedIPField.GetPointRepresentation(-1,0)
    restrictedIPField.GetIpFieldRepr()
    #restrictedIPField2 = restrictedIPField.GetRestrictedIPField(ElementFilter(dimensionality=3))
    #restrictedIPField2.Allocate()
    restrictedIPField2 = restrictedIPField.GetRestrictedIPField(ElementFilter(tag="X0"))
    #restrictedIPField2.Allocate()
    print(restrictedIPField2.data)

    restrictedIPField2 = -(2*sig22.GetRestrictedIPField(ElementFilter(tags=["X0"])))
    print(restrictedIPField2.data)
    restrictedIPField2.Allocate(1)
    print(restrictedIPField2*np.array([0.1, 0.3]))
    import BasicTools.Containers.ElementNames as EN
    restrictedIPField2.GetDataFor(EN.Quadrangle_4)

    obj = RestrictedIPField(data={})
    obj.SetRule()

    #mustfail
    error = False
    try:
        restrictedIPField+restrictedIPField2
        error = True# pragma: no cover
    except:
        pass
    if error:# pragma: no cover
        raise

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
