# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Union
import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType
from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.FE.Fields.FieldBase import FieldBase


class FEField(FieldBase):
    def __init__(self,name=None,mesh=None,space=None,numbering=None,data=None):
        super(FEField,self).__init__(name=name,mesh = mesh)
        self.data = data
        self.space = space
        self.numbering = numbering

    def Allocate(self,val=0):
        if val == 0:
            self.data = np.zeros(self.numbering["size"],dtype=PBasicFloatType)
        else:
            self.data = np.ones(self.numbering["size"],dtype=PBasicFloatType)*val

    def copy(self):
        return FEField(name=self.name,mesh=self.mesh,space=self.space,numbering=self.numbering, data=self.data.copy())

#    def GetValueAtIP(self,elemtype,el,ip):
#        sp = self.space[elemtype]
#        num = self.numbering[elemtype][el,:]
#        return sp.Eval_FieldI(ip,self.data[num],None,None,der=-1)



    def GetPointRepresentation(self,fillvalue=0):
        """
        Function to push the data from the field into a vector homogeneous to
        the mesh (for visualisation for example). Entities with no dofs are filled
        with the fillvalues (default 0)
        """

        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfNodes(),dtype=PBasicFloatType)
        else:
            res = np.ones(self.mesh.GetNumberOfNodes(),dtype=PBasicFloatType)*fillvalue

        if len(self.numbering["doftopointLeft"]) == 0:
            print("Warning : transfert vector is empty")

        res[self.numbering["doftopointLeft"]] = self.data[self.numbering["doftopointRight"]]

        return res

    def SetDataFromPointRepresentation(self,userdata, fillvalue=0.):
        if fillvalue==0.:
           self.data = np.zeros(self.numbering["size"])
        else:
           self.data = np.ones(self.numbering["size"])*fillvalue

        self.data[self.numbering["doftopointRight"]] = userdata[self.numbering["doftopointLeft"]]

    def GetCellRepresentation(self,fillvalue:PBasicFloatType=0, method:Union[str,int]='mean') -> np.ndarray:
        """
        Function to push the data from the field into a vector homogeneous to
        the mesh cell (for visualisation for example). Entities with no dofs are filled
        with the fillvalues (default 0)
        the method controls to transfer function
        """

        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)
        else:
            res = np.ones(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)*fillvalue

        cpt =0
        for name, data in self.mesh.elements.items():
            nbelems = data.GetNumberOfElements()

            numbering = self.numbering[name]
            if name is None:
                cpt += nbelems
                continue

            if method == 'mean':
                data = np.mean(self.data[numbering],axis=1)
            elif method == 'max':
                data = np.max(self.data[numbering],axis=1)
            elif method == 'min':
                data = np.min(self.data[numbering],axis=1)
            elif method == 'maxdiff' or method == "maxdifffraction":
                cols = self.data[numbering].shape[1]
                op = np.zeros( (cols,(cols*(cols-1))//2) )
                icpt = 0
                for i in range(0,cols-1):
                    for j in range(i+1,cols):
                        op[i,icpt] = 1
                        op[j,icpt] = -1
                        icpt += 1
                data = np.max(abs(self.data[numbering].dot(op)),axis=1)
                if method == "maxdifffraction":
                    data /= np.mean(self.data[numbering],axis=1)
            else:
                col = min(int(method),self.data[numbering].shape[1])
                data = self.data[numbering][:,col]

            res[cpt:cpt+nbelems] = data
            cpt += nbelems

        return res

    def CheckCompatiblility(self,B):
        if isinstance(B,type(self)):
            if id(self.mesh) != id(B.mesh):
                raise (Exception("The support of the fields are not the same"))
            if id(self.space) != id(B.space):
                raise (Exception("The space of the fields are not the same"))
            if id(self.numbering) != id(B.numbering):
                raise (Exception("The numbering of the fields are not the same"))

    def unaryOp(self,op):
        res = type(self)(name = None,mesh=self.mesh,space=self.space, numbering = self.numbering )
        res.data = op(self.data)
        return res

    def binaryOp(self,other,op):
        self.CheckCompatiblility(other)
        res = type(self)(name = None,mesh=self.mesh,space=self.space, numbering = self.numbering )
        if isinstance(other,type(self)):
            res.data = op(self.data,other.data)
            return res
        elif type(other).__module__ == np.__name__ and np.ndim(other) != 0 :
            res = np.empty(other.shape,dtype=object)
            for res_data,other_data in np.nditer([res,other],flags=["refs_ok"],op_flags=["readwrite"]):
                res_data[...] = op(self,other_data)
            return res
        elif np.isscalar(other):
            res.data = op(self.data,other)
            return res
        else:
            raise Exception(f"operator {op} not valid for types :{type(self)} and {type(other)} ")

    def GetTestField(self):
        from BasicTools.FE.SymWeakForm import  GetTestSufixChar
        tc = GetTestSufixChar()
        if len(self.name) == 0:
            raise Exception("FEField must have a name")
        elif self.name[-1] == tc:
            raise Exception("this FEField is already a test field")
        else:
            return FEField(name=self.name+tc,mesh=self.mesh,space=self.space,numbering=self.numbering)

    def __str__(self):
        return  f"<FEField object '{self.name}' ({id(self)})>"

    def __repr__(self):
        return str(self)

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube([2.,3.,4.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])

    from BasicTools.FE.FETools import PrepareFEComputation
    spaces,numberings,offset, NGauss = PrepareFEComputation(mesh,numberOfComponents=1)

    sig11 = FEField(name = "temp",mesh=mesh,space=spaces,numbering=numberings[0])
    print(sig11.GetTestField())
    sig11.Allocate()
    print(sig11)

    sig22 = sig11+0.707107
    sig12 = 2*(-sig22)*5

    A = sig11**2
    B = sig11*sig22
    C = sig22**2
    D = 1.5*sig12*2
    E = A-B+C+(D)**2
    vonMises = np.sqrt(E)

    print(vonMises.data)
    print(np.linalg.norm([sig22, sig22 ] ).data )
    print(A/C)

    dummyField = FEField()
    dummyField.data = np.arange(3)+1

    print("dummyField")
    print(dummyField*dummyField-dummyField**2/dummyField)

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
