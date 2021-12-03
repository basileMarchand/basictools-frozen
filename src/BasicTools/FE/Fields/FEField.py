# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

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

#    def GetValueAtIP(self,elemtype,el,ip):
#        sp = self.space[elemtype]
#        num = self.numbering[elemtype][el,:]
#        return sp.Eval_FieldI(ip,self.data[num],None,None,der=-1)

    def __str__(self):
        TFormat.II()
        res =  TFormat.InBlue("NodalField")+"\n"
        if self.name is not None:
          res += TFormat.GetIndent()
          res += TFormat.InGreen("Name : ") + self.name
        TFormat.DI()
        return res

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

    def GetCellRepresentation(self,fillvalue=0):
        """
        Function to push the data from the field into a vector homogeneous to
        the mesh (for visualisation for example). Entities with no dofs are filled
        with the fillvalues (default 0)
        """

        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)
        else:
            res = np.ones(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)*fillvalue

        if len(self.numbering["doftocellLeft"]) == 0:
            print("Warning : transfert vector is empty")
        res[self.numbering["doftocellLeft"]] = self.data[self.numbering["doftocellRight"]]

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
        elif type(other).__module__ == np.__name__ and np.ndim(other) != 0 :
            res = np.empty(other.shape,dtype=object)
            for res_data,other_data in np.nditer([res,other],flags=["refs_ok"],op_flags=["readwrite"]):
                res_data[...] = op(self,other_data)
            return res

        else:
            res.data = op(self.data,other)
        return res

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube([2.,3.,4.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])

    from BasicTools.FE.FETools import PrepareFEComputation
    spaces,numberings,offset, NGauss = PrepareFEComputation(mesh,numberOfComponents=1)

    sig11 = FEField(name = "temp",mesh=mesh,space=spaces,numbering=numberings[0])
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
