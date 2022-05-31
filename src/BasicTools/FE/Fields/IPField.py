# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Dict, Any, Callable, Union

import numpy as np

from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType, ArrayLike
from BasicTools.FE.IntegrationsRules import GetRule, IntegrationRulesType, ElementIntegrationRuleType
from BasicTools.Containers import ElementNames as EN
from BasicTools.FE.Fields.FieldBase import FieldBase
from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.Containers.Filters import ElementFilter, IntersectionElementFilter
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh

class IPField(FieldBase):
    """Class to hold information at the integration points
    """
    def __init__(self, name: str = None, mesh: UnstructuredMesh =None, rule: IntegrationRulesType=None, ruleName:str=None, data: Dict[str,ArrayLike] = None):
        """Constructor

        Parameters
        ----------
        name : str, optional
            The name of the field, by default None
        mesh : UnstructuredMesh, optional
            the mesh , by default None
        rule : IntegrationRulesType, optional
            the integration rule to use, by default None
        ruleName : str, optional
            the name of the integration rule to use, by default None
        data : Dict[str,ArrayLike], optional
            the data per element type, the ArrayLike must be of size (number of elements, number of integration points), by default None
        """
        super().__init__(name=name, mesh=mesh)
        if data is None:
            self.data = {}
        else:
            self.data = data

        self.rule = None
        self.SetRule(ruleName=ruleName,rule=rule)

    def SetRule(self, ruleName:str=None, rule:IntegrationRulesType=None) -> None:
        """Set the integration rule to be used by the field

        Please read the constructor documentation for more information about the arguments
        """
        self.rule = GetRule(ruleName=ruleName,rule=rule)

    def GetRuleFor(self, elementType: str) -> ElementIntegrationRuleType:
        """Helper function to get the integration rule for one type of elements

        Parameters
        ----------
        elementType : str
            the element name

        Returns
        -------
        ElementIntegrationRuleType
            The integration rule, this is a tuple with the positions and the weight of the integration rule
        """
        return self.rule[elementType]

    def GetDataFor(self, elementType:str) -> np.ndarray:
        """Get the numpy array storing the data for one type of element

        Parameters
        ----------
        elementType : str
            the element name

        Returns
        -------
        np.ndarray
            the data of size (number of elements, number of integration points)
        """
        return self.data[elementType]

    def Allocate(self, val:PBasicFloatType= 0.) :
        """Allocate the memory to store the data

        Parameters
        ----------
        val : PBasicFloatType, optional
            initial fill value, by default 0.
        """
        self.data = dict()
        for name,data in self.mesh.elements.items():
            nbIntegrationPoints = len(self.GetRuleFor(name)[1])
            nbElements = data.GetNumberOfElements()
            self.data[name] = np.zeros((nbElements,nbIntegrationPoints), dtype=PBasicFloatType)+val

    def CheckCompatibility(self, B:Any) -> None:
        """Check of two field are compatible to make arithmetic operations

        Parameters
        ----------
        B : Any
            the other IpField

        Raises
        ------
        Exception
            In the case the fields are not compatible
        """
        if isinstance(B,type(self)):
            if id(self.mesh) != id(B.mesh): # pragma: no cover
                raise Exception("The support of the fields are not the same")
            if id(self.rule) != id(B.rule): # pragma: no cover
                raise Exception("The rules of the fields are not the same")

    def unaryOp(self, op:Callable, out: IPField=None) -> IPField:
        """Internal function to apply a unary operation

        Parameters
        ----------
        op : Callable
            a function to apply to the data
        out : IPField, optional
            output IPFiled if non is provided a new one is created, by default None

        Returns
        -------
        IPField
            The output IPField
        """
        if out is None:
            res = type(self)(name = None, mesh=self.mesh,rule=self.rule  )
        else:
            res = out
        res.data = { key:op(self.data[key]) for key in self.data.keys()}
        return res

    def binaryOp(self, other: Any , op:Callable, out:IPField=None) -> IPField:
        """Internal function to apply a binary operator. A + B for example


        Parameters
        ----------
        other : Any
            a second IPField
        op : Callable
            the binary operator to apply
        out : IPField, optional
            output IPFiled if Non is provided a new one is created, by default None, by default None

        Returns
        -------
        IPField
            The output IPField

        Raises
        ------
        Exception
            if the field are not compatible
        Exception
            if the binary operator cant be apply correctly
        """
        self.CheckCompatibility(other)
        if out is None:
            res = type(self)(name = None, mesh=self.mesh,rule=self.rule  )
        else:
            res = out

        if isinstance(other,(type(self),RestrictedIPField) ):
            res.data = { key:op(self.data[key],other.data[key]) for key in set(self.data.keys()).union(other.data.keys())}
            return res
        elif type(other).__module__ == np.__name__ and np.ndim(other) != 0:
            res = np.empty(other.shape,dtype=object)
            for res_data,other_data in np.nditer([res,other],flags=["refs_ok"],op_flags=["readwrite"]):
                res_data[...] = op(self,other_data)
            return res
        elif np.isscalar(other):
            res.data = { key:op(self.data[key],other) for key in self.data.keys()}
            return res
        else:
            raise Exception(f"operator {op} not valid for types :{type(self)} and {type(other)} ") # pragma: no cover

    def GetCellRepresentation(self, fillValue:PBasicFloatType= 0., method:Union[str,int]='mean') -> np.ndarray:
        """Function to push the data from the field into a vector homogeneous to
        the mesh (for visualisation for example).

        Parameters
        ----------
        fillValue : PBasicFloatType, optional
            field value for cell without data, by default 0.
        method : Union[str,int], optional
            the method to pass the information from the integration points
            to the cells:
            'mean'   : compute the mean value in every cell
            'max'    : extract the max value for every cell
            'min'    : extract the min value for every cell
            'maxdiff': compute the maximal difference for every cell max(abs(vi-vj)) for i!=j in range(number of integration point)
            'maxdifffraction': same as before but divided by the mean value
            'int' or, int: a int select a specific integration point the value is clipped to [0,number of integration point[
            , by default 'mean'

        Returns
        -------
        np.ndarray
            np.ndarray of size (number of elements in the mesh, 1) with the extracted information
        """
        if fillValue==0.:
            res = np.zeros(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)
        else:
            res = np.ones(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)*fillValue

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

    def GetPointRepresentation(self, fillValue:PBasicFloatType= 0., elementFilter: ElementFilter = None) -> np.ndarray:
        """Generate a point representation of this field

        Parameters
        ----------
        fillValue : float, optional
            value to be used on nodes with no data available, by default 0.
        elementFilter : ElementFilter, optional
            elements to construct the point representation, by default None

        Returns
        -------
        np.ndarray
            a numpy vector of size mesh.GetNumberOfNodes with the computed values
        """

        if elementFilter is None:
            elementFilter = ElementFilter()

        elementFilter.mesh = self.mesh

        res = np.zeros(self.mesh.GetNumberOfNodes(),dtype=PBasicFloatType)
        pointCpt = np.zeros(self.mesh.GetNumberOfNodes(),dtype=PBasicIndexType)

        from BasicTools.FE.Fields.FieldTools import ElementWiseIpToFETransferOp
        interpolator = ElementWiseIpToFETransferOp(self.rule,LagrangeSpaceGeo)

        for name, data, ids in elementFilter:
            if name not in interpolator or name not in self.data:
                continue
            for i in range(data.GetNumberOfNodesPerElement()):
                res[data.connectivity[ids,i]] +=  interpolator[name][i,:].dot(self.data[name][ids,:].T)
                pointCpt[data.connectivity[ids,i]] += 1
        mask = pointCpt==0
        pointCpt[mask] = 1
        res[mask] = fillValue
        res /= pointCpt
        return res

    def Flatten(self, elementFilter: ElementFilter = None) -> np.ndarray:
        """Flatten all the data into a single ndarray

        see also: self.SetDataFromNumpy

        Parameters
        ----------
        elementFilter : ElementFilter, optional
            ElementFilter to operate, by default (None) all element are treated

        Returns
        -------
        np.ndarray
            the extracted data
        """
        if elementFilter is None:
            elementFilter = ElementFilter()
        elementFilter.mesh = self.mesh

        nbValues = 0
        for elemType,data,ids in elementFilter:
            nbValues += np.prod(self.data[elemType].shape)
        res = np.empty(nbValues,dtype=PBasicFloatType)
        cpt = 0
        for elemType,data,ids in elementFilter:
            localSize= np.prod(self.data[elemType].shape)
            res[cpt:cpt+localSize] = self.data[elemType].flatten()
            cpt += localSize
        return res

    def SetDataFromNumpy(self, inData:ArrayLike, elementFilter: ElementFilter = None) -> None:
        """Set the internal data from a one dimensional array

        see also: self.Flatten

        Parameters
        ----------
        inData : ArrayLike
            the incoming data
        elementFilter : ElementFilter, optional
            ElementFilter to operate, by default (None) all element are treated
        """
        if elementFilter is None:
            elementFilter = ElementFilter()

        elementFilter.mesh = self.mesh

        nbValues = 0
        for elemType,data,ids in elementFilter:
            nbValues += np.prod(self.data[elemType].shape)

        if inData.size != nbValues: # pragma: no cover
            raise(Exception("incompatible sizes"))

        cpt = 0
        for elemType, data, ids in elementFilter:
            localSize = np.prod(self.data[elemType].shape)

            self.data[elemType][:,:] = inData[cpt:cpt+localSize].reshape(self.data[elemType].shape)
            cpt += localSize

    def GetRestrictedIPField(self, elementFilter:ElementFilter) -> RestrictedIPField:
        """Create a RestrictedIPField only on element selected by the elementFilter

        Parameters
        ----------
        elementFilter : ElementFilter
             ElementFilter to work on

        Returns
        -------
        RestrictedIPField
            a new RestrictedIPField
        """
        res = RestrictedIPField(name=self.name,mesh=self.mesh,rule=self.rule,elementFilter=elementFilter)
        res.AllocateFromIpField(self)
        return res

    def __str__(self) -> str:
        """Nice string representations

        Returns
        -------
        str
            a string in the form of "<IPField object 'name' (id(self))>"
        """
        return  f"<IPField object '{self.name}' ({id(self)})>"

class RestrictedIPField(IPField):
    """Class to hold information at integration point on a subdomain of the mesh (ElementFilter)
    This class is experimental and for the moment not compatible with the integration module

    """
    def __init__(self, name: str=None, mesh:UnstructuredMesh=None, rule:IntegrationRulesType=None, ruleName:str=None, data:Dict[str,ArrayLike]=None, elementFilter:ElementFilter=None) -> None:
        """Constructor

        Parameters
        ----------
        name : str, optional
            The name of the field, by default None
        mesh : UnstructuredMesh, optional
            the mesh , by default None
        rule : IntegrationRulesType, optional
            the integration rule to use, by default None
        ruleName : str, optional
            the name of the integration rule to use, by default None
        data : Dict[str,ArrayLike], optional
            the data per element type, the ArrayLike must be of size (number of elements, number of integration points), by default None
        elementFilter : ElementFilter, optional
            the filter to select the element with data
        """
        super(RestrictedIPField,self).__init__(name=name, mesh=mesh,rule=rule,ruleName=ruleName,data=data)
        if elementFilter == None:
            self.elementFilter = ElementFilter(mesh=mesh)
        else:
            self.elementFilter = elementFilter

    def Allocate(self, val:PBasicFloatType=0.) -> None:
        """Allocate the memory to store the data

        Parameters
        ----------
        val : PBasicFloatType, optional
            initial fill value, by default 0.
        """
        self.elementFilter.SetMesh(self.mesh)
        self.data = dict()
        for name, data, ids in self.elementFilter :
            nbItegPoints = len(self.GetRuleFor(name)[1])
            nbElements = len(ids)
            self.data[name] = np.zeros((nbElements,nbItegPoints), dtype=PBasicFloatType)+val

    def GetRestrictedIPField(self, elementFilter:ElementFilter) -> RestrictedIPField:
        """Create a RestrictedIPField only on element intersection between elementFilter and the internal elementFilter

        Parameters
        ----------
        elementFilter : ElementFilter
             ElementFilter to work on

        Returns
        -------
        RestrictedIPField
            a new RestrictedIPField
        """
        res = RestrictedIPField(name=self.name,mesh=self.mesh,rule=self.rule,elementFilter= IntersectionElementFilter(filters=[elementFilter,self.elementFilter]) )
        res.AllocateFromIpField(self)
        return res

    def AllocateFromIpField(self, inputIPField: Union[IPField,RestrictedIPField]) -> None:
        """Fill internal data from a  external ipField (data is copied)

        Parameters
        ----------
        ipField : Union[IPField,RestrictedIPField]
            The IPField to extract the data

        Raises
        ------
        Exception
            In the case the inputIPField is not of the type IPField or RestrictedIPField
        """
        self.name = inputIPField.name
        self.mesh = inputIPField.mesh
        self.rule = inputIPField.rule
        self.elementFilter.SetMesh(self.mesh)
        self.data = dict()

        if isinstance(inputIPField, RestrictedIPField):
            for name,data,ids in self.elementFilter :
                idsII =  inputIPField.elementFilter.GetIdsToTreat(data)
                idds = np.empty(data.GetNumberOfElements(),dtype=PBasicIndexType)
                idds[idsII] = np.arange(len(idsII))
                self.data[name] = inputIPField.data[name][idds[ids],:]
        elif isinstance(inputIPField, IPField):
            for name,data,ids in self.elementFilter :
                self.data[name] = inputIPField.data[name][ids,:]
        else: # pragma: no cover
            raise Exception(f"don't know how to treat ipField of type {type(inputIPField)}")

    def SetDataFromNumpy(self, inData: ArrayLike) ->None:
        """Set the internal data from a one dimensional array

        Parameters
        ----------
        inData : ArrayLike
            the incoming data
        elementFilter : ElementFilter, optional
            ElementFilter to operate, by default (None) all element are treated
        """
        nbValues = 0
        for elemType,data,ids in self.elementFilter:
            nbValues += len(ids)*len(self.GetRuleFor(elemType)[1])

        if inData.size != nbValues: # pragma: no cover
            raise(Exception("incompatible sizes"))

        cpt = 0
        for elemType,data,ids in self.elementFilter:
            localSize = len(ids)*len(self.GetRuleFor(elemType)[1])
            self.data[elemType] = inData[cpt:cpt+localSize].reshape( (len(ids), len(self.GetRuleFor(elemType)[1]))  )
            cpt += localSize

    def GetIpFieldRepresentation(self, fillValue:PBasicFloatType=0.)->IPField:
        """Get a IPFieldRepresentation. Function to expand the representation to he hold domain

        Parameters
        ----------
        fillValue : PBasicFloatType, optional
            field value for cell without data, by default 0.

        Returns
        -------
        IPField
            IPField over the all the element with fillValue on element not covered by the RestrictedIPField
        """
        res = IPField(name=self.name,mesh=self.mesh,rule=self.rule)
        res.Allocate(fillValue)
        for name, data, ids in self.elementFilter :
            if name not in self.data:
                continue
            res.data[name][ids,:] = self.data[name]
        return res

    def CheckCompatibility(self, B: Any) -> None:
        """Check of two field are compatible to make arithmetic operations

        Parameters
        ----------
        B : Any
            the other RestrictedIPField

        Raises
        ------
        Exception
            In the case the fields are not compatible
        """
        super(RestrictedIPField,self).CheckCompatibility(B)
        if isinstance(B,type(self)):
            if not self.elementFilter.IsEquivalent(B.elementFilter):
                raise Exception("The elementFilter of the fields are not the same")

    def unaryOp(self,op:Callable, out:RestrictedIPField=None) -> RestrictedIPField:
        """Internal function to apply a unary operation

        Parameters
        ----------
        op : Callable
            a function to apply to the data
        out : RestrictedIPField, optional
            output IPFiled if non is provided a new one is created, by default None

        Returns
        -------
        RestrictedIPField
            The output RestrictedIPField
        """

        res = type(self)(name = None, mesh=self.mesh,rule=self.rule, elementFilter=self.elementFilter)
        return super(RestrictedIPField,self).unaryOp(op,out=res)

    def binaryOp(self, other:Any ,op: Callable, out:RestrictedIPField=None) -> RestrictedIPField:
        """Internal function to apply a binary operator. A + B for example

        Parameters
        ----------
        other : Any
            a second RestrictedIPField
        op : Callable
            the binary operator to apply
        out : RestrictedIPField, optional
            output RestrictedIPFiled, if None is provided a new one is created, by default None, by default None, by default None

        Returns
        -------
        RestrictedIPField
            The output RestrictedIPField
        """
        res = type(self)(name = None, mesh=self.mesh,rule=self.rule, elementFilter=self.elementFilter)
        return super(RestrictedIPField,self).binaryOp(other,op,out=res)

    def GetCellRepresentation(self, fillValue:PBasicFloatType = 0., method:Union[str,int]='mean') -> np.ndarray:
        if fillValue==0.:
            res = np.zeros(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)
        else:
            res = np.ones(self.mesh.GetNumberOfElements(),dtype=PBasicFloatType)*fillValue

        cpt =0
        self.elementFilter.SetMesh(self.mesh)
        for name, eldata in self.mesh.elements.items():
            ids = self.elementFilter.GetIdsToTreat(eldata)
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
            res[cpt+np.asarray(ids,dtype=PBasicIndexType)]=data
            cpt += nbelems

        return res

def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube([2.,3.,4.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])
    print(mesh)
    sig11 = IPField("Sig_11",mesh=mesh,rule=LagrangeP1)
    print(sig11.GetCellRepresentation())
    sig11.Allocate()
    print(sig11)
    print(sig11.GetPointRepresentation())

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
    methods = ["min","max","mean","maxdiff","maxdifffraction","-1"]
    for method in methods:
        data = sig22.GetCellRepresentation(method=method)

    data = sig22.GetCellRepresentation(1)
    sig22.SetDataFromNumpy(sig22.Flatten())
    data2 = sig22.GetCellRepresentation()

    if np.linalg.norm(data-data2) > 0  :
        raise() # pragma: no cover

    data2 = sig22.GetPointRepresentation(0)

    dummyField = IPField()
    print(str(dummyField))
    dummyField.data[None] = np.arange(3)+1
    print("dummyField")
    print(dummyField*dummyField-dummyField**2/dummyField)

    restrictedIPField = sig22.GetRestrictedIPField(ElementFilter(tag="Skin"))
    restrictedIPField = -(2*restrictedIPField)

    print(restrictedIPField )
    data2 = restrictedIPField.GetPointRepresentation(0)
    restrictedIPField.GetIpFieldRepresentation()
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

    obj = RestrictedIPField(mesh=mesh, data={})
    obj.SetRule()
    print(obj.GetIpFieldRepresentation())
    data = obj.GetCellRepresentation(fillValue=1.)
    NumberOfValues = 0
    for name,data,ids in obj.elementFilter:
        NumberOfValues += len(ids)*len(obj.GetRuleFor(name)[1])
    obj.SetDataFromNumpy(np.zeros(NumberOfValues,dtype=PBasicFloatType))

    methods = ["min","max","mean","maxdiff","maxdifffraction","-1"]
    data = obj.GetCellRepresentation(fillValue=1.)
    for method in methods:
        data = obj.GetCellRepresentation(method=method)


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
