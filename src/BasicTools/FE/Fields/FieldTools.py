# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Union, List, Tuple, Callable, Optional, Iterable, Dict

import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType, ArrayLike
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.Filters import ElementFilter,NodeFilter, IntersectionElementFilter, Filter, FilterOP
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.Fields.IPField import FieldBase, IPField, RestrictedIPField
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo, FESpaceType
from BasicTools.FE.DofNumbering import ComputeDofNumbering
from BasicTools.FE.IntegrationsRules import IntegrationRulesType


def Maximum(A,B):
    return A.binaryOp(B,np.maximum)

def Minimum(A,B):
    return A.binaryOp(B,np.minimum)

def ElementWiseIpToFETransferOp(integrationRule: IntegrationRulesType , space:FESpaceType )-> Dict[str,np.ndarray]:
    """Generate transfer operator (element wise) to pass information from integration points to a FE field. This is done
    by solving a least-squares problem.

    Parameters
    ----------
    integrationRule : IntegrationRuleType
        The integration rule of the original data
    space : FESpaceType
        The target space

    Returns
    -------
    Dict[str,np.ndarray]
        the operator to pass the data from the integration points to the FEField
    """
    res: Dict[str,np.ndarray] = dict()
    for name, ir in integrationRule.items():
        space_ipValues = space[name].SetIntegrationRule(ir[0],ir[1])
        valN = np.asarray( space_ipValues.valN, dtype=PBasicFloatType)
        sol = np.linalg.lstsq(valN, np.eye(valN.shape[0],valN.shape[0]), rcond=None)[0]
        res[name] = sol
    return res

def ElementWiseFEToFETransferOp(originSpace: FESpaceType, targetSpace: FESpaceType)-> Dict[str,np.ndarray]:
    """Generate transfer operator (element wise) to pass information from a FE Field to a different FE field. This is done
    by solving a least-squares problem.

    Parameters
    ----------
    originSpace : FESpaceType
        The initial space
    targetSpace : FESpaceType
        The target space

    Returns
    -------
    Dict[str,np.ndarray]
        the operator to pass the data from the initial space to the target FEField
    """
    res: Dict[str,np.ndarray] = dict()
    for name, ir in targetSpace.items():
        ir = (originSpace[name].posN, np.ones(originSpace[name].GetNumberOfShapeFunctions(),dtype=PBasicFloatType) )
        spaceIpValues = targetSpace[name].SetIntegrationRule(ir[0],ir[1])
        valN = np.asarray( spaceIpValues.valN, dtype=PBasicFloatType)
        sol = np.linalg.lstsq(valN, np.eye(valN.shape[0],valN.shape[1]), rcond=None)[0]
        res[name] = sol
    return res

def NodeFieldToFEField(mesh: UnstructuredMesh, nodeFields: Dict[str,ArrayLike]=None) -> Dict[str,FEField]:
    """Create FEField(P isoparametric) from the node field data.
    if nodesFields is None the mesh.nodeFields is used

    Parameters
    ----------
    mesh : UnstructuredMesh
        The support for the FEFields
    nodeFields : Dict[str,ArrayLike], optional
        the dictionary containing the nodes fields to be converted to FEFields,
        if None the mesh.nodeFields is used, by default None

    Returns
    -------
    Dict[str,FEField]
        A dictionary the keys are field names and the values are the FEFields
    """


    if nodeFields is None:
        nodeFields = mesh.nodeFields
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)
    res = {}
    for name,values in nodeFields.items():
        if len(values.shape) == 2 and values.shape[1] > 1:
            for i in range(values.shape[1]):
                res[name+"_"+str(i)] = FEField(name=name+"_"+str(i), mesh=mesh, space=LagrangeSpaceGeo, numbering=numbering, data=np.asarray(values[:,i], dtype=PBasicFloatType, order="C") )
        else:
            res[name] = FEField(name=name, mesh=mesh, space=LagrangeSpaceGeo, numbering=numbering,data=np.asarray(values, dtype=PBasicFloatType, order="C"))
    return res

def ElemFieldsToFEField(mesh: UnstructuredMesh, elemFields: Dict[str,ArrayLike]=None) -> Dict[str,FEField]:
    """Create FEField(P0) from the elements field data.
    if elemFields is None the mesh.elemFields is used

    Parameters
    ----------
    mesh : UnstructuredMesh
        The support for the FEFields
    elemFields : Dict[str,ArrayLike], optional
        a dict containing the elements fields to be converted to FEFields
        if elemFields is None the mesh.elemFields is used, by default None

    Returns
    -------
    Dict[str,FEField]
        A dictionary the keys are field names and the values are the FEFields
    """

    if elemFields is None:
        elemFields = mesh.elemFields
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0)
    res = {}
    for name,values in elemFields.items():
        if len(values.shape) == 2 and values.shape[1] > 1:
            for i in range(values.shape[1]):
                res[name+"_"+str(i)] = FEField(name=name+"_"+str(i), mesh=mesh, space=LagrangeSpaceP0, numbering=numbering, data=np.asarray(values[:,i], dtype=PBasicFloatType, order="C"))
        else:
            res[name] = FEField(name=name, mesh=mesh, space=LagrangeSpaceP0, numbering=numbering, data=np.asarray(values, dtype=PBasicFloatType, order="C"))
    return res

def FEFieldsDataToVector(listOfFields: List[FEField], outVector: Optional[ArrayLike]=None) -> np.ndarray:
    """Put FEFields data into a vector format

    Parameters
    ----------
    listOfFields : List[FEField]
        list of FEFields (the order is important)
    outVector : Optional[ArrayLike], optional
        if provided the output vector, by default None

    Returns
    -------
    np.ndarray
        _description_
    """

    if outVector is None:
        s = sum([f.numbering.size for f in listOfFields])
        outVector = np.zeros(s,dtype=PBasicFloatType)
    offset = 0
    for f in listOfFields:
        outVector[offset:offset+f.numbering.size] = f.data
        offset += f.numbering.size
    return outVector

def VectorToFEFieldsData(inVector: ArrayLike, listOfFields: List[FEField]):
    """Put vector data in FEFields

    Parameters
    ----------
    inVector : ArrayLike
        vector with the data
    listOfFields : List[FEField]
        list of FEFields to store the data
    """

    offset = 0
    for f in listOfFields:
        f.data = inVector[offset:offset+f.numbering.size]
        offset += f.numbering.size

def GetPointRepresentation(listOfFEFields: List[FEField], fillValue: PBasicFloatType=0.) -> np.ndarray:
    """Get a np.ndarray compatible with the points of the mesh
    for a list of FEFields. Each field in the list is treated
    as a scalar component for the output field

    Parameters
    ----------
    listOfFEFields : List[FEField]
        list of FEFields to extract the information
    fillValue : PBasicFloatType, optional
        value to use in the case some field are not defined
        over all the nodes, by default 0.

    Returns
    -------
    np.ndarray
        the output and nd array of size (nb points, nb of FEField)
    """

    nbFields= len(listOfFEFields)
    res = np.empty((listOfFEFields[0].mesh.GetNumberOfNodes(), nbFields) )
    for i,f in enumerate(listOfFEFields):
        res[:,i] = f.GetPointRepresentation(fillvalue= fillValue)
    return res

def GetCellRepresentation(listOfFEFields: List[FEField], fillValue: PBasicFloatType=0., method="mean") -> np.ndarray:
    """ Get a np.ndarray compatible with the points of the mesh
    for a list of FEFields. Each field in the list is treated
    as a scalar component for the output field

    Parameters
    ----------
    listOfFEFields : list
        list of FEFields to extract the information
    fillValue: float
        value to use in the case some field are not defined
        over all the nodes, by default 0.
    method: str | mean
        Read FEField.GetCellRepresentation for more information, by default "mean".

    Returns
    -------
    np.array
        Array  of size (number of elements, number of fields)

    """
    nbFields= len(listOfFEFields)
    res = np.empty((listOfFEFields[0].mesh.GetNumberOfElements(), nbFields))
    for i,f in enumerate(listOfFEFields):
        res[:,i] = f.GetCellRepresentation(fillvalue=fillValue, method=method)
    return res

class IntegrationPointWrapper(FieldBase):
    """ Class to generate a FEField at the integration points
        Two important function are available : diff and GetIpField.


    Parameters
    ----------
    field : FEField
        the field to evaluate
    rule : _type_
        the integration rule
    elementFilter : ElementFilter, optional
        in this case a RestrictedIPField is generated, by default None

    """

    def __init__(self, field: FEField, rule, elementFilter: Optional[ElementFilter]=None):
        if not isinstance(field,FEField):
            raise(Exception("IntegrationPointWrapper work only on FEFields"))

        self.feField = field
        self.rule = rule
        self._ipCache = None
        self._diffIpCache = {}
        self.elementFilter = elementFilter

    @property
    def name(self):
        return self.feField.name

    def ResetCache(self):
        self._ipCache = None
        self._diffIpCache = {}

    def diff(self, compName:str) -> Union[IPField, RestrictedIPField] :
        """Get the IPField representation of this field with a derivation in the direction of compName

        Parameters
        ----------
        compName : str
            compName are available in : from BasicTools.FE.SymWeakForm import space

        Returns
        -------
        (IPField | RestrictedIPField)
            Depending on the presence of an elementFilter
        """
        from BasicTools.FE.SymWeakForm import space

        for cm in range(3):
            if space[cm] == compName:
                break
        else:
            cm = compName

        if cm not in self._diffIpCache:
            from BasicTools.FE.Fields.FieldTools import TransferFEFieldToIPField
            self._diffIpCache[cm] = TransferFEFieldToIPField(self.feField,der=cm,rule=self.rule, elementFilter= self.elementFilter)
        return self._diffIpCache[cm]

    def GetIpField(self)-> Union[IPField, RestrictedIPField]:
        """Get the integration point representation of the field

        Returns
        -------
        Union[IPField, RestrictedIPField]
            Depending on the presence of an elementFilter
        """

        if self._ipCache is None:
            from BasicTools.FE.Fields.FieldTools import TransferFEFieldToIPField
            self._ipCache =  TransferFEFieldToIPField(self.feField,der=-1,rule=self.rule, elementFilter= self.elementFilter)
        return self._ipCache

    def unaryOp(self,op):
        return self.GetIpField().unaryOp(op)

    def binaryOp(self,other,op):
        return self.GetIpField().binaryOp(other,op)

    @property
    def data(self):
        return self.GetIpField().data

DescriptionValue = Union[float,Callable[[np.ndarray],float]]
DescriptionUnit = Tuple[Filter,DescriptionValue]

def CreateFieldFromDescription(mesh: UnstructuredMesh,
                               fieldDefinition: Iterable[DescriptionUnit],
                               fieldType: str = "FE") -> Union[FEField,IPField]:
    """
    Create a FEField from a field description

    Parameters
    ----------
    mesh : UnstructuredMesh
    fieldDefinition : List[Tuple[Union[ElementFilter,NodeFilter],Union[float,Callable[[np.ndarray],float]]]]
        A field definition is a list of Tuples (with 2 element each ).
        The first element of the tuple is a ElementFilter or a NodeFilter
        The second element of the tuple is a float or float returning callable (argument is a point in the space)
    fieldType : str, optional
        type of field created FEField or IPField
        ("FE", "FE-P0", "IP") by default "FE" (isoparametric)

    Returns
    -------
    Union[FEField,IPField]
        created Field with values coming from the field description
    """
    if fieldType == "FE":
        from BasicTools.FE.FETools import PrepareFEComputation
        spaces,numberings,offset, NGauss = PrepareFEComputation(mesh,numberOfComponents=1)
        res = FEField(mesh=mesh,space=spaces,numbering=numberings[0])
        res.Allocate()
        FillFEField(res,fieldDefinition)

    elif fieldType == "FE-P0":
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
        numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0)
        res = FEField(mesh=mesh,space=LagrangeSpaceP0,numbering=numbering)
        res.Allocate()
        FillFEField(res,fieldDefinition)

    elif fieldType == "IP":
        res = IPField(mesh=mesh,ruleName="LagrangeIsoParam")
        res.Allocate()
        FillIPField(res,fieldDefinition)

    else:
        raise(Exception("fieldType not valid"))

    return res

def GetTransferOpToIPField(inField: FEField,
                           ruleName: Optional[str] = None,
                           rule: Optional[Tuple[np.ndarray,np.ndarray]]=None,
                           der: int=-1,
                           elementFilter:Optional[ElementFilter]=None) -> FeToIPOp:
    """Compute the transfer operator (.dot()) to construct a IPField (of RestrictedIPField)
        IPField = FeToIPOp.dot(FEField)

    Parameters
    ----------
    inField : FEField
        the input FEField
    ruleName : Optional[str], optional
        the ruleName of the integration rule to use, by default None
        (see BasicTools.FE.IntegrationsRules:GetRule for more info)
    rule : Optional[Tuple[np.ndarray,np.ndarray]], optional
        the integration rule to use, by default None
        (see BasicTools.FE.IntegrationsRules:GetRule for more info)
    der : int, optional
        the coordinate to be derivated
        -1 no derivations only evaluation at IP
        [0,1,2] the coordinate number to compute derivative of the FEField
        by default -1
    elementFilter : Optional[ElementFilter], optional
        An element filter to restrict the operator.
        In this case a RestrictedIPField is generated
        _description_, by default None

    Returns
    -------
    FeToIPOp
        An instance of an object with the .dot operator
    """
    from BasicTools.FE.IntegrationsRules import GetRule
    from BasicTools.FE.Spaces.IPSpaces import GenerateSpaceForIntegrationPointInterpolation
    from BasicTools.FE.Integration import IntegrateGeneral
    from BasicTools.FE.DofNumbering import ComputeDofNumbering

    iRule = GetRule(ruleName=ruleName,rule=rule)
    gaussSpace = GenerateSpaceForIntegrationPointInterpolation(iRule)
    mesh = inField.mesh

    class FeToIPOp(dict):
        """Helper class to hold the transfer matrices between the FEField and the IPField

        Parameters
        ----------
        iRule :
            The integration rule to be used
        elementFilter : ElementFilter
            The ElementFilter to generate the RestrictedIPField
        """
        def __init__(self,iRule,elementFilter):
            super(FeToIPOp,self).__init__()
            self.iRule = iRule
            self.elementFilter = elementFilter

        def dot(self, inField: FEField) -> Union[IPField,RestrictedIPField]:
            """Implementation of the transfer

            Parameters
            ----------
            inField : FEField
                input FEField

            Returns
            -------
            Union[IPField,RestrictedIPField]
                the type depends on the presence or not of a self.elementFilter
            """
            return TransferFEFieldToIPField(inField,rule=self.iRule,elementFilter=self.elementFilter,op=self)

    res = FeToIPOp(iRule,elementFilter=elementFilter)

    for elemType,d in mesh.elements.items():

        if elementFilter is not None:
            eF = IntersectionElementFilter(mesh,(ElementFilter(inField.mesh,elementTypes=[elemType]) ,elementFilter) )
        else:
            eF = ElementFilter(inField.mesh,elementTypes=[elemType])

        idsToTreat = eF.GetIdsToTreat(d)
        if len(idsToTreat) == 0:
            continue

        numberingRight = ComputeDofNumbering(mesh, Space=gaussSpace, elementFilter=eF)

        rightField = FEField(name="Gauss'",numbering=numberingRight,mesh=mesh,space=gaussSpace)

        from BasicTools.FE.SymWeakForm import GetField,GetTestField,space
        LF = GetField(inField.name,1)
        RF = GetTestField("Gauss",1)

        if der == -1:
            symForm = LF.T*RF
        else:
            symForm = LF.diff(space[der]).T*RF

        transferMatrixMatrix,_ = IntegrateGeneral(mesh=mesh,constants={},fields=[],wform=symForm,
                                                unkownFields= [inField],testFields=[rightField],
                                                onlyEvaluation=True,integrationRule=iRule,
                                                elementFilter=eF)
        res[elemType] = transferMatrixMatrix

    return res

def TransferFEFieldToIPField(inFEField: FEField, ruleName:str=None, rule=None,der:int=-1, elementFilter:Optional[ElementFilter]=None, op =None)-> Union[IPField,RestrictedIPField]:
    """Transfer FEField to a IPField

    Parameters
    ----------
    inFEField : FEField
        the input FEField
    ruleName : Optional[str], optional
        the ruleName of the integration rule to use, by default None
        (see BasicTools.FE.IntegrationsRules:GetRule for more info)
    rule : Optional[Tuple[np.ndarray,np.ndarray]], optional
        the integration rule to use, by default None
        (see BasicTools.FE.IntegrationsRules:GetRule for more info)
    der : int, optional
        the coordinate to be derivated
        -1 no derivations only evaluation at IP
        [0,1,2] the coordinate number to compute derivative of the FEField
        by default -1
    elementFilter : Optional[ElementFilter], optional
        An element filter to restrict the operator.
        In this case a RestrictedIPField is generated
        _description_, by default None
    op : _type_, optional
        an object returned by the GetTransferOpToIPField function
        if not given then this op is calculated internally
        by default None

    Returns
    -------
    Union[IPField,RestrictedIPField]
        The field evaluated at the integration points
    """
    if op is None:
        op = GetTransferOpToIPField(inField=inFEField, ruleName=ruleName, rule=rule, der=der, elementFilter=elementFilter)

    from BasicTools.FE.IntegrationsRules import GetRule

    iRule = GetRule(ruleName=ruleName,rule=rule)
    if elementFilter is None:
        outField = IPField(name=inFEField.name,mesh=inFEField.mesh,rule=iRule)
        outField.Allocate()
    else:
        outField = RestrictedIPField(name=inFEField.name, mesh=inFEField.mesh, rule=iRule, elementFilter=elementFilter)
        outField.Allocate()

    mesh = inFEField.mesh
    for elemType,d in mesh.elements.items():
        if elemType not in op :
            continue
        nbElements = op[elemType].shape[0]//len(iRule[elemType][1])

        data = op[elemType].dot(inFEField.data)
        outField.data[elemType] = np.reshape(data,(nbElements, len(iRule[elemType][1]) ),'F')
        outField.data[elemType] = np.ascontiguousarray(outField.data[elemType])

    return outField

def TransferPosToIPField(mesh: UnstructuredMesh,
                        ruleName: Optional[str] = None,
                        rule: Optional[Tuple[np.ndarray,np.ndarray]] = None,
                        elementFilter: Optional[ElementFilter]= None) -> List[IPField] :
    """
    Create (2 or 3) IPFields with the values of the space coordinates


    Parameters
    ----------
    mesh : UnstructuredMesh
    ruleName : str, optional
        The Integration rule name, by default None ("LagrangeIsoParam")
    rule : Tuple[np.ndarray,np.ndarray], optional
        The Integration Rule, by default None ("LagrangeIsoParam")
    elementFilter : ElementFilter, optional
        the zone where the transfer must be applied, by default None (all the elements)

    Returns
    -------
    List[IPField]
        The IPFields containing the coordinates
    """

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)

    inField = FEField(mesh=mesh,space=LagrangeSpaceGeo,numbering=numbering,data=None)

    op = GetTransferOpToIPField(inField=inField,ruleName=ruleName,rule=rule,der=-1,elementFilter=elementFilter)

    d = mesh.nodes.shape[1]
    outField = []

    for c,name in enumerate(["posx","posy","posz"][0:d]):
        inField.data = mesh.nodes[:,c]
        f = op.dot(inField)
        f.name = name
        outField.append(f)
    return outField

def FillIPField(field : IPField, fieldDefinition) -> None:
    """
    Fill a IPField using a definition

    Parameters
    ----------
    field : IPField
        IPField to treat
    fieldDefinition : List[Tuple[ElementFilter,Union[float,Callable[[np.ndarray],float]]]]
        A field definition is a list of Tuples (with 2 element each ).
        The first element of the tuple is a ElementFilter or a NodeFilter
        The second element of the tuple is a float or float returning callable (argument is a point in the space)
    """

    for f,val in fieldDefinition:
        if callable(val):
            fValue = val
        else:
            fValue = lambda x: val

        f.mesh = field.mesh
        if isinstance(f,ElementFilter):
            for name,elements,ids in f:
                geoSpace = LagrangeSpaceGeo[name]
                rule = field.GetRuleFor(name)
                nbIp = len(rule[1])
                geoSpaceIpValues = geoSpace.SetIntegrationRule(rule[0],rule[1] )

                for elementsId in ids:
                    for i in range(nbIp):
                        valN = geoSpaceIpValues.valN[i]
                        elementNodesCoordinates = field.mesh.nodes[elements.connectivity[elementsId,:],:]
                        pos = np.dot(valN ,elementNodesCoordinates).T
                        field.data[name][elementsId,i] = fValue(pos)
            continue
        raise ValueError(f"Cant use this type of filter to fill an IPField : {type(f)}")

def FillFEField(field : FEField, fieldDefinition) -> None:
    """ Fill a FEField using a definition

    Parameters
    ----------
    field : FEField
        FEField to treat
    fieldDefinition : List[Tuple[Union[ElementFilter, NodeFilter], Union[float, Callable[[np.ndarray], float]]]]
        A field definition is a list of Tuples (with 2 element each ).
        The first element of the tuple is a ElementFilter or a NodeFilter
        The second element of the tuple is a float or float returning callable (argument is a point in the space)
    """
    for f,val in fieldDefinition:
        needPos = False
        if callable(val):
            fValue = val
            needPos = True
        else:
            fValue = lambda x: val

        f.mesh = field.mesh
        if isinstance(f,(ElementFilter,FilterOP)):
            treated = np.zeros(field.numbering.size, dtype=bool)
            nodes = field.mesh.nodes
            for name,elements,ids in f:
                numbering = field.numbering[name]
                if needPos == False:
                    idDofs = numbering[ids,:].flatten()
                    field.data[idDofs] = fValue(None)
                else:
                    sp = field.space[name]
                    nbShapeFunctions = sp.GetNumberOfShapeFunctions()
                    geoSpace = LagrangeSpaceGeo[name]
                    geoSpaceIpValues  = geoSpace.SetIntegrationRule(sp.posN,np.ones(nbShapeFunctions) )
                    valN = geoSpaceIpValues.valN
                    connectivity = elements.connectivity
                    for elementId in ids:
                        elementNodesCoordinates = nodes[connectivity[elementId,:],:]
                        dofs  = numbering[elementId,:]
                        pos = np.dot(valN ,elementNodesCoordinates)
                        for i in range(nbShapeFunctions):
                            if treated[dofs[i]]:
                                continue
                            field.data[dofs[i]] = fValue(pos[i,:])
                            treated[dofs[i]] = True
                    #for i in range(nbShapeFunctions):
                    #    dofId = field.numbering[name][elementId,i]
                    #    if treated[dofId]:
                    #        continue
                    #    treated[dofId] =True
                    #    valN = geoSpaceIpValues.valN[i]
                    #    pos = np.dot(valN ,elementNodesCoordinates).T
                    #    field.data[dofId] = fValue(pos)

        elif isinstance(f,NodeFilter):
            ids = f.GetIdsToTreat()
            for pid in ids:
                dofId = field.numbering.GetDofOfPoint(pid)
                pos = field.mesh.nodes[pid,:]
                field.data[dofId] = fValue(pos)
        else:
            raise ValueError(f"Don't know how to treat this type of field {type(f)}")

def FieldsAtIp(listOfFields: List[Union[FEField,IPField,RestrictedIPField]], rule, elementFilter:Optional[ElementFilter]=None) -> List[Union[IntegrationPointWrapper,IPField,RestrictedIPField]]:
    """Convert a list of Field (FEFields,IPField,RestrictedIPField) to a list of IPFields

    Parameters
    ----------
    listOfFields : List[Union[IPField,RestrictedIPField]]
        the list of fields to be treated
    rule : _type_
        the integration rule to be used for the evaluation
    elementFilter : Optional[ElementFilter], optional
        the filter to select the element to treat, by default None

    Returns
    -------
    List[Union[IntegrationPointWrapper,IPField,RestrictedIPField]]
        the list of IPFilters

    Raises
    ------
    Exception
        in the case of a internal error
    """
    from BasicTools.FE.Fields.FieldTools import IntegrationPointWrapper
    res = []
    for f in listOfFields:
        if isinstance(f,FEField):
            res.append(IntegrationPointWrapper(f, rule, elementFilter=elementFilter))
        elif isinstance(f,IPField):
            if f.rule == rule:
                if elementFilter is None:
                    res.append(f)
                else:
                    res.append(f.GetRestrictedIPField(elementFilter) )
            else:
                print (f.rule)
                print (rule)
                print(f"skipping IPField {f.name} because it has a not compatible IP rule type {str(type(f))}")
        else:
            raise Exception(f"Don't know how to treat this type of field {type(f)}")
    return res

class FieldsMeshTransportation():
    """Class to help the transfer of field (FEField and IPField) between old and new meshes
    """
    def __init__(self):
        self.cache_numbering = {}

    def ResetCacheData(self):
        self.cache_numbering = {}

    def GetNumbering(self, mesh, space, fromConnectivity=False,discontinuous=False):
        meshId = str(id(mesh))
        spaceId = str(id(space))
        key = (meshId,spaceId,fromConnectivity,discontinuous)

        if key not in self.cache_numbering:
            self.cache_numbering[key] = ComputeDofNumbering(mesh, space, fromConnectivity=True,discontinuous=False)
        return self.cache_numbering[key]

    def TransportFEFieldToOldMesh(self, oldMesh: UnstructuredMesh, inFEField: FEField, fillValue: PBasicFloatType=0.) -> FEField:
        """Function to define a FEField on the old mesh, the inFEField mesh must be a
        transformation of the old mesh. This means the infield mesh originalIds
        (for nodes and elements) must be with respect to the old mesh.

        Parameters
        ----------
        oldMesh : UnstructuredMesh
            the old mesh
        inFEField : FEField
            a FEField
        fillValue : PBasicFloatType, optional
            Value to fill the values of the dofs not present in the inFEField , by default 0.

        Returns
        -------
        FEField
            a FEField on the old mesh filled with the values of the inFEField
        """

        if id(inFEField.mesh) == id(oldMesh):
            return inFEField

        name = inFEField.name
        space = inFEField.space
        if inFEField.numbering.fromConnectivity:
            numbering = ComputeDofNumbering(oldMesh, space, fromConnectivity=True)
            res = FEField(name = name, mesh=oldMesh, space=space, numbering=numbering)
            res.Allocate(fillValue)
            res.data[inFEField.mesh.originalIDNodes] = inFEField.data
        else :
            numbering =self.GetNumbering(oldMesh,space)
            res = FEField(name = name, mesh=oldMesh, space=space, numbering=numbering)
            res.Allocate(fillValue)
            for name,data in oldMesh.elements.items():
                if name  not in inFEField.mesh.elements:
                    continue
                newData = inFEField.mesh.elements[name]
                res.data[numbering[name][newData.originalIds,:].flatten()]= inFEField.data[inFEField.numbering[name].flatten()]
        return res

    def TransportFEFieldToNewMesh(self, inFEField: FEField, newMesh: UnstructuredMesh)-> FEField:
        """Function to define a FEField on the new mesh, the new mesh must be a
        transformation of the mesh in the infield. This means the new mesh originalIds
        (for nodes and elements) must be with respect to the mesh of the infield

        Parameters
        ----------
        inFEField : FEField
            the FEField to be transported to the new mesh
        newMesh : UnstructuredMesh
            the target mesh

        Returns
        -------
        FEField
            a FEField defined on the newMesh
        """

        if id(inFEField.mesh) == id(newMesh):
            return inFEField

        name = inFEField.name
        space = inFEField.space
        if  inFEField.numbering.fromConnectivity:
            numbering = ComputeDofNumbering(newMesh,space,fromConnectivity=True)
            res = FEField(name = name, mesh=newMesh, space=space, numbering=numbering)
            res.data = inFEField.data[newMesh.originalIDNodes]
        else:
            numbering = ComputeDofNumbering(newMesh, space)
            res = FEField(name = name, mesh=newMesh, space=space, numbering=numbering)
            res.Allocate()
            for name,data in newMesh.elements.items():
                res.data[numbering.numbering[name].flatten()] = inFEField.data[inFEField.numbering[name][data.originalIds,:].flatten()]
        return res

    def TransportIPFieldToOldMesh(self, oldMesh:UnstructuredMesh, inIPField:IPField)->IPField:
        """function to define a IPField on the old mesh, the old mesh must be the origin
        of a transformation of the mesh in the IPField. This means the IPField mesh originalIds
        (for nodes and elements) must be with respect to the old mesh

        Parameters
        ----------
        oldMesh : UnstructuredMesh
            the old mesh
        inIPField : IPField
            the IPField to be transported to the old mesh

        Returns
        -------
        IPField
            a IPField defined on the old Mesh
        """
        if id(inIPField.mesh) == id(oldMesh):
            return inIPField

        outData = {}
        for elemType,data in inIPField.mesh.elements.items():
            inData = inIPField.data[elemType]
            outData[elemType] = np.zeros((oldMesh.elements[elemType].GetNumberOfElements(),inData.shape[1]))
            outData[elemType][data.originalIds,:] = inIPField.data[elemType]
        res = IPField(name=inIPField.name,mesh=oldMesh,rule=inIPField.rule,data=outData)
        return res

    def TransportIPFieldToNewMesh(self, inIPField:IPField, newMesh:UnstructuredMesh)->IPField:
        """Function to define a IPField on the new mesh, the new mesh must be a
        transformation of the mesh in the infield. This means the new mesh originalIds
        (for nodes and elements) must be with respect to the mesh of the infield

        Parameters
        ----------
        inIPField : IPField
            the FEField to be transported to the new mesh
        newMesh : UnstructuredMesh
            the target mesh

        Returns
        -------
        IPField
            a IPField defined on the newMesh
        """
        if id(inIPField.mesh) == id(newMesh):
            return inIPField

        outputData = {}
        for elemType,data in newMesh.elements.items():
            outputData[elemType] = inIPField.data[elemType][data.originalIds,:]
        res = IPField(name=inIPField.name,mesh=newMesh,rule=inIPField.rule,data=outputData)
        return res

class FieldsEvaluator():
    """helper to evaluate expression using FEField and IPFields
    The user can add fields and constants
    Then this class can be used to compute an expression (given using an callable object)
    the callable can capture some of the argument but the last must be **args (so extra
    argument are ignored).

    """
    def __init__(self,fields=None):
        self.originals = {}
        from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac
        self.rule = IntegrationRulesAlmanac['LagrangeIsoParam']
        self.ruleAtCenter = IntegrationRulesAlmanac['ElementCenterEval']
        self.atIp = {}
        self.atCenter = {}
        if fields is not None:
            for f in fields:
                self.AddField(f)
        self.constants = {}
        self.elementFilter = None
        self.modified = True

    def AddField(self, field: FieldBase):
        """Add a field to the internal almanac

        Parameters
        ----------
        field : FieldBase
            A field to be added to the almanac
        """
        self.originals[field.name] = field

    def AddConstant(self, name: str, val: PBasicFloatType):
        """Add a constant value to the almanac

        Parameters
        ----------
        name : str
            name of the value
        val : PBasicFloatType
            Value
        """
        self.constants[name] = val

    def Update(self, what: str="all"):
        """Function to update the field calculated at different support
        (a FEField can be internally evaluated at the integration points)
        The objective of this function is to update this evaluation

        Parameters
        ----------
        what : str, optional
            ["all", "IPField", "Centroids"], by default "all"
        """
        for name,field in self.originals.items():
            if what=="all" or what =="IPField":
                self.atIp[name] = FieldsAtIp([field], self.rule, elementFilter=self.elementFilter)[0]

            if what=="all" or what =="Centroids":
                self.atCenter[name] = FieldsAtIp([field], self.ruleAtCenter, elementFilter=self.elementFilter)[0]

    def GetFieldsAt(self, on:str) -> Dict:
        """Get the internal representations of the user fields ant different support

        Parameters
        ----------
        on : str
            ["IPField", "Centroids", "FEField", "Nodes"]

        Returns
        -------
        Dict
            A dictionary containing the different representation of the user fields

        Raises
        ------
        ValueError
            In the case the "on" parameter is not valid
        """

        if on == "IPField":
            res =  self.atIp
        elif on == "Nodes":
            res = {}
            for f in self.originals.values():
                res[f.name] = f.GetPointRepresentation()
        elif on == "Centroids":
            res =  self.atCenter
        elif on == "FEField":
            res = self.originals
        else:
            raise ValueError(f"Target support not supported ({on})" )
        result = dict(self.constants)
        result.update(res)
        return result

    def GetOptimizedFunction(self,func):
        import sympy
        import inspect
        from BasicTools.FE.SymWeakForm import space

        class OptimizedFunction():
            """ Internal function to convert a function written in python into a
            compiled sympy function.

            """
            def __init__(self,func,constants):
                self.constants = dict(constants)
                # get arguments names
                args = inspect.getfullargspec(func).args
                # clean self, in the case of a member function
                if args[0] == "self":
                    args.pop(0)

                symbolicSymbols = {}

                self.args = args
                for n  in args:
                    if n in self.constants:
                        symbolicSymbols[n] = self.constants[n]
                    else:
                        symbolicSymbols[n] = sympy.Function(n)(*space)
                #symbolically evaluation of the function
                funcValues = func(**symbolicSymbols)

                repl, residual = sympy.cse(funcValues)

                restKeys = list(symbolicSymbols.keys())
                restKeys.extend(["x","y","z"])
                self.auxFunctions = []
                self.auxNames = []
                for i, v in enumerate(repl):
                    funcLam = sympy.lambdify(restKeys,v[1])
                    self.auxFunctions.append(funcLam)
                    funcName = str(v[0])
                    self.auxNames.append(funcName)
                    restKeys.append(str(v[0]))
                self.mainFunc = sympy.lambdify(restKeys,residual)

            def __call__(self,**args):
                numericFields = {"x":0,"y":0,"z":0}

                class Toto():
                    def __init__(self,obj):
                        self.obj=obj
                    def __call__(self,*args,**kwds):
                        return self.obj

                for n  in self.args:
                    if n not in args:
                        numericFields[n] = self.constants[n]
                    else:
                        numericFields[n] = Toto(args[n])

                for name,func in zip(self.auxNames, self.auxFunctions):
                    numericFields[name] = func(**numericFields)
                return self.mainFunc(**numericFields)[0]

        return OptimizedFunction(func,self.constants)

    def Compute(self, func: Callable, on:str, useSympy:bool=False):
        """Compute the function func using the user fields on the support "on"


        Parameters
        ----------
        func : Callable
            the function to evaluate
        on : str
             ["IPField", "Centroids", "FEField", "Nodes"]
        useSympy : bool, optional
            if true the function is compiled using sympy , by default False

        Returns
        -------
        Any
            The callable evaluated
        """

        if useSympy:
            func = self.GetOptimizedFunction(func)

        fields = self.GetFieldsAt(on)
        return func(**fields)

def ComputeTransferOp(field1: FEField, field2: FEField, force: bool = False)-> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the transfer Operator from field1 to field2
    Both fields must be compatibles fields: same mesh, same space
    but with different numberings

    - example of use:

        lIndex,rIndex = ComputeTransferOp(field1,field2):

        field2.data[lIndex] = field1.data[rIndex]

    Parameters
    ----------
    field1 : FEField
        field from where the information is extracted
    field2 : FEField
        field to receive the information
    force : bool, optional
        true to bypass the compatibility check, by default False

    Returns
    -------
    Tuple(np.ndarray,np.ndarray)
        index vector for field2, index vector for field1
    """

    if field1.mesh != field2.mesh and not force:
        raise(Exception("The fields must share the same mesh"))

    if field1.space != field2.space and not force:
        raise(Exception("The fields must share the same space"))

    _extractorLeft  = []
    _extractorRight = []

    for name in field1.mesh.elements:
        a = field1.numbering.get(name,None)
        b = field2.numbering.get(name,None)
        if a is not None and b is not None:
            mask = np.logical_and(a>=0,b>=0)
            _extractorLeft.extend(b[mask])
            _extractorRight.extend(a[mask])

    extractorLeft = np.array(_extractorLeft, dtype=PBasicIndexType)
    extractorRight = np.array(_extractorRight, dtype=PBasicIndexType)

    v,index = np.unique(extractorLeft, return_index=True)
    extractorLeft = v
    extractorRight = extractorRight[index]

    return extractorLeft, extractorRight

def CheckIntegrity(GUI=False):

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars

    mesh =CreateUniformMeshOfBars(0, 1, 10)
    bars = mesh.GetElementsOfType(EN.Bar_2)
    bars.tags.CreateTag("first3").SetIds([0,1,2])
    bars.tags.CreateTag("next3").SetIds([3,4,5])
    bars.tags.CreateTag("Last").SetIds([8])
    mesh.nodesTags.CreateTag("FirstPoint").SetIds([0])
    mesh.nodesTags.CreateTag("LastPoint").SetIds([9])
    print(mesh)

    fieldDefinition =  [ (ElementFilter(), 5) ]
    fieldDefinition.append( (ElementFilter(tags=["first3" ]), 3) )
    fieldDefinition.append( (ElementFilter(tags=["Last" ]), -1) )
    fieldDefinition.append( (ElementFilter(tags=["next3" ]), lambda x : x[0]) )

    field = CreateFieldFromDescription(mesh, fieldDefinition, fieldType="IP" )
    print(field)
    print(field.data)

    fieldDefinition.append( (NodeFilter(tags=["FirstPoint" ]), -10) )
    fieldDefinition.append( (NodeFilter(tags=["LastPoint" ]), lambda x : x[0]+1.2) )

    field = CreateFieldFromDescription(mesh, fieldDefinition )
    print(field.data)
    field.Allocate(val=0)
    FillFEField(field, fieldDefinition )


    print("--*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
    nodalTransferredField = TransferFEFieldToIPField(field,ruleName="LagrangeIsoParam",elementFilter=ElementFilter(dimensionality=1,tag="next3"))
    print("--*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--")
    print("input")
    print(field.data)
    print("output")
    print(nodalTransferredField.data )

    print(nodalTransferredField.GetIpFieldRepresentation(1).data)

    res = TransferPosToIPField(mesh,ruleName="LagrangeIsoParam",elementFilter=ElementFilter(dimensionality=1,tag="next3"))
    print(res[0].data)
    print(res[0].GetIpFieldRepresentation().data)

    FE = FieldsEvaluator()
    field.name = "E"
    FE.AddField(field)
    FE.AddConstant("alpha",0)

    def op(E,alpha,**args):
        return (2*E+1+alpha)+(2*E+1+alpha)**2+(2*E+1+alpha)**3

    print("--------")
    print(field.data)
    print("--------")
    res = FE.Compute(op,"FEField")
    print(res.data)

    FE.Update("IPField")
    res = FE.Compute(op,"IPField")
    print(res.data)


    res = FE.Compute(op,"FEField", useSympy=True)
    print(res.data)
    mesh.nodeFields["FEField_onPoints"] = GetPointRepresentation((res,))

    GetCellRepresentation((res,))

    print(NodeFieldToFEField(mesh))
    print(ElemFieldsToFEField(mesh))

    vector = FEFieldsDataToVector([res])

    VectorToFEFieldsData(vector,[res])
    res = FE.GetOptimizedFunction(op)

    obj = FieldsMeshTransportation()


    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True)) #pragma no cover
