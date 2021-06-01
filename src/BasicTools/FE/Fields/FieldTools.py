# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Containers.Filters import ElementFilter,NodeFilter, IntersectionElementFilter
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.Fields.IPField import IPField, RestrictedIPField
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.DofNumbering import ComputeDofNumbering
from BasicTools.FE.Fields.IPField import FieldBase


def NodeFieldToFEField(mesh):
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)
    res = {}
    for name,values in mesh.nodeFields.items():
        if len(values.shape) == 2:
            for i in range(values.shape[1]):
                res[name+"_"+str(i)] = FEField(name=name+"_"+str(i), mesh=mesh, space=LagrangeSpaceGeo, numbering=numbering,data=values[:,i])
        else:
            res[name] = FEField(name=name, mesh=mesh, space=LagrangeSpaceGeo, numbering=numbering,data=values)
    return res

def ElemFieldsToFEField(mesh):
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0)
    res = {}
    for name,values in mesh.nodeFields.items():
        if len(values.shape) == 2:
            for i in range(values.shape[1]):
                res[name+"_"+str(i)] = FEField(name=name+"_"+str(i), mesh=mesh, space=LagrangeSpaceP0, numbering=numbering, data=values[:,i])
        else:
            res[name] = FEField(name=name, mesh=mesh, space=LagrangeSpaceP0, numbering=numbering, data=values)
    return res

def FEFieldsDataToVector(listOfFields,outvec=None):
    if outvec is None:
        s = sum([f.numbering["size"] for f in listOfFields])
        outvec = np.zeros(s)
    offset = 0
    for f in listOfFields:
        outvec[offset:offset+f.numbering["size"]] = f.data
        offset += f.numbering["size"]
    return outvec

def VectorToFEFieldsData(invec,listOfFields):
    offset = 0
    for f in listOfFields:
        f.data = invec[offset:offset+f.numbering["size"]]
        offset += f.numbering["size"]

def GetPointRepresentation(listOfFEFields,fillvalue=0):

    nbfields= len(listOfFEFields)
    res = np.empty((listOfFEFields[0].mesh.GetNumberOfNodes(), nbfields) )
    for i,f in enumerate(listOfFEFields):
        res[:,i] = f.GetPointRepresentation(fillvalue=fillvalue)
    return res

def GetCellRepresentation(listOfFEFields,fillvalue=0):

    nbfields= len(listOfFEFields)
    res = np.empty((listOfFEFields[0].mesh.GetNumberOfElements(), nbfields))
    for i,f in enumerate(listOfFEFields):
        res[:,i] = f.GetCellRepresentation(fillvalue=fillvalue)
    return res

class IntegrationPointWrapper(FieldBase):
    def __init__(self,field,rule,efmask=None):
        if not isinstance(field,FEField):
            raise(Exception("IntegrationPointWrapper work only on FEFields"))

        self.feField = field
        self.rule = rule
        self._ipcache = None
        self._diffipcache = {}
        self.efmask = efmask

    @property
    def name(self):
        return self.feField.name

    def ResetCache(self):
        self._ipcache = None
        self._diffipcache = {}

    def diff(self,compName):
        from BasicTools.FE.SymWeakForm import space

        for cm in range(3):
            if space[cm] == compName:
                break
        else:
            cm = compName

        if cm not in self._diffipcache:
            from BasicTools.FE.Fields.FieldTools import TransferFEFieldToIPField
            self._diffipcache[cm] = TransferFEFieldToIPField(self.feField,der=cm,rule=self.rule, elementFilter= self.efmask)
        return self._diffipcache[cm]

    def GetIpField(self):
        if self._ipcache is None:
            from BasicTools.FE.Fields.FieldTools import TransferFEFieldToIPField
            self._ipcache =  TransferFEFieldToIPField(self.feField,der=-1,rule=self.rule, elementFilter= self.efmask)
        return self._ipcache

    def unaryOp(self,op):
        innerself = self.GetIpField()
        return innerself.unaryOp(op)

    def binaryOp(self,other,op):
        innerself = self.GetIpField()
        return innerself.binaryOp(other,op)

    @property
    def data(self):
        return self.GetIpField().data

def CreateFieldFromDescription(mesh, fieldDefinition,ftype="FE"):

    if ftype == "FE":
        from BasicTools.FE.FETools import PrepareFEComputation
        spaces,numberings,offset, NGauss = PrepareFEComputation(mesh,numberOfComponents=1)
        res = FEField(mesh=mesh,space=spaces,numbering=numberings[0])
        res.Allocate()
        FillFEField(res,fieldDefinition)

    elif ftype == "FE-P0":
        from BasicTools.FE.DofNumbering import ComputeDofNumbering
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
        numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0)
        res = FEField(mesh=mesh,space=LagrangeSpaceP0,numbering=numbering)
        res.Allocate()
        FillFEField(res,fieldDefinition)

    elif ftype == "IP":
        res = IPField(mesh=mesh,ruleName="LagrangeIsoParam")
        res.Allocate()
        FillIPField(res,fieldDefinition)

    else:
        raise(Exception("ftype not valid"))

    return res

def GetTransferOpToIPField(inField,ruleName=None,rule=None,der=-1,elementFilter=None):
    from BasicTools.FE.IntegrationsRules import GetRule
    from BasicTools.FE.Spaces.IPSpaces import GenerateSpaceForIntegrationPointInterpolation
    from BasicTools.FE.Integration import IntegrateGeneral
    from BasicTools.FE.DofNumbering import ComputeDofNumbering

    irule = GetRule(ruleName=ruleName,rule=rule)
    gaussSpace = GenerateSpaceForIntegrationPointInterpolation(irule)
    mesh = inField.mesh

    class FeToIPOp(dict):
        def __init__(self,irule,elementFilter):
            super(FeToIPOp,self).__init__()
            self.irule = irule
            self.elementFilter = elementFilter
        def dot(self,inField):
            return TransferFEFieldToIPField(inField,rule=self.irule,elementFilter=self.elementFilter,op=self)


    res = FeToIPOp(irule,elementFilter=elementFilter)

    for elemType,d in mesh.elements.items():

        eF = ElementFilter(inField.mesh,elementTypes=[elemType])

        if elementFilter is not None:
            eF = IntersectionElementFilter(mesh,(eF,elementFilter) )

        idsToTreat = eF.GetIdsToTreat(d)
        if len(idsToTreat) == 0:
            continue

        numberingRight = ComputeDofNumbering(mesh,Space=gaussSpace,elementFilter=eF)

        rightField = FEField(name="Gauss'",numbering=numberingRight,mesh=mesh,space=gaussSpace)

        from BasicTools.FE.SymWeakForm import GetField,GetTestField,space
        LF = GetField(inField.name,1)
        RF = GetTestField("Gauss",1)

        if der == -1:
            symForm = LF.T*RF
        else:
            symForm = LF.diff(space[der]).T*RF

        interpMatrixMatrix,_ = IntegrateGeneral(mesh=mesh,constants={},fields=[],wform=symForm,
                                                unkownFields= [inField],testFields=[rightField],
                                                onlyEvaluation=True,integrationRule=irule,
                                                elementFilter=eF)
        res[elemType] = interpMatrixMatrix

    return res


def TransferFEFieldToIPField(inField,ruleName=None,rule=None,der=-1,elementFilter=None,op=None):

    if op is None:
        op = GetTransferOpToIPField(inField=inField,ruleName=ruleName,rule=rule,der=der,elementFilter=elementFilter)

    from BasicTools.FE.IntegrationsRules import GetRule

    irule = GetRule(ruleName=ruleName,rule=rule)
    if elementFilter is None:
        outField = IPField(name=inField.name,mesh=inField.mesh,rule=irule)
        outField.Allocate()
    else:
        outField = RestrictedIPField(name=inField.name,mesh=inField.mesh,rule=irule,efmask=elementFilter)
        outField.Allocate()

    mesh = inField.mesh
    for elemType,d in mesh.elements.items():
        if elemType not in op :
            continue
        nbelements = op[elemType].shape[0]//len(irule[elemType][1])

        data = op[elemType].dot(inField.data)
        outField.data[elemType] = np.reshape(data,(nbelements, len(irule[elemType][1]) ),'F')
        outField.data[elemType] = np.ascontiguousarray(outField.data[elemType])

    return outField


def TranferPosToIPField(mesh,ruleName=None,rule=None,elementFilter=None):

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

def FillIPField(field,fieldDefinition):

    for f,val in fieldDefinition:
        if callable(val):
            fval = val
        else:
            fval = lambda x: val

        f.mesh = field.mesh
        if isinstance(f,ElementFilter):
            for name,elements,ids in f:
                geoSpace = LagrangeSpaceGeo[name]
                rule = field.GetRuleFor(name)
                nbip = len(rule[1])
                geospace_ipvalues = geoSpace.SetIntegrationRule(rule[0],rule[1] )

                for elid in ids:
                    for i in range(nbip):
                        valN = geospace_ipvalues.valN[i]
                        xcoor = field.mesh.nodes[elements.connectivity[elid,:],:]
                        pos = np.dot(valN ,xcoor).T
                        field.data[name][elid,i] = fval(pos)
            continue
        raise(Exception("Cant use this type of filter to fill an IPField : {}".format(str(type(f)))))

def FillFEField(field,fieldDefinition):
    for f,val in fieldDefinition:
        if callable(val):
            fval = val
        else:
            fval = lambda x: val

        f.mesh = field.mesh
        if isinstance(f,ElementFilter):

            for name,elements,ids in f:
                geoSpace = LagrangeSpaceGeo[name]
                sp = field.space[name]
                nbsf = sp.GetNumberOfShapeFunctions()
                geospace_ipvalues  = geoSpace.SetIntegrationRule(sp.posN,np.ones(nbsf) )

                for elid in ids:
                    for i in range(nbsf):
                        dofid = field.numbering[name][elid,i]
                        valN = geospace_ipvalues.valN[i]
                        xcoor = field.mesh.nodes[elements.connectivity[elid,:],:]
                        pos = np.dot(valN ,xcoor).T
                        field.data[dofid] = fval(pos)

        elif isinstance(f,NodeFilter):
            ids = f.GetIdsToTreat()
            for pid in ids:
                dofid = field.numbering.GetDofOfPoint(pid)
                pos = field.mesh.nodes[pid,:]
                field.data[dofid] = fval(pos)


def FieldsAtIp(listOfFields,rule,efmask=None):
    from BasicTools.FE.Fields.FieldTools import IntegrationPointWrapper
    res = []
    for f in listOfFields:
        if isinstance(f,FEField):
            res.append(IntegrationPointWrapper(f,rule,efmask=efmask))
        elif isinstance(f,IPField):
            if f.rule == rule:
                if efmask is None:
                    res.append(f)
                else:
                    res.append(f.GetRestrictedIPField(efmask) )
            else:
                print (f.rule)
                print (rule)
                print(f"skiping ipfield {f.name} because it has a not compatible IP rule type {str(type(f))}")
        else:
            raise(Exception("Dont know how to treat this type of field {}".format(str(type(f)) )))
    return res

class FieldsMeshTransportation():
    def __init__(self):
        self.cache_numbering = {}

    def ResetCacheData(self):
        self.cache_numbering = {}

    def GetNumbering(self, mesh, space, fromConnectivity=False,discontinuous=False):
        meshid = str(id(mesh))
        spaceid = str(id(space))
        key = (meshid,spaceid,fromConnectivity,discontinuous)

        if key not in self.cache_numbering:
            self.cache_numbering[key] = ComputeDofNumbering(mesh, space, fromConnectivity=True,discontinuous=False)
        return self.cache_numbering[key]

    def TransportFEFieldToOldMesh(self,oldmesh, infield, fillvalue=0.):
        """ function to define a FEField on the oldmesh, the infield mesh must be a
        tranformation of the oldmesh. This means the infield mesh originalids
        (for nodes and elements) must be with respect to the oldmesh

        if fillvalue is None, a partial field is generated on the old mesh.
        if fillvalue is not None, a fill over the full mesh is generated with fillvalues
        on dofs not availables on the infield
        """

        if id(infield.mesh) == id(oldmesh):
            return infield

        name = infield.name
        space = infield.space
        if infield.numbering.fromConnectivity:
            numbering = ComputeDofNumbering(oldmesh, space, fromConnectivity=True)
            res = FEField(name = name, mesh=oldmesh, space=space, numbering=numbering)
            res.Allocate(fillvalue)
            res.data[infield.mesh.originalIDNodes] = infield.data
        else :
            numbering =self.GetNumbering(oldmesh,space)
            #numbering = ComputeDofNumbering(oldmesh, space)
            res = FEField(name = name, mesh=oldmesh, space=space, numbering=numbering)
            res.Allocate(fillvalue)
            for name,data in oldmesh.elements.items():
                if name  not in infield.mesh.elements:
                    continue
                newdata = infield.mesh.elements[name]
                res.data[numbering[name][newdata.originalIds,:].flatten()]= infield.data[infield.numbering[name].flatten()]
        return res

    def TransportFEFieldToNewMesh(self,infield,newmesh):
        """ function to define a FEField on the newmesh, the new mesh must be a
        tranformation of the mesh in the infield. This means the newmesh originalids
        (for nodes and elements) must be with respect tot the mesh of the infield
        """

        if id(infield.mesh) == id(newmesh):
            return infield

        name = infield.name
        space = infield.space
        if  infield.numbering.fromConnectivity:
            numbering = ComputeDofNumbering(newmesh,space,fromConnectivity=True)
            res = FEField(name = name, mesh=newmesh, space=space, numbering=numbering)
            res.data = infield.data[newmesh.originalIDNodes]
        else:
            numbering = ComputeDofNumbering(newmesh, space)
            res = FEField(name = name, mesh=newmesh, space=space, numbering=numbering)
            res.Allocate()
            for name,data in newmesh.elements.items():
                res.data[numbering.numbering[name].flatten()] = infield.data[infield.numbering[name][data.originalIds,:].flatten()]
        return res

    def TransportIPFieldToOldMesh(self,oldmesh,ipfield):
        if id(ipfield.mesh) == id(oldmesh):
            return ipfield

        outdata = {}
        for elemType,data in ipfield.mesh.elements.items():
            indata = ipfield.data[elemType]
            outdata[elemType] = np.zeros((oldmesh.elements[elemType].GetNumberOfElements(),indata.shape[1]))
            outdata[elemType][data.originalIds,:] = ipfield.data[elemType]
        res = IPField(name=ipfield.name,mesh=oldmesh,rule=ipfield.rule,data=outdata)
        return res


    def TransportIPFieldToNewMesh(self,ipfield,newmesh):
        if id(ipfield.mesh) == id(newmesh):
            return ipfield

        outputdata = {}
        for elemType,data in newmesh.elements.items():
            outputdata[elemType] = ipfield.data[elemType][data.originalIds,:]
        res = IPField(name=ipfield.name,mesh=newmesh,rule=ipfield.rule,data=outputdata)
        return res


class FieldsEvaluator():
    def __init__(self,fields=None):
        self.originals = {}
        from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac
        self.rule = IntegrationRulesAlmanac['LagrangeIsoParam']
        self.ruleAtCenter = IntegrationRulesAlmanac['ElementCenterEval']
        self.atIp = {}
        self.atCenter = {}
        if fields is not None:
            for f in fields:
                self.Addfield(f)
        self.constants = {}
        self.efmask = None
        self.modified = True

    def AddField(self,field):
        self.originals[field.name] = field

    def AddConstant(self,name,val):
        self.constants[name] = val

    def Update(self,what="all"):
        for name,field in self.originals.items():
            if what=="all" or what =="IPField":
                self.atIp[name] = FieldsAtIp([field],self.rule,efmask=self.efmask)[0]

            if what=="all" or what =="Centroids":
                self.atCenter[name] = FieldsAtIp([field],self.ruleAtCenter,efmask=self.efmask)[0]

    def GetFieldsAt(self,on):
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
            raise Exception("Target support not supported (" + str(on) + ")" )
        result = dict(self.constants)
        result.update(res)
        return result

    def GetOptimizedFunction(self,func):
        import sympy
        import inspect
        from BasicTools.FE.SymWeakForm import space

        class OptimizedFunction():
            def __init__(self,func,constants):
                self.constants = dict(constants)
                # get arguments names
                args = inspect.getargspec(func).args
                # clean self, in the case of a member function
                if args[0] == "self":
                    args.pop(0)

                symbsimbols = {}

                self.args = args
                for n  in args:
                    if n in self.constants:
                        symbsimbols[n] = self.constants[n]
                    else:
                        symbsimbols[n] = sympy.Function(n)(*space)
                #symbolicaly evaluation of the function
                funcValues = func(**symbsimbols)

                repl, redu = sympy.cse(funcValues)

                restkeys = list(symbsimbols.keys())
                restkeys.extend(["x","y","z"])
                self.auxFuncs = []
                self.auxNames = []
                for i, v in enumerate(repl):
                    funclam = sympy.lambdify(restkeys,v[1])
                    self.auxFuncs.append(funclam)
                    funname = str(v[0])
                    self.auxNames.append(funname)
                    restkeys.append(str(v[0]))
                self.mainFunc = sympy.lambdify(restkeys,redu)

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

                for name,func in zip(self.auxNames,self.auxFuncs):
                    numericFields[name] = func(**numericFields)
                return self.mainFunc(**numericFields)[0]

        return OptimizedFunction(func,self.constants)

    def Compute(self,func, on,usesympy=False,ef=None):

        if usesympy:
            func = self.GetOptimizedFunction(func)

        fields = self.GetFieldsAt(on)
        return func(**fields)

        if on == "IPField":
            return func(**self.atIp)
        elif on == "Nodes":
            fields = {}
            for f in self.originals.values():
                if isinstance(f,FEField) and f.space == LagrangeSpaceGeo :
                    fields[f.name] = f
            return func(**fields)
        elif on == "Centroids":
            return func(**self.atCenter).data
        elif on == "FEField":
            return func(**self.originals)
        elif on == "IPField":
            return func(**self.atCenter)
        else:
            raise


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

    field = CreateFieldFromDescription(mesh, fieldDefinition, ftype="IP" )
    print(field)
    print(field.data)

    fieldDefinition.append( (NodeFilter(tags=["FirstPoint" ]), -10) )
    fieldDefinition.append( (NodeFilter(tags=["LastPoint" ]), lambda x : x[0]+1.2) )

    field = CreateFieldFromDescription(mesh, fieldDefinition )
    print(field.data)
    field.Allocate(val=0)
    FillFEField(field, fieldDefinition )


    print("--*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
    nodalTransferedField = TransferFEFieldToIPField(field,ruleName="LagrangeIsoParam",elementFilter=ElementFilter(dimensionality=1,tag="next3"))
    print("--*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--")
    print("input")
    print(field.data)
    print("output")
    print(nodalTransferedField.data )

    print(nodalTransferedField.GetIpFieldRepr(1).data)

    res = TranferPosToIPField(mesh,ruleName="LagrangeIsoParam",elementFilter=ElementFilter(dimensionality=1,tag="next3"))
    print(res[0].data)
    print(res[0].GetIpFieldRepr().data)

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


    res = FE.Compute(op,"FEField",usesympy=True)
    print(res.data)
    mesh.nodeFields["FEField_onPoints"] = GetPointRepresentation((res,)) 

    GetCellRepresentation((res,))

    print(NodeFieldToFEField(mesh))
    print(ElemFieldsToFEField(mesh))

    vect = FEFieldsDataToVector([res])

    VectorToFEFieldsData(vect,[res])
    res = FE.GetOptimizedFunction(op)

    obj = FieldsMeshTransportation()


    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True)) #pragma no cover
