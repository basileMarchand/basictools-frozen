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
from BasicTools.FE.Fields.IPField import FieldBase


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
    res = np.empty(listOfFEFields[0].mesh.GetNumberOfElements(), nbfields)
    for i,f in enumerate(listOfFEFields):
        res[:,i] = f.GetCellRepresentation(fillvalue=fillvalue)
    return res

class IntegrationPointWrapper(FieldBase):
    def __init__(self,field,rule):
        self.feField = field
        self.rule = rule
        self._ipcache = None
        self._diffipcache = {}

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
            from BasicTools.FE.Fields.FieldTools import TransferFEFieldToIPFieldDer
            self._diffipcache[cm] = TransferFEFieldToIPFieldDer(self.feField,der=cm,rule=self.rule)
        return self._diffipcache[cm]

    def GetIpField(self):
        if self._ipcache is None:
            from BasicTools.FE.Fields.FieldTools import TransferFEFieldToIPFieldDer
            self._ipcache =  TransferFEFieldToIPFieldDer(self.feField,der=-1,rule=self.rule)
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
        outField.data[elemType] =np.reshape(data,(nbelements, len(irule[elemType][1]) ),'F')

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
                dofid = field.numbering["almanac"][("P",pid,None)]
                pos = field.mesh.nodes[pid,:]
                field.data[dofid] = fval(pos)

def FieldsAtIp(listOfFields,rule):
    from BasicTools.FE.Fields.FieldTools import IntegrationPointWrapper
    res = []
    for f in listOfFields:
        if isinstance(f,FEField):
            res.append(IntegrationPointWrapper(f,rule))
        elif isinstance(f,IPField):
            if f.rule == rule:
                res.append(f)
            else:
                print (f.rule)
                print (rule)
                print(f"skiping ipfield {f.name} because it has a not compatible IP rule type {str(type(f))}")
        else:
            raise(Exception("Dont know how to treat this type of field {}".format(str(type(f)) )))
    print([f.name for f in res])
    return res


class FieldsEvaluator():
    def __init__(self,fields=None):
        self.originals = {}
        from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac
        self.rule = IntegrationRulesAlmanac['LagrangeIsoParam']
        self.atIp = {}
        self.atCenter = {}
        if fields is not None:
            for f in fields:
                self.Addfield(f)
        self.constants = {}


    def Addfield(self,field):
        from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac
        self.originals[field.name] = field
        self.atIp.update( {f.name:f for f in FieldsAtIp([field],self.rule)})
        self.atCenter.update( {f.name:f for f in FieldsAtIp([field],IntegrationRulesAlmanac['ElementCenterEval'])})

    def GetFieldsAt(self,on):
        if on == "IntegrationPoints":
            res =  self.atIp
        elif on == "Nodes":
            res = {}
            for f in self.originals.values():
                if isinstance(f,FEField) and f.space == LagrangeSpaceGeo :
                    res[f.name] = f
        elif on == "Centroids":
            res =  self.atCenter
        elif on == "FEField":
            res = self.originals
        elif on == "IPField":
            res= self.atCenter
        else:
            raise
        result = dict(self.constants)
        result.update(res)
        return result

    def Compute(self,func, on,usesympy=False,ef=None):

        if usesympy :
            import sympy
            import inspect
            spec = inspect.getargspec(func)
            print("-------------spec")
            print(spec)
            print(spec.args)
            args = spec.args
            if args[0] == "self":
                args.pop(0)
            print(args)
            symbsimbols = {v:sympy.symbols(v) for v  in args}
            from BasicTools.FE.SymWeakForm import space

            symbsimbols = {v:sympy.Function(v)(*space) for v  in args}

            funcValues = func(**symbsimbols)

            print("-----------------------------------------------------------")
            print(funcValues )
            repl,redu = sympy.cse(funcValues)
            print(repl)
            print("----REDU-------------------------------------------------------")
            print(redu)
            fields = self.GetFieldsAt(on)
            fields = {sympy.Function(k)(*space):v for k,v in fields.items()}
            restkeys = []
            restvalues = []
            for i, v in enumerate(repl):
                print("-----")
                print(v)
                print("with subms", fields)
                res = v[1].subs(fields).doit()
#                sympy.lambdify(restkeys,redu)(fields)
                print("subs resulr " , res)
                fields[v[0]] = res
                restkeys.append( str(v[0]) )
                restvalues.append( res)

                for i,v in enumerate(restvalues):
                    print( i)
                    print( v)
                raise

            print("----REDU-------------------------------------------------------")
            if type(redu):
                redu = [x.subs(fields).doit() for x in redu]
            else:
                redu = redu.subs(fields).doti()
            return redu
            funclam = sympy.lambdify(restkeys,redu)
            print(redu)
            print(restkeys)
            print(restvalues)
            return funclam(*restvalues)



        fields = self.GetFieldsAt(on)
        return func(**fields)

        if on == "IntegrationPoints":
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
            return func(**self.atCenter).data
        else:
            raise


def CheckIntegrity(GUI=False):

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars
    import BasicTools.Containers.ElementNames as EN

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
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True)) #pragma no cover
