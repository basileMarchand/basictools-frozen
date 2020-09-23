# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Containers.Filters import ElementFilter,NodeFilter
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.Fields.IPField import IPField
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo


def CreateFieldFromDescription(mesh, fieldDefinition,ftype="FE"):

    if ftype == "FE":
        from BasicTools.FE.FETools import PrepareFEComputation
        spaces,numberings,offset, NGauss = PrepareFEComputation(mesh,numberOfComponents=1)

        res = FEField(mesh=mesh,space=spaces,numbering=numberings[0])
        res.Allocate()

        FillFEField(res,fieldDefinition)
    elif ftype == "IP":

        res = IPField(mesh=mesh,ruleName="LagrangeIsoParam")
        res.Allocate()
        FillIPField(res,fieldDefinition)
    else:
        raise()
    return res

def TransferFEFieldToIPField(inField,ruleName=None,rule=None):

    outField = IPField(name=inField.name,mesh=inField.mesh,ruleName=ruleName,rule=rule)
    outField.Allocate()

    for name,elements in outField.mesh.elements.items():
        space = inField.space[name]
        rule = outField.GetRuleFor(name)
        nbip = len(rule[1])
        space.SetIntegrationRule(rule[0],rule[1] )
        for elid in range(elements.GetNumberOfElements()):
            for i in range(nbip):
                valN = space.valN[i]
                valField = inField.data[inField.numbering[name][elid,:]]
                val  = np.dot(valN ,valField).T
                outField.data[name][elid,i] = val

    return outField

def TransferFEFieldToIPFieldDer(inField,der=-1,ruleName=None,rule=None,elementFilter=None):

    if elementFilter is None:
        elementFilter = ElementFilter(mesh=inField.mesh)

    outField = IPField(name=inField.name,mesh=inField.mesh,ruleName=ruleName,rule=rule)
    outField.Allocate()

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
    geospaces = LagrangeSpaceGeo

    mesh = inField.mesh
    for name, data, elids in elementFilter:
        NnodeperEl = EN.numberOfNodes[name]

        geospace = geospaces[name]
        space = inField.space[name]
        rule = outField.GetRuleFor(name)
        nbip = len(rule[1])

        geospace.SetIntegrationRule(rule[0],rule[1] )
        space.SetIntegrationRule(rule[0],rule[1] )
        numbering = inField.numbering[name]
        outputHeavyData = outField.data[name]
        for elid in elids:
            xcoor = np.array([mesh.nodes[data.connectivity[elid,j]] for j in range(NnodeperEl)])

            for iip in range(nbip):
                if der == -1:
                    valN = space.valN[iip]
                else:
                    _Jack, _Jdet, Jinv = geospace.GetJackAndDetI(iip,xcoor)
                    valN = Jinv(space.valdphidxi[iip])[der]

                dofsids = numbering[elid]
                valField = inField.data[dofsids]
                val  = np.dot(valN, valField).T
                outputHeavyData[elid,iip] = val

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
                geoSpace.SetIntegrationRule(rule[0],rule[1] )

                for elid in ids:
                    for i in range(nbip):
                        valN = geoSpace.valN[i]
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

                geoSpace.SetIntegrationRule(sp.posN,np.ones(nbsf) )

                for elid in ids:
                    for i in range(nbsf):
                        dofid = field.numbering[name][elid,i]
                        valN = geoSpace.valN[i]
                        xcoor = field.mesh.nodes[elements.connectivity[elid,:],:]
                        pos = np.dot(valN ,xcoor).T
                        field.data[dofid] = fval(pos)

        elif isinstance(f,NodeFilter):
            ids = f.GetIdsToTreat()
            for pid in ids:
                dofid = field.numbering["almanac"][("P",pid,None)]
                pos = field.mesh.nodes[pid,:]
                field.data[dofid] = fval(pos)

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

    print(field.data)

    nodalTransferedField = TransferFEFieldToIPField(field,ruleName="LagrangeIsoParam")
    print(nodalTransferedField.data )

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True)) #pragma no cover
