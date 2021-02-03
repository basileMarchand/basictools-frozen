# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

from BasicTools.Containers import Filters

__cache__ = {}
__cacheSize__ = 1
numberingAlgorithm = "NumpyBase"
#numberingAlgorithm = "DictBase"

def GetNumberingFromCache(mesh,space,elementFilter=None,discontinuous=False,fromConnectivity=False):
    return None
    key = (id(mesh),id(space),id(elementFilter),discontinuous,fromConnectivity)
    return __cache__.get(key,None)

def SetNumberingToCache(obj, mesh,space,elementFilter=None,discontinuous=False,fromConnectivity=False):
    if len(__cache__) >=  __cacheSize__:
        for k in __cache__:
            del __cache__[k]
            break
    key = (id(mesh),id(space),id(elementFilter),discontinuous,fromConnectivity)
    __cache__[key] = obj

def ComputeDofNumbering(mesh,Space,dofs=None,fromConnectivity=False,elementFilter=None,discontinuous=False):

    if dofs is None:
        cachedData = GetNumberingFromCache(mesh=mesh, space=Space, fromConnectivity=fromConnectivity, elementFilter=elementFilter, discontinuous=discontinuous )
    if cachedData is not None:
        return cachedData

    if numberingAlgorithm == "NumpyBase":
        from BasicTools.FE.Numberings.DofNumberingNumpy import DofNumberingNumpy
        res = DofNumberingNumpy()
    elif numberingAlgorithm == "DictBase":
        from BasicTools.FE.Numberings.DofNumberingDict import DofNumberingDict
        res = DofNumberingDict()
    else:
        raise(Exception(f"Numbering algorithm of type {numberingAlgorithm} not available "))


    if fromConnectivity:
        if dofs is not None or elementFilter is not None or discontinuous:
           raise(Exception("cant take dofs, sign, discontinuous or elementFilter different from the default values"))
        res.ComputeNumberingFromConnectivity(mesh,Space)
        SetNumberingToCache(res,mesh=mesh, space=Space, fromConnectivity=fromConnectivity, elementFilter=elementFilter, discontinuous=discontinuous )
        return res
    else:
        if dofs is not None:
            res = dofs

        res.ComputeNumberingGeneral(mesh=mesh, space=Space, elementFilter=elementFilter, discontinuous=discontinuous )
        SetNumberingToCache(res,mesh=mesh, space=Space, fromConnectivity=fromConnectivity, elementFilter=elementFilter, discontinuous=discontinuous )

        return res

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    res2 = CreateCube([2.,2.,2.],[-1.0,-1.0,-1.0],[2./46, 2./46,2./46])

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo, LagrangeSpaceP0, LagrangeSpaceP1, LagrangeSpaceP2, ConstantSpaceGlobal

    spacesToTest = { "ConstantSpaceGlobal":ConstantSpaceGlobal,"LagrangeSpaceGeo":LagrangeSpaceGeo, "LagrangeSpaceP0":LagrangeSpaceP0, "LagrangeSpaceP1":LagrangeSpaceP1, "LagrangeSpaceP2":LagrangeSpaceP2 }


    from BasicTools.FE.Spaces.IPSpaces import GenerateSpaceForIntegrationPointInterpolation
    from BasicTools.FE.IntegrationsRules import GetRule
    irule = GetRule(ruleName="LagrangeP1")
    gaussSpace = GenerateSpaceForIntegrationPointInterpolation(irule)
    spacesToTest["gaussSpace"] =gaussSpace
    import time
    for spacename, space in spacesToTest.items():
        print("******************************************************************************* {} **********************************************************".format(spacename))
        # on a tag
        print("----------------------{} tag -----------------------------".format(spacename))
        st = time.time()
        numbering = ComputeDofNumbering(res2,space,elementFilter= Filters.ElementFilter(mesh=res2,tag="X0") )
        print("----------------------{} tag -----------------------------".format(spacename))
        print(time.time()-st)

        # all
        print("----------------------{} all -----------------------------".format(spacename))
        st = time.time()
        numbering = ComputeDofNumbering(res2,space)
        print("----------------------{} all -----------------------------".format(spacename))
        print(time.time()-st)


        # from connectivity
        print("----------------------{} from connectivity -----------------------------".format(spacename))
        st = time.time()
        numbering = ComputeDofNumbering(res2, space,fromConnectivity=True)
        print("----------------------{} from connectivity -----------------------------".format(spacename))
        print(time.time()-st)

        # 3D using filter
        print("----------------------{} 3D filter-----------------------------".format(spacename))
        st = time.time()
        numbering = ComputeDofNumbering(res2,space,elementFilter=Filters.ElementFilter(mesh=res2,dimensionality=3) )
        print("----------------------{} 3D filter-----------------------------".format(spacename))
        print(time.time()-st)
    return "ok"

def CheckIntegrityUsingAlgo(algo,GUI=False):

    import BasicTools.FE.DofNumbering  as DN
    tmpAlgo = DN.numberingAlgorithm
    DN.numberingAlgorithm = algo
    try:
        res = DN.CheckIntegrity(GUI)
        DN.numberingAlgorithm = tmpAlgo
        return res
    except:
        DN.numberingAlgorithm = tmpAlgo
        raise


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
