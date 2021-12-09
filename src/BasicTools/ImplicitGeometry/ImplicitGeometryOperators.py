# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

import BasicTools.Helpers.ParserHelper as PH

from BasicTools.ImplicitGeometry.ImplicitGeometryFactory import RegisterClass
from BasicTools.ImplicitGeometry.ImplicitGeometryBase import ImplicitGeometryBase, smoothmin, smoothmax

class ImplicitGeometryUnion(ImplicitGeometryBase):
    """ImplicitGeometry Union of zones
        z : zones
    """
    def __init__(self,a=None):
        super(ImplicitGeometryUnion,self).__init__()
        if a is None:
            self.Zones = []
        else:
            self.Zones = a
        self.smoothControl = 0.0


    def GetDistanceToPoint(self, pos):

        op = ImplicitGeometryUnion.ApplyUnionOnLevelset

        if len(self.Zones) == 1:
            return self.ApplyInsideOut(self.Zones[0].GetDistanceToPoint(pos))
        res = op(self.Zones[0].GetDistanceToPoint(pos),self.Zones[1].GetDistanceToPoint(pos),self.smoothControl)
        for n in range(2,len(self.Zones)):
            res = op(res,self.Zones[n].GetDistanceToPoint(pos),self.smoothControl)

        return self.ApplyInsideOut(res)

    @classmethod
    def ApplyUnionOnLevelset(self,ls1,ls2,smoothControl=0):
        if smoothControl == 0:
            op = np.minimum
        else:
            def op(a,b):
                return smoothmin(a,b,smoothControl)
        return op(ls1,ls2)

    def __str__(self):
        res = "ImplicitGeometryUnion:\n"
        for z in self.Zones:
            res += "    "+ str(z) + "\n"
        return res

RegisterClass("Union",ImplicitGeometryUnion)

def CreateImplicitGeometryIntersection(ops):
    obj = ImplicitGeometryIntersection()
    obj.Zones = ops["Zones"]
    return obj

class ImplicitGeometryIntersection(ImplicitGeometryBase):
    """ImplicitGeometry Intersection of zones
        z : zones
    """

    def __init__(self,z=None):
        super(ImplicitGeometryIntersection,self).__init__()
        self.Zones = z
        self.smoothControl = 0.0

    def GetDistanceToPoint(self, pos):

        op = ImplicitGeometryIntersection().ApplyIntersectionOnLevelset

        if len(self.Zones) == 1:
            return self.ApplyInsideOut(self.Zones[0].GetDistanceToPoint(pos))

        res = op(self.Zones[0].GetDistanceToPoint(pos),self.Zones[1].GetDistanceToPoint(pos),self.smoothControl)
        for n in range(2,len(self.Zones)):
            res = ImplicitGeometryIntersection().ApplyIntersectionOnLevelset(res,self.Zones[n].GetDistanceToPoint(pos),self.smoothControl)

        return self.ApplyInsideOut(res)

    @classmethod
    def ApplyIntersectionOnLevelset(self,ls1,ls2,smoothControl=0):
        if smoothControl == 0:
            op = np.maximum
        else:
            def op(a,b):
                return smoothmax(a,b,smoothControl)
        return op(ls1,ls2)


RegisterClass("Intersection",ImplicitGeometryIntersection)


class ImplicitGeometryDifference(ImplicitGeometryBase):
    """ImplicitGeometry Difference of zones (z1-z2)
        z1 : zone
        z2 : zone
    """

    def __init__(self,Zone1=None,Zone2=None):
        super(ImplicitGeometryDifference,self).__init__()
        self.Zone1 = Zone1
        self.Zone2 = Zone2
        self.smoothControl = 0

    def GetDistanceToPoint(self, pos):
        op = ImplicitGeometryDifference.ApplyDifferenceOnLevelset
        res = op(self.Zone1.GetDistanceToPoint(pos),self.Zone2.GetDistanceToPoint(pos),self.smoothControl)
        return self.ApplyInsideOut(res)

    @classmethod
    def ApplyDifferenceOnLevelset(self,ls1,ls2,smoothControl=0):
        if smoothControl == 0:
            def op(a,b):
                return np.maximum(a,-b)
        else:
            def op(a,b):
                return smoothmax(a,-b,smoothControl)

        return op(ls1,ls2)


RegisterClass("Difference",ImplicitGeometryDifference)

def CreateImplicitGeometryOffset(ops):
    res = ImplicitGeometryOffset()
    if "Zones" in ops:
        if len(ops["Zones"]) > 1 :
            raise Exception("Cant treat more than one zone")
        res.Zone1 = ops["Zones"][0]
    PH.ReadProperties(ops, ["offset"], res)

    return res


class ImplicitGeometryOffset(ImplicitGeometryBase):
    """ImplicitGeometry offset of the levelset
        z : zone
        offset : offset value (possitive shrink, negative grow )
    """
    def __init__(self,a=None,val=0.):
        super(ImplicitGeometryOffset,self).__init__()
        self.Zone1 = a
        self.offset = val

    def GetDistanceToPoint(self, pos):
        res = self.Zone1.GetDistanceToPoint(pos) + self.offset
        return self.ApplyInsideOut(res)

RegisterClass("Offset",ImplicitGeometryOffset,CreateImplicitGeometryOffset)

def CreateImplicitGeometryInsideOut(ops):
    res = ImplicitGeometryInsideOut()
    if "Zones" in ops:
        if len(ops["Zones"]) > 1 :
            raise Exception("Cant treat more than one zone")
        res.Zone1 = ops["Zones"][0]
    return res


class ImplicitGeometryInsideOut(ImplicitGeometryBase):
    """ImplicitGeometry InsideOut  of the levelset
        z : zone
    """
    def __init__(self,a=None,val=0.):
        super(ImplicitGeometryInsideOut,self).__init__()
        self.Zone1 = a

    def GetDistanceToPoint(self, pos):
        res = -self.Zone1.GetDistanceToPoint(pos)
        return self.ApplyInsideOut(res)

RegisterClass("InsideOut",ImplicitGeometryInsideOut,CreateImplicitGeometryInsideOut)



def CreateImplicitGeometrySymmetric(ops):
    res = ImplicitGeometrySymmetric()
    if "Zones" in ops:
        if len(ops["Zones"]) > 1 :
            raise Exception("Cant treat more than one zone")
        res.Zone1 = ops["Zones"][0]
    PH.ReadProperties(ops, ["center"], res)

    return res

class ImplicitGeometrySymmetric(ImplicitGeometryBase):
    """ImplicitGeometry Central Symmetriy of a zone
        all point will be mapped to the first quadrant
        z : zone
        center : central point to compute the symetry
    """
    def __init__(self,a=None):
        super(ImplicitGeometrySymmetric,self).__init__()
        self.Zone1 = a
        self.center = np.array([0.,0.,0.])

    def GetDistanceToPoint(self, pos):
        res = self.Zone1.GetDistanceToPoint(np.abs(pos - self.center)+ self.center)
        return self.ApplyInsideOut(res)

RegisterClass("Symmetric",ImplicitGeometrySymmetric,CreateImplicitGeometrySymmetric)


def CreateImplicitGeometryShell(ops):
    res = ImplicitGeometryShell()
    if "Zones" in ops:
        if len(ops["Zones"]) > 1 :
            raise Exception("Cant treat more than one zone")
        res.Zone1 = ops["Zones"][0]
    PH.ReadProperties(ops, ["thickness"], res)

    return res


class ImplicitGeometryShell(ImplicitGeometryBase):
    """ImplicitGeometry Shell arround a zone
        z : zone
        thickness : thickness of the shell
    """
    def __init__(self,a=None,val=0.):
        super(ImplicitGeometryShell,self).__init__()
        self.Zone1 = a
        self.thickness = val

    def GetDistanceToPoint(self, pos):
        res = np.abs(self.Zone1.GetDistanceToPoint(pos) - self.thickness/2.) - self.thickness/2.
        return self.ApplyInsideOut(res)

RegisterClass("Shell",ImplicitGeometryShell,CreateImplicitGeometryShell)

def CheckIntegrity(GUI=False):
    def MustFail(func):
        try:
            func()
            raise #pragma: no cover
        except:
            pass

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([5,6,7]);
    myMesh.SetSpacing(1./(myMesh.GetDimensions()-1)*2);
    myMesh.SetOrigin([-1.,-1.,-1.]);
    myMesh.elements["hex8"].tags.CreateTag("2elems").SetIds([0,1])
    myMesh.nodesTags.CreateTag("3points").SetIds([2,3,4])
    print(myMesh)

    OnePoint3D = np.array([1,2,3])
    TwoPoints3D = np.array([[0.,0.,0.],[1.,2.,3.]], dtype=float)

    import BasicTools.TestData as TestData
    from functools import partial
    from BasicTools.Helpers.Tests import TestTempDir

    testDataPath = TestData.GetTestDataPath()

    from BasicTools.ImplicitGeometry.ImplicitGeometryObjects import ImplicitGeometrySphere

    SP1 = ImplicitGeometrySphere(center=[0,0,0],radius=1.)
    SPX = ImplicitGeometrySphere(center=[0.5,0,0],radius=1.)
    SPY = ImplicitGeometrySphere(center=[0,0.5,0],radius=1.)

    ########################### ImplicitGeometryUnion #######################
    IGUnion =ImplicitGeometryUnion([SP1,SPX,SPY])
    IGUnion(myMesh)
    IGUnion.smoothControl = 0.1
    IGUnion(myMesh)
    print(IGUnion)

    IGUnion =ImplicitGeometryUnion([SPY])
    IGUnion(myMesh)



    ########################### ImplicitGeometryIntersection #######################
    IGIntersection =ImplicitGeometryIntersection([SP1,SPX,SPY])
    IGIntersection(myMesh)
    IGIntersection.smoothControl = 0.1
    IGIntersection(myMesh)
    print(IGIntersection)

    IGIntersection = CreateImplicitGeometryIntersection({"Zones":[SP1]})
    IGIntersection(myMesh)

    ########################### ImplicitGeometryDifference #######################
    IGDifference =ImplicitGeometryDifference(Zone1=SP1,Zone2=SPX)
    IGDifference(myMesh)
    IGDifference.smoothControl = 0.1
    IGDifference(myMesh)
    print(IGDifference)
    ########################### ImplicitGeometryOffset #######################
    IGOffset = CreateImplicitGeometryOffset({"Zones":[IGUnion],"offset":0.1})
    IGOffset(myMesh)

    MustFail(partial(CreateImplicitGeometryOffset, {"Zones":[IGUnion,IGUnion],"offset":0.1}))
    ########################### ImplicitGeometrySymmetric #######################
    IGSymmetric = CreateImplicitGeometrySymmetric({"Zones":[SPX]})
    IGSymmetric(myMesh)

    MustFail(partial(CreateImplicitGeometrySymmetric, {"Zones":[IGUnion,IGUnion]}))

    ########################### ImplicitGeometryShell #######################
    IGShell = CreateImplicitGeometryShell({"Zones":[IGUnion],"thickness":0.1})
    IGShell(myMesh)

    MustFail(partial(CreateImplicitGeometryShell, {"Zones":[IGUnion,IGUnion],"thickness":0.1}))
    ########################### ImplicitGeometryInsideOut #######################
    IGIO = CreateImplicitGeometryInsideOut({"Zones":[IGUnion]})
    IGIO(myMesh)

    MustFail(partial(CreateImplicitGeometryInsideOut, {"Zones":[IGUnion,IGUnion]}))

    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=False))

