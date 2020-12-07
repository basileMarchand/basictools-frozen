# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import math

from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from BasicTools.FE.SymWeakForm import Gradient,Divergence, GetField,GetTestField, GetScalarField


class Physics(BOO):
    def __init__(self):
        self.integrationRule = None
        self.spaces = [None]
        self.bilinearWeakFormulations = []
        self.linearWeakFormulations = []
        self.numberings = None
        self.spaceDimension = 3

    def Reset(self):
        self.numberings = None

    def ExpandNames(self,data):
        if data[1] == 1:
            return data[0]
        return [ data[0] + "_" +str(d)  for d in range(data[1]) ]

    def GetBulkMassFormulation(self,alpha=1.):
        u = self.primalUnknown
        ut = self.primalTest

        a = GetScalarField(alpha)

        ener = u.T*ut*a
        return ener

    def SetSpaceToLagrange(self,P=None,isoParam=None):
        if P is None and isoParam is None:
            raise(ValueError("Please set the type of integration P=1,2 or isoParam=True"))

        if P is not None and isoParam is not None:
            raise(ValueError("Please set the type of integration P=1,2 or isoParam=True"))

        if isoParam is None or isoParam == False :
            if P == 1 :
                from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
                space = LagrangeSpaceP1
                self.integrationRule =  "LagrangeP1"
            elif  P == 2:
                from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP2
                space = LagrangeSpaceP2
                self.integrationRule =  "LagrangeP2"
            else:
                raise(ValueError("I dont understand"))
        else:
            from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
            space = LagrangeSpaceGeo
            self.integrationRule =  None

        self.spaces = [space]*len(self.GetPrimalNames())

    def AddBFormulation(self, zoneOrElementFilter, data ) :

        if type(zoneOrElementFilter) == str:
            ef = ElementFilter(tag=zoneOrElementFilter)
        else:
            ef  = zoneOrElementFilter

        self.bilinearWeakFormulations.append((ef,data) )

    def AddLFormulation(self, zoneOrElementFilter, data ) :


        if type(zoneOrElementFilter) == str:
            ef = ElementFilter(tag=zoneOrElementFilter)
        else:
            ef  = zoneOrElementFilter

        self.linearWeakFormulations.append((ef,data) )

    def GetNumberOfUnkownFields(self):
        return len(self.GetPrimalNames())

    def ComputeDofNumbering(self,mesh,tagsToKeep=None,fromConnectivity=False):
        from BasicTools.FE.DofNumbering import ComputeDofNumbering
        if self.numberings is None:
            self.numberings = [None]*self.GetNumberOfUnkownFields()
        else:
            return

        from BasicTools.Containers.Filters import ElementFilter,UnionElementFilter
        allFilters = UnionElementFilter(mesh)

        ff = ElementFilter(mesh,tags=tagsToKeep)
        allFilters.filters.append(ff)
        allFilters.filters.extend([f for f, form in self.linearWeakFormulations] )
        allFilters.filters.extend([f for f, form in self.bilinearWeakFormulations] )

        for d in range(self.GetNumberOfUnkownFields()):
            for i in range(0,d):
                if self.spaces[d] == self.spaces[i]:
                    self.numberings[d] = self.numberings[i]
                    break
            else:
                if fromConnectivity:
                    self.numberings[d] = ComputeDofNumbering(mesh,self.spaces[d],fromConnectivity = True ,dofs=self.numberings[d])
                else:
                    self.numberings[d] = ComputeDofNumbering(mesh,self.spaces[d],fromConnectivity =False,elementFilter=allFilters,dofs=self.numberings[d])

    def ComputeDofNumberingFromConnectivity(self,mesh):
        self.ComputeDofNumbering(mesh, fromConnectivity=True)


class MecaPhysics(Physics):
    def __init__(self,dim=3):
        super(MecaPhysics,self).__init__()
        self.dim = dim

        self.mecaPrimalName = ("u",self.dim)
        self.SetMecaPrimalName(self.mecaPrimalName[0])

        self.mecaSpace = None

        self.young = 1.
        self.poisson = 0.3
        self.density = 1.

        self.planeStress = True

    def SetMecaPrimalName(self,name):
        self.mecaPrimalName = (name,self.dim)
        self.primalUnknown = GetField(self.mecaPrimalName[0],self.mecaPrimalName[1])
        self.primalTest = GetTestField(self.mecaPrimalName[0],self.mecaPrimalName[1])

    def GetPrimalNames(self):
        return self.ExpandNames(self.mecaPrimalName)


    def GetBulkFormulation(self,young=None, poisson=None,alpha=None ):
        from BasicTools.FE.SymWeakForm import ToVoigtEpsilon,Strain
        from BasicTools.FE.MaterialHelp import HookeLaw

        u = self.primalUnknown
        ut = self.primalTest

        if young is None:
            young = self.young
        if poisson is None:
            poisson = self.poisson

        young = GetScalarField(young)*GetScalarField(alpha)

        op = HookeLaw()
        op.Read({"E":young, "nu":poisson})
        self.HookeLocalOperator = op.HookeIso(dim=self.dim,planeStress=self.planeStress)
        Symwfb = ToVoigtEpsilon(Strain(u,self.dim)).T*self.HookeLocalOperator*ToVoigtEpsilon(Strain(ut,self.dim))
        return Symwfb

    def GetPressureFormulation(self,pressure):
        ut = self.primalTest

        p = GetScalarField(pressure)

        from BasicTools.FE.SymWeakForm import GetNormal
        Normal = GetNormal(self.dim)

        Symwfp = p*Normal.T*ut

        return Symwfp

    def GetForceFormulation(self,direction,flux="f"):

        ut = self.primalTest
        f = GetScalarField(flux)

        from sympy.matrices import Matrix
        if not isinstance(direction,Matrix):
            direction = Matrix([direction]).T

        return  f*direction.T*ut

    def GetAccelerationFormulation(self,direction,density=None):

        if density is None:
            density = self.density

        ut = self.primalTest
        density = GetScalarField(density)
        from sympy.matrices import Matrix
        if not isinstance(direction,Matrix):
            direction = [GetScalarField(d) for d in direction]
            direction = Matrix([direction]).T

        return  density*direction.T*ut



    def PostTraitementFormulations(self):
        import BasicTools.FE.SymWeakForm as wf
        symdep = self.primalUnknown


        nodalEnergyT = GetTestField("elastic_energy",1)
        symEner = 0.5*wf.ToVoigtEpsilon(wf.Strain(symdep)).T*self.HookeLocalOperator*wf.ToVoigtEpsilon(wf.Strain(symdep))*nodalEnergyT


        trStrainT = GetTestField("tr_strain_",1)
        symTrStrain = wf.Trace(wf.Strain(symdep))*trStrainT

        trStressT = GetTestField("tr_stress_",1)
        symTrStress = wf.Trace(wf.FromVoigtSigma(wf.ToVoigtEpsilon(wf.Strain(symdep)).T*self.HookeLocalOperator))*trStressT

        postQuantities = {"elastic_energy" : symEner,
                          "tr_strain_": symTrStrain,
                          "tr_stress_": symTrStress }

        return postQuantities


class MecaPhysicsAxi(MecaPhysics):
    def __init__(self):
        super(MecaPhysics,self).__init__(dim=2)

    def GetFieldR(self):
        return GetScalarField("r")

    def GetBulkFormulation(self,young=None, poisson=None,alpha=None ):
        from BasicTools.FE.MaterialHelp import HookeLaw

        u = self.primalUnknown
        ut = self.primalTest

        if young is None:
            young = self.young
        if poisson is None:
            poisson = self.poisson

        young = GetScalarField(young)*GetScalarField(alpha)

        op = HookeLaw()
        op.Read({"E":young, "nu":poisson})
        self.HookeLocalOperator = op.HookeIso(dim=self.dim,planeStress=self.planeStress, axisymetric=self.axisymetric)

        r = self.GetFieldR()

        from BasicTools.FE.SymWeakform import StrainAxyCol
        epsilon_u = StrainAxyCol(u,r)
        epsilon_ut = StrainAxyCol(ut,r)
        Symwfb = 2*math.pi*epsilon_u*self.HookeLocalOperator*epsilon_ut*r
        return Symwfb

    def GetPressureFormulation(self,pressure):
        return super().GetPressureFormulation(pressure)*self.GetFieldR()

    def GetForceFormulation(self,direction,flux="f"):
        return super().GetForceFormulation(direction,flux)*self.GetFieldR()

    def GetAccelerationFormulation(self,direction,density=None):
        return super().GetForceFormulation(direction,density)*self.GetFieldR()

    def PostTraitementFormulations(self):
        import BasicTools.FE.SymWeakForm as wf
        symdep = self.primalUnknown

        pir2 = 2*math.pi*self.GetFieldR()

        nodalEnergyT = GetTestField("strain_energy",1)
        symEner = pir2*0.5*wf.ToVoigtEpsilon(wf.Strain(symdep)).T*self.HookeLocalOperator*wf.ToVoigtEpsilon(wf.Strain(symdep))*nodalEnergyT

        from sympy import prod

        trStrainT = GetTestField("tr(strain)",1)
        strain = wf.StrainAxyCol(symdep)
        symTrStrain = prod(strain[0:3]).T*trStrainT

        trStressT = GetTestField("tr(stress)",1)
        stress = strain*self.HookeLocalOperator
        symTrStress = prod(stress[0:3]).T*trStressT

        postQuantities = {"strain_energy" : symEner,
                          "tr(strain)": symTrStrain,
                          "tr(stress)": symTrStress }

        return postQuantities

class BasicPhysics(Physics):
    def __init__(self):
        super(BasicPhysics,self).__init__()
        self.PrimalNameTrial = ("u",1)
        self.PrimalNameTest = ("u",1)
        self.Space = None
        self.SetPrimalName(self.PrimalNameTrial[0])

    def SetPrimalName(self,unknowName,testName=None,unknowDim=1,testDim=1):
        self.PrimalNameTrial = (unknowName,unknowDim)
        if testName is None:
            testName = unknowName
        self.PrimalNameTest = (testName,testDim)
        self.primalUnknown = GetField(self.PrimalNameTrial[0],self.PrimalNameTrial[1])
        self.primalTest = GetTestField(self.PrimalNameTest[0],self.PrimalNameTest[1])


    def GetPrimalNames(self):
        return [self.PrimalNameTrial[0]]

    def GetPrimalDims(self):
        return [self.PrimalNameTrial[1]]


    def GetBulkFormulation_dudi_dtdj(self,u=0,t=0,i=0,j=0,alpha=1.):

        a = GetScalarField(alpha)

        unk = self.primalUnknown

        if self.PrimalNameTrial[1] > 1:
            dtestdj = Gradient(unk,self.spaceDimension)[i,u]
        else:
            dtestdj = Gradient(unk,self.spaceDimension)[i]

        ut = self.primalTest
        if self.PrimalNameTest[1] > 1:
            dtrialdi = Gradient(ut,self.spaceDimension)[j,t]
        else:
            dtrialdi = Gradient(ut,self.spaceDimension)[j]

        Symwfb = dtrialdi*(a)*dtestdj
        return Symwfb

    def GetBulkLaplacian(self,alpha=1):
        from BasicTools.FE.SymWeakForm import Gradient
        a = GetScalarField(alpha)
        u = self.primalUnknown
        ut = self.primalTest
        Symwfb = Gradient(u,self.spaceDimension).T*(a)*Gradient(ut,self.spaceDimension)
        return Symwfb

    def GetFlux(self,flux="f"):
        tt = self.primalTest
        f = GetScalarField(flux)

        return f*tt

class ThermalPhysics(Physics):
    def __init__(self):
        super(ThermalPhysics,self).__init__()
        self.thermalPrimalName = ("t",1)
        self.SetPrimalNames(self.thermalPrimalName)
        self.thermalSpace = None

    def GetPrimalNames(self):
        return [ self.thermalPrimalName[0]]

    def SetPrimalNames(self,data):
        self.thermalPrilamName = data
        self.primalUnknown = GetField(self.thermalPrimalName[0],1)
        self.primalTest = GetTestField(self.thermalPrimalName[0],1)

    def SetThermalPrimalName(self,name):
        self.thermalPrimalName = name

    def GetBulkFormulation(self, alpha=1. ):
        t = self.primalUnknown
        tt = self.primalTest

        if hasattr(alpha, '__iter__'):
            from sympy import diag
            K = diag(*alpha)
            Symwfb = Gradient(t,self.spaceDimension).T*K*Gradient(tt,self.spaceDimension)
        else:
            alpha = GetScalarField(alpha)
            Symwfb = Gradient(t,self.spaceDimension).T*(alpha)*Gradient(tt,self.spaceDimension)

        return Symwfb

    def GetNormalFlux(self,flux="f"):

        tt = self.primalTest
        f = GetScalarField(flux)

        return f*tt

class StokesPhysics(Physics):
    def __init__(self):
        super(StokesPhysics,self).__init__()
        self.velocityPrimalName = ("v",3)
        self.pressurePrimalName = ("p",1)
        self.SetPrimalNames(self.velocityPrimalName[0])


    def GetPrimalNames(self):
        res = [self.velocityPrimalName[0] + "_" + str(c) for c in range(self.velocityPrimalName[1]) ]
        res.append(self.pressurePrimalName)
        return res

    def SetMecaPrimalName(self,vname,pname):
        self.velocityPrimalName = (vname,self.dim)
        self.pressurePrimalName = (pname,1)
        self.primalUnknownV = GetField(self.velocityPrimalName[0],self.velocityPrimalName[1])
        self.primalUnknownP = GetField(self.pressurePrimalName[0],self.pressurePrimalName[1])
        self.primalTestV = GetTestField(self.velocityPrimalName[0],self.velocityPrimalName[1])
        self.primalTestP = GetTestField(self.pressurePrimalName[0],self.pressurePrimalName[1])


    def SetSpaceToLagrange(self,P=None,isoParam=None):
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP2

        self.spaces = [LagrangeSpaceP2]*self.spaceDimension
        self.spaces.append(LagrangeSpaceP1)
        self.integrationRule =  "LagrangeP2"


    def GetBulkFormulation(self,mu=1.):

        mu = GetScalarField(mu)
        v  = self.primalUnknownV
        vt = self.primalTestV
        p  = self.primalUnknownP
        pt = self.primalTestP

        res = Gradient(v,self.spaceDimension).T*mu*Gradient(vt,self.spaceDimension)  -  Divergence(vt,self.spaceDimension)*p + pt*Divergence(v,self.spaceDimension)

        return res



class ThermoMecaPhysics(Physics):
    def __init__(self):
        super(ThermoMecaPhysics,self).__init__()
        self.mecaPhys = MecaPhysics()
        self.thermalPhys =ThermalPhysics()

    def GetPrimalNames(self):
        res = self.mecaPhys.GetPrimalNames()
        res.extend(self.thermalPhys.GetPrimalNames())
        return res

    def GetBulkFormulation(self,young=1., poisson=0.3, alpha=1.):
        res = self.mecaPhys.GetBulkFormulation(young=young, poisson=poisson,alpha=alpha)
        res += self.thermalPhys.GetBulkFormulation(alpha=alpha)

        # need to add the clouplig terms

        #res += self.HookeLocalOperator

        return res


def CheckIntegrity(GUI=False):
    res = ThermoMecaPhysics()
    print(res.GetBulkFormulation())
    print(res.GetPrimalNames())

    print(BasicPhysics().GetBulkFormulation_dudi_dtdj())
    t = BasicPhysics()
    t.spaceDimension = 3
    t.SetPrimalName("U","V",3,3)
    print(t.primalUnknown)
    print(t.primalTest)
    print(t.GetBulkFormulation_dudi_dtdj(u=0,i=1,t=1,j=2) )
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))# pragma: no cover
