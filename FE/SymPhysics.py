# -*- coding: utf-8 -*-


from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from BasicTools.FE.WeakForm import SymWeakToNumWeak
from BasicTools.FE.WeakForm import Gradient,Divergence
from BasicTools.FE.WeakForm import GetField,GetTestField
from sympy import Symbol

class Physics(BOO):
    def __init__(self):
        self.integrationRule = None
        self.spaces = [None]
        self.bilinearWeakFormulations = []
        self.linearWeakFormulations = []
        self.numberings=[]
        self.spaceDimension = 3

    def ExpandNames(self,data):
        if data[1] == 1:
            return data[0]
        return [ data[0] + "_" +str(d)  for d in range(data[1]) ]

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

    def AddBFormulation(self, zone, data ) :
        self.bilinearWeakFormulations.append((zone,data) )

    def AddLFormulation(self, zone, data ) :
        self.linearWeakFormulations.append((zone,data) )


    def GetNumberOfUnkownFields(self):
        return len(self.GetPrimalNames())

    def ComputeDofNumbering(self,mesh,tagsToKeep=None):
        from BasicTools.FE.DofNumbering import ComputeDofNumbering
        self.numberings = []
        for d in range(self.GetNumberOfUnkownFields()):
            self.numberings.append(None)
            if tagsToKeep is not None:
                for tag in tagsToKeep:
                    self.numberings[d] = ComputeDofNumbering(mesh,self.spaces[d],fromConnectivity =False,tag=tag,dofs=self.numberings[d])
            for tag,form in self.linearWeakFormulations:
                self.numberings[d] = ComputeDofNumbering(mesh,self.spaces[d],fromConnectivity =False,tag=tag,dofs=self.numberings[d])
            for tag,form in self.bilinearWeakFormulations:
                self.numberings[d] = ComputeDofNumbering(mesh,self.spaces[d],fromConnectivity =False,tag=tag,dofs=self.numberings[d])
            #print("size of numbering", len(self.numberings))
            #print("(mesh.elements['quad4'].connectivity[0,:]", mesh.elements["quad4"].connectivity[0,:])
            #return

class MecaPhysics(Physics):
    def __init__(self):
        super(MecaPhysics,self).__init__()
        self.mecaPrimalName = ("u",3)
        self.pressureName = ("p",1)
        self.mecaSpace = None

        self.young = 1.
        self.poisson = 0.3


    def SetMecaPrimalName(self,name):
        self.mecaPrimalName = name

    def GetPrimalNames(self):
        return self.ExpandNames(self.mecaPrimalName)

    def GetBulkFormulation(self,young=None, poisson=None,factor=None ):
        from BasicTools.FE.WeakForm import GetMecaElasticProblem
        from BasicTools.FE.MaterialHelp import HookeIso

        if young is None:
            young = self.young
        if poisson is None:
            poisson = self.poisson
        if factor is not None:
            young *=factor

        self.HookeLocalOperator = HookeIso(young,poisson,dim=self.mecaPrimalName[1])
        Symwfb = GetMecaElasticProblem(self.mecaPrimalName[0],dim=self.mecaPrimalName[1],K=self.HookeLocalOperator)
        return Symwfb

    def GetPressureFormulation(self,pressureName):
        from BasicTools.FE.WeakForm import GetMecaNormalPressure
        if pressureName is None:
            pressureName = self.pressureName
        Symwfp = GetMecaNormalPressure(self.pressureName,name=self.mecaPrimalName)
        wfp = SymWeakToNumWeak(Symwfp)
        return wfp

    def GetForceFormulation(self,direction,flux="f"):
        from BasicTools.FE.WeakForm import GetTestField
        from sympy import Symbol

        ut = GetTestField(self.mecaPrimalName[0],self.mecaPrimalName[1])
        if isinstance(flux,str):
            f = Symbol(flux)
        else:
            f = float(flux)

        from sympy.matrices import Matrix
        if not isinstance(direction,Matrix):
            direction = Matrix([direction]).T

        wflux = f*direction.T*ut
        print(wflux)
        wfp = SymWeakToNumWeak(wflux)
        return wfp

    def PostTraitementFormulations(self):
        import BasicTools.FE.WeakForm as wf
        symdep = GetField("u",3)
        nodalEnergy = GetField("ElasticEnergy",1)
        nodalEnergyT = GetTestField("ElasticEnergy",1)
        symEner = wf.ToVoigtEpsilon(wf.Strain(symdep)).T*self.HookeLocalOperator*wf.ToVoigtEpsilon(wf.Strain(symdep))*nodalEnergyT
        symMass = nodalEnergy.T*nodalEnergyT
        post1 = ("ElasticEnergy","nodal",symEner + symMass)
        post2 = ("ElasticEnergy","global",symEner )
        #post3 = ("ElasticEnergy","globalPreCalculated",3.14159 )

        return (post1,post2)




class ThermalPhysics(Physics):
    def __init__(self):
        super(ThermalPhysics,self).__init__()
        self.thermalPrimalName = ("t",1)
        self.thermalSpace = None

    def GetPrimalNames(self):
        return self.thermalPrimalName

    def SetThermalPrimalName(self,name):
        self.thermalPrimalName = name

    def GetBulkFormulation(self,alpha=1 ):
        from BasicTools.FE.WeakForm import Gradient
        from BasicTools.FE.WeakForm import GetField,GetTestField
        #from sympy import Identity

        t  = GetField(self.thermalPrimalName[0],1)
        tt = GetTestField(self.thermalPrimalName[0],1)
        Symwfb = Gradient(t,self.spaceDimension).T*(alpha)*Gradient(tt,self.spaceDimension)
        #wfb = SymWeakToNumWeak(Symwfb)
        return Symwfb

    def GetNormalFlux(self,flux="f"):
        from BasicTools.FE.WeakForm import GetTestField
        from sympy import Symbol

        tt = GetTestField(self.thermalPrimalName,1)
        if isinstance(flux,str):
            f = Symbol(flux)
        else:
            f = float(flux)

        wflux = f*tt
        wfp = SymWeakToNumWeak(wflux)
        return wfp

class StokesPhysics(Physics):
    def __init__(self):
        super(StokesPhysics,self).__init__()
        self.velocityPrimalName = "v"
        self.pressurePrimalName = "p"

    def GetPrimalNames(self):
        res = [self.velocityPrimalName + "_" + str(c) for c in [0,1,2] ]
        res.append(self.pressurePrimalName)
        return res

    def SetSpaceToLagrange(self,P=None,isoParam=None):
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP2

        self.spaces = [LagrangeSpaceP2]*self.spaceDimension
        self.spaces.append(LagrangeSpaceP1)
        self.integrationRule =  "LagrangeP2"


    def GetBulkFormulation(self,mu=1.):
        from BasicTools.FE.WeakForm import Gradient,Divergence
        from BasicTools.FE.WeakForm import GetField,GetTestField
        #from sympy import Identity

        v  = GetField(self.velocityPrimalName,self.spaceDimension)
        vt = GetTestField(self.velocityPrimalName,self.spaceDimension)
        p  = GetField(self.pressurePrimalName,1)
        pt = GetTestField(self.pressurePrimalName,1)

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
        res = self.mecaPhys.GetBulkFormulation(young=young, poisson=poisson)
        res += self.thermalPhys.GetBulkFormulation(alpha=alpha)

        # need to add the clouplig terms

        #res += self.HookeLocalOperator

        return res


def CheckIntegrity(GUI=False):
    res = ThermoMecaPhysics()
    print(res.GetBulkFormulation())
    print(res.GetPrimalNames())

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))# pragma: no cover
