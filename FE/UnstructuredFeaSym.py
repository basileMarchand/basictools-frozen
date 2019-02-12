# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.FE.FeaBase import FeaBase

from BasicTools.FE.WeakForm import SymWeakToNumWeak
import BasicTools.FE.WeakForm as wf

from BasicTools.FE.Integration import IntegrateGeneral
from BasicTools.FE.Integration import Integrate
from BasicTools.FE.Fields.FEField import FEField


class UnstructuredFeaSym(FeaBase):
    def __init__(self,spaceDim=3 ):
        super(UnstructuredFeaSym,self).__init__(spaceDim=spaceDim)

        self.constants = {}
        self.fields = {}

        self.physics = []

        self.spaces = []
        self.numberings = []

    def ComputeDofNumbering(self,tagsToKeep=None):

        self.spaces = []
        self.numberings = []

        for phy in self.physics:
            self.spaces.extend(phy.spaces)
            phy.ComputeDofNumbering(self.mesh,tagsToKeep)
            self.numberings.extend(phy.numberings)

        self.totalNumberOfDof = 0
        for num in self.numberings:
            self.totalNumberOfDof += num["size"]

        self.unkownFields = []
        for phy in self.physics:
            for dim in range(phy.GetNumberOfUnkownFields()):
                field = FEField()
                field.numbering = phy.numberings[dim]
                field.name = phy.GetPrimalNames()[dim]
                field.data = None
                field.mesh = self.mesh
                field.space = phy.spaces[dim]
                self.unkownFields.append(field)

        self.solver.constraints.SetNumberOfDofs(self.totalNumberOfDof)

    def GetLinearProblem(self,computeK=True, computeF=True, linearWeakFormulations = None):

        rhsRes = None
        lhsRes = None

        if linearWeakFormulations is not None:
            for zone,form in linearWeakFormulations:
                if form is None:
                    continue
                self.PrintDebug("integration of f "+ str(zone) )
                _,f = IntegrateGeneral(mesh=self.mesh,wform=form, tag=zone, constants=self.constants, fields=list(self.fields.values()),unkownFields= self.unkownFields)
                if rhsRes is None:
                    rhsRes = f
                else:
                    rhsRes += f
            return (lhsRes,rhsRes)


        if (computeF):
          self.PrintDebug("In Integration F")
          for phy in self.physics:
              linearWeakFormulations = phy.linearWeakFormulations

              for zone,form in linearWeakFormulations:
                if form is None:
                    continue
                self.PrintDebug("integration of f "+ str(zone) )
                _,f = IntegrateGeneral(mesh=self.mesh,wform=form, tag=zone, constants=self.constants, fields=list(self.fields.values()),unkownFields= self.unkownFields,
                                integrationRuleName=phy.integrationRule)
                if rhsRes is None:
                    rhsRes = f
                else:
                    rhsRes += f

        if (computeK):
          self.PrintDebug("In Integration K")
          for phy in self.physics:
              bilinearWeakFormulations = phy.bilinearWeakFormulations

              for zone,form in bilinearWeakFormulations:
                if form is None:
                    continue
                self.PrintDebug("Integration of bilinear formulation on : " + str(zone))
                k,f = IntegrateGeneral(mesh=self.mesh,wform=form, tag=zone, constants=self.constants, fields=list(self.fields.values()), unkownFields= self.unkownFields,
                                integrationRuleName=phy.integrationRule)
                if not (f is None):
                    if rhsRes is None:
                        rhsRes = f
                    else:
                        rhsRes += f

                if lhsRes is None:
                    lhsRes = k
                else:
                    lhsRes += k

        return (lhsRes,rhsRes)

def CheckIntegrity(GUI=False):
    for P in [1,2]:
        for tetra in [False,True]:
           print("in CheckIntegrityFlexion P="+str(P)+" tetra="+str(tetra))
           res = CheckIntegrityFlexion( P = P,tetra = tetra,GUI=GUI);
           if res.lower()!="ok": return res + " " + str(P) + " " + str(tetra)
    return "ok"

def CheckIntegrityFlexion(P,tetra,GUI=False):

    # the main class
    problem = UnstructuredFeaSym()
    problem.SetGlobalDebugMode(True)

    # the mecanical problem
    from BasicTools.FE.SymPhysics import MecaPhysics
    mecaPhysics = MecaPhysics()
    mecaPhysics.SetSpaceToLagrange(P=P)
    mecaPhysics.AddBFormulation( "3D",mecaPhysics.GetBulkFormulation(1.0,0.3)  )
    mecaPhysics.AddLFormulation( "Z1", mecaPhysics.GetForceFormulation([1,0,0],0.002)  )
    mecaPhysics.AddLFormulation( "Z0", None  )
    problem.physics.append(mecaPhysics)

    # the boundary conditions
    from BasicTools.FE.KR.KRBlock import KRBlock
    dirichlet = KRBlock()
    dirichlet.AddArg("u").On('Z0').Fix0().Fix1().Fix2().To(offset=[1,2,3],first=[1,0,1] )
    dirichlet.constraintDiretions= "Global"

    problem.solver.constraints.AddConstraint(dirichlet)

    # the mesh

    from BasicTools.Containers.UnstructuredMeshTools import CreateCube

    nx = 11; ny = 12; nz = 13;
    mesh = CreateCube(dimensions=[nx,ny,nz],origin=[0,0,0.], spacing=[1./(nx-1),1./(ny-1), 10./(nz-1)], ofTetras=tetra )
    problem.SetMesh(mesh)
    print(mesh)
    # we compute the numbering
    problem.ComputeDofNumbering()
    #print(mecaPhysics.numberings[0])


    from BasicTools.Helpers.Timer import Timer
    with Timer("Assembly "):
        k,f = problem.GetLinearProblem()


    #problem.solver = LinearProblem()
    #problem.solver.SetAlgo("EigenCG")
    #problem.solver.SetAlgo("EigenLU")

    problem.solver.SetAlgo("Direct")
    problem.ComputeConstraintsEquations()
    print("k.shape", k.shape)
    print("f.shape",f.shape)

    with Timer("Solve"):
        problem.Solve(k,f)

    problem.PushVectorToUnkownFields()

    print("done solve")

    symdep = wf.GetField("u",3)
    from BasicTools.FE.MaterialHelp import HookeIso
    K = HookeIso(1,0.3,dim=3)
    symCellData = wf.GetField("cellData",1)
    symCellDataT = wf.GetTestField("cellData",1)

    print("Post process")

    EnerForm = wf.ToVoigtEpsilon(wf.Strain(symdep)).T*K*wf.ToVoigtEpsilon(wf.Strain(symdep))*symCellDataT + symCellData.T*symCellDataT

    symEner = wf.ToVoigtEpsilon(wf.Strain(symdep))[0,0]*symCellDataT + symCellData.T*symCellDataT

    #from BasicTools.Actions.OpenInParaView import OpenInParaView
    #problem.PushVectorToMesh(True,f,"normalFlux")
    #problem.PushVectorToMesh(True,problem.sol,"sol")
    #OpenInParaView(mesh)
    #return()
    print("Post process Eval")

    m,energyDensity = Integrate(mesh=problem.mesh, wform=EnerForm, tag="3D", constants={},
                        fields={f.name:f for f in problem.unkownFields}, dofs=["cellData"], spaces=[problem.spaces[0] ],
                        numbering=[problem.numberings[0]], integrationRuleName="NodalEvalP"+str(P),
                        onlyEvaluation=True)
    print("energyDensity",energyDensity)
    energyDensity /= m.diagonal()

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    P0Numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0,tag="3D")

    m,P0energyDensity = Integrate(mesh=problem.mesh, wform=EnerForm, tag="3D", constants={},
                        fields={f.name:f for f in problem.unkownFields}, dofs=["cellData"], spaces=[LagrangeSpaceP0 ],
                        numbering=[P0Numbering], integrationRuleName="ElementEvalGeo",
                        onlyEvaluation=True)

    P0energyDensity /= m.diagonal()


    if GUI  :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        problem.PushVectorToMesh(True,f,"normalFlux",problem.numberings)
        problem.PushVectorToMesh(True,problem.sol,"sol",problem.numberings)
        problem.PushVectorToMesh(True,energyDensity,"PEnergy",[problem.numberings[0]])
        problem.PushVectorToMesh(False,energyDensity,"CEnergy_FromP",[problem.numberings[0]])
        problem.PushVectorToMesh(False,P0energyDensity,"CEnergy",[P0Numbering])

        OpenInParaView(mesh,filename="UnstructuredFeaSym_Sols_P"+str(P)+("Tetra"if tetra else "Hexa")+".xmf",run=True)

    print(Timer())
    return("ok")



if __name__ == '__main__':

    import time
    starttime = time.time()
    print(CheckIntegrity(True))#pragma: no cover

    stoptime = time.time()
    print("Total Time {0}s".format(stoptime-starttime))
    print("Done")
