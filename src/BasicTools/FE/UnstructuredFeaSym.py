# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


from BasicTools.FE.FeaBase import FeaBase
from BasicTools.Containers.Filters import ElementFilter

from BasicTools.FE.WeakForms.NumericalWeakForm import SymWeakToNumWeak
import BasicTools.FE.SymWeakForm as SWF

from BasicTools.FE.Integration import IntegrateGeneral
from BasicTools.FE.Fields.FEField import FEField


class UnstructuredFeaSym(FeaBase):
    def __init__(self,spaceDim=3 ):
        super(UnstructuredFeaSym,self).__init__(spaceDim=spaceDim)

        self.constants = {}
        self.fields = {}

        self.physics = []

        self.spaces = []
        self.numberings = []

    def SetMesh(self,support):
        super(UnstructuredFeaSym,self).SetMesh(support)
        for phy in self.physics:
            phy.Reset()

    def ComputeDofNumbering(self,tagsToKeep=None,fromConnectivity=False):

        self.spaces = []
        self.numberings = []

        for phy in self.physics:
            self.spaces.extend(phy.spaces)
            phy.ComputeDofNumbering(self.mesh,tagsToKeep,fromConnectivity=fromConnectivity)
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
                field.mesh = self.mesh
                field.space = phy.spaces[dim]
                field.Allocate()
                self.unkownFields.append(field)

        self.solver.constraints.SetNumberOfDofs(self.totalNumberOfDof)

    def GetLinearProblem(self,computeK=True, computeF=True, lform = None ,bform = None,unkownFields=None):

        rhsRes = None
        lhsRes = None

        if lform is not None:
            if unkownFields is None:
               unkownFields =  self.unkownFields
            for ff,form in lform:
                if form is None:
                    continue
                self.PrintDebug("integration of lform "+ str(ff) )
                #self.PrintDebug("integration of lform "+ str(form) )
                _,f = IntegrateGeneral(mesh=self.mesh,
                                       wform=form,
                                       constants=self.constants,
                                       fields=list(self.fields.values()),
                                       unkownFields= unkownFields,
                                       elementFilter=ff)
                if rhsRes is None:
                    rhsRes = f
                else:
                    rhsRes += f
            return (lhsRes,rhsRes)

        if bform is not None:
            if unkownFields is None:
               unkownFields =  self.unkownFields
            for ff,form in bform:
                if form is None:
                    continue
                self.PrintDebug("integration of bform " + str(ff) )
                #self.PrintDebug("integration of bform " + str(form) )
                k,f = IntegrateGeneral(mesh=self.mesh,
                                       wform=form,
                                       constants=self.constants,
                                       fields=list(self.fields.values()),
                                       unkownFields= unkownFields,
                                       elementFilter=ff)
                if rhsRes is None:
                    rhsRes = f
                else:
                    rhsRes += f

                if lhsRes is None:
                    lhsRes = k
                else:
                    lhsRes += k

            return (lhsRes,rhsRes)

        if (computeF):
          self.PrintDebug("In Integration F")
          for phy in self.physics:
              linearWeakFormulations = phy.linearWeakFormulations

              for ff,form in linearWeakFormulations:
                if form is None:
                    continue
                self.PrintDebug("integration of f "+ str(ff) )
                _,f = IntegrateGeneral(mesh=self.mesh,wform=form,  constants=self.constants, fields=list(self.fields.values()),unkownFields= self.unkownFields,
                                integrationRuleName=phy.integrationRule,elementFilter=ff)
                if rhsRes is None:
                    rhsRes = f
                else:
                    rhsRes += f

        if (computeK):
          self.PrintDebug("In Integration K")
          for phy in self.physics:
              bilinearWeakFormulations = phy.bilinearWeakFormulations

              for ff,form in bilinearWeakFormulations:
                if form is None:
                    continue
                self.PrintDebug("Integration of bilinear formulation on : " + str(ff))
                k,f = IntegrateGeneral(mesh=self.mesh,wform=form,  constants=self.constants, fields=list(self.fields.values()), unkownFields= self.unkownFields,
                                integrationRuleName=phy.integrationRule,elementFilter=ff)
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
           if res.lower()!="ok": return "ERROR: "+res + " " + str(P) + " " + str(tetra)
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
    from BasicTools.FE.KR.KRBlock import KRBlockVector
    dirichlet = KRBlockVector()
    dirichlet.AddArg("u").On('Z0').Fix0().Fix1().Fix2().To(offset=[1,2,3],first=[1,0,1] )
    dirichlet.constraintDiretions= "Global"

    problem.solver.constraints.AddConstraint(dirichlet)

    # the mesh

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    nx = 11; ny = 12; nz = 13;
    nx = nx//2
    ny = ny//2
    nz = nz//2

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
    problem.Print("k.shape " + str(k.shape) )
    problem.Print("f.shape "+ str(f.shape))

    with Timer("Solve"):
        problem.Solve(k,f)

    problem.PushSolutionVectorToUnkownFields()
    from BasicTools.FE.Fields.FieldTools import GetPointRepresentation
    problem.mesh.nodeFields["sol"] = GetPointRepresentation(problem.unkownFields)


    #Creation of a fake fields to export the rhs member
    rhsFields = [ FEField(mesh=mesh,space=None,numbering=problem.numberings[i]) for i in range(3) ]
    from BasicTools.FE.Fields.FieldTools import VectorToFEFieldsData
    VectorToFEFieldsData(f,rhsFields)
    problem.mesh.nodeFields["RHS"] = GetPointRepresentation(rhsFields)

    print("done solve")

    symdep = SWF.GetField("u",3)
    from BasicTools.FE.MaterialHelp import HookeIso
    K = HookeIso(1,0.3,dim=3)
    symCellData = SWF.GetField("cellData",1)
    symCellDataT = SWF.GetTestField("cellData",1)

    print("Post process")

    EnerForm = SWF.ToVoigtEpsilon(SWF.Strain(symdep)).T*K*SWF.ToVoigtEpsilon(SWF.Strain(symdep))*symCellDataT + symCellData.T*symCellDataT

    print("Post process Eval")
    ff = ElementFilter(mesh=problem.mesh, tag="3D")
    energyDensityField = FEField(name="cellData",mesh=problem.mesh,numbering=problem.numberings[0],space=problem.spaces[0])
    m,energyDensity = IntegrateGeneral(mesh=problem.mesh, wform=EnerForm,  constants={},
                        fields=problem.unkownFields, unkownFields = [ energyDensityField ], integrationRuleName="NodalEvalP"+str(P),
                        onlyEvaluation=True,
                        elementFilter=ff)
    print("energyDensity",energyDensity)
    energyDensity /= m.diagonal()
    energyDensityField.data = energyDensity

    problem.mesh.nodeFields["PEnergy"] = energyDensityField.GetPointRepresentation()
    problem.mesh.elemFields["PEnergy_fromP"+str(P) ] = energyDensityField.GetCellRepresentation()


    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
    from BasicTools.FE.DofNumbering import ComputeDofNumbering

    P0Numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0,elementFilter=ElementFilter(mesh=mesh,dimensionality=mesh.GetDimensionality()))

    P0energyDensityField = FEField(name="cellData",mesh=problem.mesh,numbering=P0Numbering,space=LagrangeSpaceP0)
    m,P0energyDensity = IntegrateGeneral(mesh=problem.mesh, wform=EnerForm, constants={},
                        fields=problem.unkownFields, unkownFields = [ P0energyDensityField ],
                        integrationRuleName="ElementCenterEval",
                        onlyEvaluation=True,
                        elementFilter=ff)

    P0energyDensity /= m.diagonal()
    P0energyDensityField.data = P0energyDensity

    problem.mesh.elemFields["CEnergy"] = P0energyDensityField.GetCellRepresentation()

    if GUI  :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
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
