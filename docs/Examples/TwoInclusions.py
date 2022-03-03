
from BasicTools.FE.UnstructuredFeaSym import UnstructuredFeaSym

#Main class to
problem = UnstructuredFeaSym()

# the mecanical problem
from BasicTools.FE.SymPhysics import MecaPhysics
mecaPhysics = MecaPhysics()

# Definition of the degree of the spaces [1 or 2]
mecaPhysics.SetSpaceToLagrange(P=1)

# add weak form terms to the tanget matrix
mecaPhysics.AddBFormulation( "Bulk",mecaPhysics.GetBulkFormulation(1.0,0.3)  )
mecaPhysics.AddBFormulation( "Inclusion1",mecaPhysics.GetBulkFormulation(5.0,0.3)  )
youngModulusInclusionII = 0.5
mecaPhysics.AddBFormulation( "Inclusion2",mecaPhysics.GetBulkFormulation(0.5,0.3)  )

# add weak form term to the rhs
mecaPhysics.AddLFormulation( "Right", mecaPhysics.GetForceFormulation([1,0,0],1)  )

problem.physics.append(mecaPhysics)

# the boundary conditions
from BasicTools.FE.KR.KRBlock import KRBlockVector
dirichlet = KRBlockVector()
dirichlet.AddArg("u").On('Left').Fix0().Fix1().Fix2().To()

problem.solver.constraints.AddConstraint(dirichlet)

# Read The mesh
from BasicTools.IO.GmshReader import ReadGmsh
mesh = ReadGmsh("TwoInclussions.msh")
mesh.ConvertDataForNativeTreatment()
print(mesh)

problem.SetMesh(mesh)

# we compute the numbering (maping from the mesh to the linear system)
problem.ComputeDofNumbering()

from BasicTools.Helpers.Timer import Timer
with Timer("Assembly "):
    k,f = problem.GetLinearProblem()

#compute the constraints to add to the system
problem.ComputeConstraintsEquations()


with Timer("Solve"):
    problem.Solve(k,f)

problem.PushSolutionVectorToUnkownFields()

from BasicTools.FE.Fields.FieldTools import GetPointRepresentation
# Recover a point representation of the displacement
problem.mesh.nodeFields["sol"] = GetPointRepresentation(problem.unkownFields)

from BasicTools.FE.Fields.FEField import FEField
#Creation of a fake fields to export the rhs member
rhsFields = [ FEField(mesh=mesh,space=None,numbering=problem.numberings[i]) for i in range(3) ]
from BasicTools.FE.Fields.FieldTools import VectorToFEFieldsData
VectorToFEFieldsData(f,rhsFields)
problem.mesh.nodeFields["RHS"] = GetPointRepresentation(rhsFields)

print("Done solve")
print("Compute of the strain energy only on the second inclusion (integral in each element) ")

import BasicTools.FE.SymWeakForm as SWF
from BasicTools.FE.MaterialHelp import HookeIso
from BasicTools.Containers.Filters import ElementFilter

symdep = SWF.GetField("u",3)
K = HookeIso(youngModulusInclusionII,0.3,dim=3)
symCellData = SWF.GetField("cellData",1)
symCellDataT = SWF.GetTestField("cellData",1)

print("Post process")

EnerForm = SWF.ToVoigtEpsilon(SWF.Strain(symdep)).T*K*SWF.ToVoigtEpsilon(SWF.Strain(symdep))*symCellDataT

print("Post process Eval")
ff = ElementFilter(mesh=problem.mesh, tag="Inclusion2")
from BasicTools.FE.DofNumbering import  ComputeDofNumbering
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP0
from BasicTools.FE.Integration import IntegrateGeneral
p0Numbering = ComputeDofNumbering(mesh,LagrangeSpaceP0)

energyDensityField = FEField(name="cellData",mesh=problem.mesh,numbering=p0Numbering,space=LagrangeSpaceP0)

m,energyDensity = IntegrateGeneral(mesh=problem.mesh, wform=EnerForm,  constants={},
                        fields=problem.unkownFields, unkownFields = [ energyDensityField ], elementFilter=ff)
print("energyDensity",energyDensity)
energyDensityField.data = energyDensity

problem.mesh.elemFields["Energy"] = energyDensity

import numpy as np
print("Strain energy on the second inclusion:", np.sum(energyDensity) )
from BasicTools.IO import XdmfWriter as XW
writer = XW.XdmfWriter('TwoInclussions_Output.xdmf')
writer.SetHdf5(False)
writer.Open()
writer.Write(mesh,PointFields=list(mesh.nodeFields.values()), PointFieldsNames=list(mesh.nodeFields.keys()),
                CellFields=list(mesh.elemFields.values()), CellFieldsNames=list(mesh.elemFields.keys()))
writer.Close()


