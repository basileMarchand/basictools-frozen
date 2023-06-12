# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Union, Dict, Optional, List

import numpy as np
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh

from BasicTools.NumpyDefs import ArrayLike, PBasicFloatType
from BasicTools.FE.KR.KRBase import KRBase
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.Linalg.ConstraintsHolder import ConstraintsHolder

class KRFromDistribution(KRBase):
    """Class to impose values in the primary unknown of the system
    for the moment this only work with iso-parametric fields
    """
    def __init__(self):
        super().__init__()
        self.type = "RFromVector",
        self.fieldsData = []

    def FixUsingNodeFields(self, mesh: UnstructuredMesh, nodalField:Dict[str,ArrayLike] ) -> KRFromDistribution :
        """Set the fields to be used to impose the dof values.
        Internally we convert the numpy array into fields with the function :
        BasicTools.FE.Fields.FieldTools.NodeFieldToFEField()

        Parameters
        ----------
        mesh : UnstructuredMesh
            the mesh
        nodalField : Dict[str,ArrayLike]
            a dict with key:value in the form of "u":numpy(3,NbNodes)

        Returns
        -------
        KRFromDistribution
            self
        """

        from BasicTools.FE.Fields.FieldTools import NodeFieldToFEField
        dataFields = NodeFieldToFEField(mesh, nodalField)
        return self.FixUsingFEFields(dataFields)

    def FixUsingFEFields(self, inputFEFields:List[FEField]) -> KRFromDistribution:
        """ Set the fields to be used to impose the dof values

        Parameters
        ----------
        inputFEFields : List[FEField]
            a list of FEField (the order is not important)

        Returns
        -------
        KRFromDistribution
            self
        """
        self.fieldsData = inputFEFields
        return self

    def GenerateEquations(self, mesh: UnstructuredMesh, fields:List[FEField], CH:Optional[ConstraintsHolder]=None, solutionFields:Optional[List[FEField]]=None, reactionFields:Optional[List[FEField]]=None) -> ConstraintsHolder :
        """Generate the equation to impose the constraints defined by this class.
        More information about this

        Parameters
        ----------
        mesh : UnstructuredMesh
            The mesh
        fields : List[FEField]
            The list of fields to be used to compute the numbering  (in the order on the linear system to be solved)
        CH : Optional[ConstraintsHolder], optional
            the ConstraintsHolder to store the equations, if None is provided a new object is created, by default None
        solutionFields : Optional[List[FEField]], optional
            In the case of an iterative solver (i.e. in a newton loop) the current solution, by default None
        reactionFields : Optional[List[FEField]], optional
            The user can provide the reactions fields ( the integral of the internal forces), by default None

        Returns
        -------
        ConstraintsHolder
            the ConstraintsHolder with the new equations

        Raises
        ------
        RuntimeError
            if no data is provided. the user must call FixUsingFEFields or FixUsingNodeField to provide the fields
        LookupError
            If the user impose a restriction on a field (using .AddArg()) but not available on the problem
        """

        if len(self.fieldsData) == 0:
            raise Exception("Please set data to use (FixUsingFEFields,FixUsingNodeField) ")

        offsets, fieldOffsets, totalNumberOfDofs  = self._ComputeOffsets(fields)

        CH = self._GetConstraintHolder(CH)
        if CH.nbdof == 0:
            CH.nbdof = totalNumberOfDofs


        unknownFieldDic = {f.name:f for f in fields }

        ef = ElementFilter(mesh=mesh, tags=self.on )

        dataFieldDic = {f.name:f for f in self.fieldsData.values() }

        class FieldData():
            def __init__(self, name:str, unknownField, dataField, offset, solutionFields=None, reactionFields=None  ):
                self.name = name
                self.unknownField = unknownField[name]
                self.dataField = dataField[name]
                self.offset = offset[name]

                if solutionFields is None:
                    self.solutionField = None
                else:
                    solutionFields = {f.name:f for f in solutionFields }
                    self.solutionField = solutionFields[name[1:]]


                if reactionFields is None:
                    reactionFields= None
                else:
                    reactionFields = {f.name:name for f in reactionFields }
                    self.reactionField = reactionFields["R"+name[1:]]


        fieldsToTreat = []
        # sanity check and prep of dofs to impose
        for arg in self.args:
            if "_" in arg:
                if arg in unknownFieldDic.keys() and arg in dataFieldDic.keys() :
                    fieldsToTreat.append( FieldData(arg, unknownFieldDic, dataFieldDic, fieldOffsets, solutionFields,reactionFields  ) )
                else:
                    raise LookupError(f"Field ({arg}) not found.")

            else:
                if arg in unknownFieldDic.keys() and arg in dataFieldDic.keys() :
                    fieldsToTreat.append( FieldData(arg, unknownFieldDic, dataFieldDic, fieldOffsets, solutionFields,reactionFields  ) )
                else:
                    for i in range(3):
                        composedName = f"{arg}_{i}"
                        if composedName in fieldOffsets:
                            fieldsToTreat.append( FieldData(composedName, unknownFieldDic, dataFieldDic, fieldOffsets, solutionFields,reactionFields  ) )
                        else:
                            break

        for fieldToTreat in fieldsToTreat:
            totalNumberOfShapeFuntions = fieldToTreat.unknownField.numbering.size
            treated = np.zeros(totalNumberOfShapeFuntions, dtype=bool)
            offset = fieldToTreat.offset
            for name, data, elementIds in ef:
                numberOfShapeFuntions = fieldToTreat.unknownField.space[name].GetNumberOfShapeFunctions()
                for elemId in elementIds:
                    for i in range(numberOfShapeFuntions):
                        dofIid = fieldToTreat.unknownField.numbering[name][elemId,i]
                        if treated[dofIid]:
                            continue
                        treated[dofIid] = True
                        dataValue = fieldToTreat.dataField.data[fieldToTreat.dataField.numbering[name][elemId,i]]
                        if fieldToTreat.solutionField is not None:
                            solValue = fieldToTreat.solutionField.data[fieldToTreat.solutionField.numbering[name][elemId,i]]
                        else:
                            solValue = 0.

                        CH.AddFactor(dofIid+offset,1)
                        CH.AddConstant(dataValue-solValue)
                        CH.NextEquation()
        return CH

    def __str__(self) -> str:
        res = self.type + " "
        if len(self.args) > 1:
            res += "("
        res += " and ".join(str(x) for x in self.args)
        if len(self.args) > 1:
            res += ")"

        return res

def CheckIntegrity(GUI=False):

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    mesh = CreateSquare()
    print(mesh)

    imposedU = mesh.nodes + 3

    obj = KRFromDistribution()
    obj.AddArg("u").On("X0")
    obj.FixUsingNodeFields(mesh, {"u":imposedU})

    from BasicTools.FE.Fields.FieldTools import NodeFieldToFEField
    uFields = list(NodeFieldToFEField(mesh, {"u":imposedU}).values())

    ch = obj.GenerateEquations(mesh,uFields)
    mat, dofs = ch.ToSparse()
    print(mesh.nodes)
    print(imposedU.T.ravel()[dofs[:-1]])
    print(mat.toarray())
    print(dofs)

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    nx = 5; ny = 5; nz = 5;
    mesh = CreateCube(dimensions=[nx,ny,nz], origin=[0,0,0.], spacing=[1/(nx-1),1/(ny-1), 1/(nz-1)], ofTetras=False )

    from BasicTools.FE.UnstructuredFeaSym import UnstructuredFeaSym
    problem = UnstructuredFeaSym()

    from BasicTools.FE.SymPhysics import MecaPhysics
    mecaPhysics = MecaPhysics()
    mecaPhysics.SetSpaceToLagrange(1)
    mecaPhysics.AddBFormulation( "3D",mecaPhysics.GetBulkFormulation(1.0,0.3)  )
    problem.physics.append(mecaPhysics)

    from BasicTools.FE.KR.KRBlock import KRBlockVector
    dirichlet = KRBlockVector()
    dirichlet.AddArg("u").On('Z0').Fix0().Fix1().Fix2()
    problem.solver.constraints.AddConstraint(dirichlet)

    obj = KRFromDistribution()
    # this to impose only displacement in the z direction
    #obj.AddArg("u_2").On("Z1")
    obj.AddArg("u").On("Z1")
    imposedU = np.full_like(mesh.nodes,0.1)
    imposedU[:,2] += np.cos(3*mesh.nodes[:,0])*0.1
    obj.FixUsingNodeFields(mesh, {"u":imposedU})
    problem.solver.constraints.AddConstraint(obj)

    #problem.solver.constraints.SetConstraintsMethod("Penalisation")
    problem.solver.SetAlgo("Direct")

    problem.SetMesh(mesh)

    problem.ComputeDofNumbering()

    k,f = problem.GetLinearProblem()
    problem.ComputeConstraintsEquations()
    problem.Solve(k,f)

    problem.PushSolutionVectorToUnkownFields()
    from BasicTools.FE.Fields.FieldTools import GetPointRepresentation
    problem.mesh.nodeFields["sol"] = GetPointRepresentation(problem.unkownFields)

    from BasicTools.Actions.OpenInParaView import OpenInParaView
    OpenInParaView(mesh,filename="Test_KRFromVector.xmf",run=GUI)
    from BasicTools.Containers.Filters import NodeFilter
    ids = NodeFilter(mesh=mesh,etags=["Z1"]).GetIdsToTreat()
    error = np.abs(problem.mesh.nodeFields["sol"][ids,:] - imposedU[ids,:])

    if np.linalg.norm(error/imposedU[ids,:]) > 1e-8:
        RuntimeError("Error in the CheckIntegrity ")

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
