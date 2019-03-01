# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.FE.Integration import IntegrateGeneral
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.DofNumbering import ComputeDofNumbering

def GetElementaryMatrixForFormulation(elemName,wform,unknownNames,space = LagrangeSpaceP1):

    from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
    mesh = UnstructuredMesh()

    mesh.nodes = np.asarray(space[elemName].posN,dtype=float)
    mesh.originalIDNodes = np.arange(0,mesh.GetNumberOfNodes(),dtype=np.int)

    elements = mesh.GetElementsOfType(elemName)
    elements.connectivity = np.arange(space[elemName].GetNumberOfShapeFunctions(),dtype=np.int)

    elements.connectivity.shape = (1,space[elemName].GetNumberOfShapeFunctions())
    elements.GetTag("3D").AddToTag(0)

    elements.originalIds = np.arange(0,1,dtype=np.int)
    elements.cpt = elements.connectivity.shape[0]

    mesh.PrepareForOutput()
    print(mesh)
    numbering = ComputeDofNumbering(mesh,space,)

    unkownFields = []
    for name in unknownNames:
        print(name)
        unkownFields.append(FEField(name,mesh,space,numbering))


    M,f = IntegrateGeneral(mesh=mesh,wform=wform, unkownFields= unkownFields,constants={},fields=[])

    return M


def CheckIntegrity(GUI=False):
    from BasicTools.FE.SymPhysics import MecaPhysics

    import BasicTools.Containers.ElementNames as EN


    mecaPhysics = MecaPhysics()
    wform = mecaPhysics.GetBulkFormulation(1.0,0.3)

    res = GetElementaryMatrixForFormulation(EN.Hexaedron_8,wform, unknownNames =mecaPhysics.GetPrimalNames() )

    for line in res.toarray().tolist():
        print(line)

    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity(GUI=True))

    print("Done")
