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


def ComputeMassMatrix(mesh):

    from scipy.sparse import coo_matrix
    from BasicTools.FE.IntegrationsRules import Lagrange as Lagrange
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
    from BasicTools.Containers import Filters

    nbNodes = mesh.GetNumberOfNodes()
    dim     = mesh.GetDimensionality()
    
    spaces = LagrangeSpaceGeo
    for name, data in mesh.elements.items():
        p,w =  Lagrange(name)
        spaces[name].SetIntegrationRule(p,w)
      
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)
    numberings = [numbering]*dim
    
    offset = []
    totaldofs = 0
    for n in numberings:
        offset.append(totaldofs)
        totaldofs += n["size"]    
      
    ev = []
    ei = []
    ej = []

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dim)
    
    for name,data,ids in ff:
        p,w =  Lagrange(name)
        lenNumbering = len(numberings[0][name][0,:])
        ones = np.ones(lenNumbering,dtype=int)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]

            for ip in range(len(w)):
                Jack, Jdet, Jinv = spaces[name].GetJackAndDetI(ip,xcoor)
                for j in range(dim):
                    leftNumbering = numberings[j][name][el,:] + offset[j]
                    left = spaces[name].valN[ip]
                    ev.extend(((w[ip]*Jdet)*np.outer(left, left)).ravel())
                    for i in leftNumbering:
                        ei.extend(i*ones)
                        ej.extend(leftNumbering.ravel())

    return coo_matrix((ev, (ei,ej)), shape=(dim*nbNodes,dim*nbNodes)).tocsr()



def CheckIntegrity(GUI=False):
    from BasicTools.FE.SymPhysics import MecaPhysics

    import BasicTools.Containers.ElementNames as EN


    mecaPhysics = MecaPhysics()
    wform = mecaPhysics.GetBulkFormulation(1.0,0.3)

    res = GetElementaryMatrixForFormulation(EN.Hexaedron_8,wform, unknownNames =mecaPhysics.GetPrimalNames() )

    for line in res.toarray().tolist():
        print(line)
        
    import BasicTools.TestData as BasicToolsTestData
    from BasicTools.IO import GeofReader as GR
    mesh = GR.ReadGeof(BasicToolsTestData.GetTestDataPath()+"cube2.geof")
    ComputeMassMatrix(mesh)
    
    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity(GUI=True))

    print("Done")
