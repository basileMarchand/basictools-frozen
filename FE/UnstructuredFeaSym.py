# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.FE.FeaBase import FeaBase

class FeaSymMeca(FeaBase):
    def __init__(self,dim=3 ):
        super(FeaSymMeca,self).__init__(dim=dim)

        #tuple  (zone,formulation)
        self.bilinearWeakFormulations = []
        #tuple  (zone,formulation)
        self.linearWeakFormulations = []
        self.constants = {}
        self.fields = {}

    def SetSpaceToLagrange(self):
        from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
        self.space = LagrangeSpaceGeo

    def SetBulkFormulationToMeca(self,young=1,poisson=0.3,zone="3D"):
        from BasicTools.FE.WeakForm import GetMecaElasticProblem
        from BasicTools.FE.MaterialHelp import HookeIso
        K = HookeIso(1,0.3,dim=3)
        Symwfb = GetMecaElasticProblem("u",3,K)
        from BasicTools.FE.WeakForm import SymWeakToNumWeak
        wfb = SymWeakToNumWeak(Symwfb)
        self.bilinearWeakFormulations = [(zone,wfb)]

    def AddPressure(self,zone,pressureName="p"):
        from BasicTools.FE.WeakForm import GetMecaNormalPressure
        Symwfp = GetMecaNormalPressure(pressureName)
        from BasicTools.FE.WeakForm import SymWeakToNumWeak
        wfp = SymWeakToNumWeak(Symwfp)
        self.linearWeakFormulations.append((zone,wfp))

    def AddForce(self,zone,direction,flux="F"):
        from BasicTools.FE.WeakForm import GetTestField
        from sympy import Symbol

        ut = GetTestField("u",3)
        if isinstance(flux,str):
            f = Symbol(flux)
        else:
            f = float(flux)

        from sympy.matrices import Matrix
        if not isinstance(direction,Matrix):
            direction = Matrix([direction]).T

        wflux = f*direction.T*ut
        from BasicTools.FE.WeakForm import SymWeakToNumWeak
        wfp = SymWeakToNumWeak(wflux)
        self.linearWeakFormulations.append((zone,wfp))

    def GetLinearProblem(self):

        dim = 3
        dofs= ["u_"+ str(i) for i in range(dim)]
        spaces = [self.space]*dim
        numberings = [self.numbering]*dim
        from BasicTools.FE.Integration import Integrate

        self.PrintDebug("In Integration F")
        rhsRes = None
        lhsRes = None
        for zone,form in self.linearWeakFormulations:
            f = Integrate(mesh=self.mesh,wform=form, tag=zone, constants=self.constants, fields=self.fields, dofs=dofs,spaces=spaces,numbering=numberings)[1]
            if rhsRes is None:
                rhsRes = f
            else:
                rhsRes += f

        self.PrintDebug("In Integration K")
        for zone,form in self.bilinearWeakFormulations:
            self.PrintDebug("Integration of bilinear formulation on : " + str(zone))
            k,f = Integrate(mesh=self.mesh,wform=form, tag=zone, constants=self.constants, fields=self.fields, dofs=dofs,spaces=spaces,numbering=numberings)
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

    problem = FeaSymMeca()

    problem.SetBulkFormulationToMeca(1.0,0.3,zone="3D")
    problem.AddPressure("Z1","p")

    problem.SetSpaceToLagrange()



    from BasicTools.FE.UnstructuredMeshTools import CreateMeshFromConstantRectilinearMesh
    import BasicTools.FE.ConstantRectilinearMesh as ConstantRectilinearMesh

    nx = 11; ny = 12; nz = 13;
    CRMesh = ConstantRectilinearMesh.ConstantRectilinearMesh()
    CRMesh.SetDimensions([nx,ny,nz]);
    CRMesh.SetSpacing([0.1, 0.1, 10./(nz-1)]);
    CRMesh.SetOrigin([0, 0, 0]);

    mesh = CreateMeshFromConstantRectilinearMesh(CRMesh)
    from BasicTools.FE.UnstructuredMeshTools import ComputeSkin
    ComputeSkin(mesh,inplace=True)

    # add tags to generate X1
    class ElementIteratorByZone(object):
        def __init__(self,zone,name):
            self.zone = zone
            self.tol = 1E-6
            self.tagname = name

        def applyOnElements(self,mesh,elements,el_id, allPoints=True):
            nodes = mesh.nodes[elements.connectivity[el_id,:],:]
            total = np.sum(self.zone(nodes) < self.tol)
            if allPoints:
                if total == elements.GetNumberOfNodesPerElement():
                    elements.tags.CreateTag(self.tagname,errorIfAlreadyCreated=False).AddToTag(el_id)
            else:
                if total > 0:
                    elements.tags.CreateTag(self.tagname,errorIfAlreadyCreated=False).AddToTag(el_id)


        def applyOnNodes(self,mesh):
            mask = self.zone(mesh.nodes) < self.tol
            for i in range(len(mask)):
                if mask[i]:
                    mesh.nodesTags.CreateTag(self.tagname,errorIfAlreadyCreated=False).AddToTag(i)


    mesh.ComputeBoundingBox()

    op0 = ElementIteratorByZone( lambda p: (p[:,2]-mesh.boundingMin[2]),"Z0" )
    op1 = ElementIteratorByZone( lambda p: (mesh.boundingMax[2]-p[:,2]),"Z1" )


    for name,el in mesh.elements.items():
        for i in range(el.GetNumberOfElements()):
            op0.applyOnElements(mesh,el,i)
            op1.applyOnElements(mesh,el,i)

    import BasicTools.FE.ElementNames as EN
    hexs = mesh.GetElementsOfType(EN.Hexaedron_8)
    hexs.GetTag("3D").SetIds(range(hexs.GetNumberOfElements()))

    print(mesh)
    op0.applyOnNodes(mesh)
    op1.applyOnNodes(mesh)
    print(mesh)


    problem.SetMesh(mesh)
    problem.ComputeDofNumbering()

    problem.constants["p"] = 1.0
    k,f = problem.GetLinearProblem()


    problem.AddDirichlet("Z0",0,0.)
    problem.AddDirichlet("Z0",1,0.)
    problem.AddDirichlet("Z0",2,0.)

    kk,ff = problem.ApplyBC(k,f)

    from BasicTools.FE.LinearSolver import LinearProblem
    """
    import BasicTools.FE.ElementNames  as EN
    dirichletids = mesh.GetElementsOfType(EN.Quadrangle_4).tags["Z0"].GetIds()
    dirichletdof = np.sort(np.unique(numbering[EN.Quadrangle_4][dirichletids,:].flatten() ))
    dirichletdofs = np.concatenate((dirichletdof,dirichletdof+numbering["size"], dirichletdof+2*numbering["size"]))

    #print(list(dirichletdofs))
    #print(list(f))

    dofs = np.ones(len(f), dtype=bool)
    dofs[dirichletdofs] = False

    import BasicTools.FE.FeaBase as FeaBase
    fixedValues = np.zeros(len(f),dtype=np.float)
    [K, rhsfixed] = FeaBase.deleterowcol(k, np.logical_not(dofs), np.logical_not(dofs), fixedValues)
    print(len(rhsfixed))
    rhs = f[dofs]-rhsfixed[dofs]

    prob = LinearProblem()
    prob.SetAlgo("Direct")


    prob.SetOp(K.tocsc())
    #print(list(dirichletdofs))
    print(k.shape )
    print(K.shape )
    print(rhs.shape )
    res = prob.Solve(rhs)
    #print(list(res))
    print("Done")
"""
    problem.solver = LinearProblem()
    problem.solver.SetAlgo("Direct")
    print(f)
    print(ff)

    problem.Solve(kk,ff)


    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        problem.PushVectorToMesh(True,f,"normalFlux")
        problem.PushVectorToMesh(True,problem.sol,"sol")
        """
        F = f.view()
        F.shape = (3,mesh.GetNumberOfNodes())
        F = F.T
        mesh.nodeFields["normalflux"] =  F[numbering["permutation"],:]
        fixedValues[dofs] = res
        fixedValues.shape = (3,mesh.GetNumberOfNodes())
        fixedValues= fixedValues.T
        mesh.nodeFields["sol"] =  fixedValues[numbering["permutation"],:]
        """
        OpenInParaView(mesh)
    return("ok")



if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
