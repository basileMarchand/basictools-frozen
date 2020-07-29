# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" Class to treat Constants Rectilinear Finit Element Problems

"""

import numpy as np
from scipy.sparse import coo_matrix
import scipy.sparse.linalg as linalg
import scipy.linalg as denselinalg
import  scipy.sparse as sps

import BasicTools.Containers.ElementNames as EN
import BasicTools.FE.FeaBase as FeaBase
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Linalg.MatOperations import deleterowcol


class BundaryCondition(BaseOutputObject):
    def __init__(self,dim=3, size= 1):
        super(BundaryCondition,self).__init__()
        self.sz = size
        self.nodes = np.empty((self.sz,dim),dtype=int)
        self.dofs = np.empty((self.sz,),dtype=int)
        self.vals = np.empty((self.sz,1),dtype=float)
        self.dim = dim
        self.cpt = 0

    def reserve(self, size):
        self.nodes = np.resize(self.nodes, (size,self.dim))
        self.dofs = np.resize(self.dofs, ( size,))
        self.vals = np.resize(self.vals, ( size,1))

        self.sz = size

    def tighten(self):
        self.reserve(self.cpt)

    def eliminate_double(self, overwrite=True):

        if len(self.nodes) == 0:
            return
        m = np.amax(self.nodes,axis=0)+1
        #print(m)

        pointers = [None]*(np.prod(m)*3)
        m[1:] *= m[0]
        m *= 3
        cachesize = np.insert(np.delete(m,len(m)-1),0,3)

        i = 0 ;
        while (i < self.cpt):

          node = self.nodes[i]
          pp = sum(node*cachesize) + self.dofs[i]

          if pointers[pp] is None:
             pointers[pp] = i
             i +=1
          else:
             p = pointers[pp]

             # we copy the value of the last encountered val
             if overwrite :
                 self.vals[p] = self.vals[i]

             #then we move the last value to the current place
             self.nodes[i] = self.nodes[self.cpt-1]
             self.dofs[i] = self.dofs[self.cpt-1]
             self.vals[i] = self.vals[self.cpt-1]
             self.cpt -= 1


        self.tighten()

    def append(self, nodes, dof,val):
        if self.cpt >= self.sz:
            self.reserve(self.sz*2)
        self.nodes[self.cpt,:] = nodes
        self.dofs[self.cpt] = dof
        self.vals[self.cpt] = val
        self.cpt += 1

    def __str__(self):
       res = ""
       i = 0
       while (i < self.cpt):
           res += str(self.nodes[i])+ " "
           res += str(self.dofs[i]) + " "
           res += str(self.vals[i])+ " \n"
           i +=1
       return res

class ElementaryMatrix():
    def __init__(self,dim=3, physics="disp"):
        self.dim = dim
        self.geoFactor = [1,1,1]
        self.physics = physics
        self.young = 1
        self.poisson = 0.3
        self.density = 1
        self.thermalK = 1

    def GetMassMatrix(self):
        if self.physics == "disp":
            from BasicTools.FE.SymPhysics import MecaPhysics
            phys = MecaPhysics(self.dim)
            wform = phys.GetBulkMassFormulation()
            return self.Integrate(["u_"+str(x) for x in range(self.dim)],wform)
        elif self.physics == "thermal":
            from BasicTools.FE.SymPhysics import  ThermalPhysics
            wform = ThermalPhysics().GetBulkMassFormulation()
            return self.Integrate(["t"],wform)

    def GetTangetMatrix(self):
        from BasicTools.FE.SymPhysics import MecaPhysics,ThermalPhysics

        if self.physics == "disp":
            physics = MecaPhysics(dim=self.dim)
            wform = physics.GetBulkFormulation(self.young,self.poisson)
        elif self.physics == "thermal":
            physics = ThermalPhysics()
            physics.spaceDimension = self.dim
            wform = physics.GetBulkFormulation(alpha=self.thermalK)
        else:
            raise(Exception("Physics not suppported"))
        return self.Integrate(physics.GetPrimalNames() ,wform)

    def Integrate(self,primalNames,wform):
        from BasicTools.FE.FETools import GetElementaryMatrixForFormulation
        if self.dim == 3:
            el = EN.Hexaedron_8
        else:
            el = EN.Quadrangle_4

        KGeneric,rhs = GetElementaryMatrixForFormulation(el,wform, unknownNames = primalNames, geoFactor= self.geoFactor)

        n = len(primalNames)
        permutation = np.array([x*EN.numberOfNodes[el]+y for y in range(EN.numberOfNodes[el]) for x in range(n)])
        rhs = rhs[permutation]
        KGeneric = KGeneric[permutation,:]
        KGeneric = KGeneric[:,permutation]
        return KGeneric

class Fea(FeaBase.FeaBase):

    def __init__(self):
        super(Fea,self).__init__()
        self.linearSolver = "EigenCG"
        self.writer = None
        self.minthreshold = 0.9e-3
        self.tol = 1.e-6
        self.dofpernode = 1
        self.init = False
        self.dirichlet_bcs=None
        self.neumann_bcs=None
        self.neumann_nodal=None
        self.young = 1.
        self.poisson = 0.3
        self.density = 1.

    def BuildProblem(self,support=None, dofpernode = None, dirichlet_bcs = None, neumann_bcs= None, KOperator= None, MOperator= None, neumann_nodal= None):
        if self.init == True:
            return

        self.init = True

        if support is not None:
            self.support = support ;

        if self.support.IsConstantRectilinear() == False :
            raise Exception("Must be a ConstantRectilinear mesh type ") #pragma: no cover

        self.outer_v = []

        if dofpernode is not None:
            self.dofpernode = dofpernode

        if dirichlet_bcs is not None:
            self.dirichlet_bcs = dirichlet_bcs

        if neumann_bcs is not None:
            self.neumann_bcs = neumann_bcs

        if neumann_bcs is not None:
            self.neumann_bcs = neumann_bcs

        if neumann_nodal is not None:
            self.neumann_nodal = neumann_nodal


        self.nodesPerElement = 2**self.support.GetDimensionality()

        if KOperator is not None:
            self.KE    = KOperator
            self.ME    = MOperator
        else:
            self.myElem = ElementaryMatrix(dim = self.support.GetDimensionality())
            self.myElem.young = self.young
            self.myElem.poisson = self.poisson
            self.myElem.geoFactor = self.support.GetSpacing()
            self.KE = self.myElem.GetTangetMatrix()
            self.ME = self.myElem.GetMassMatrix()

        # dofs:
        self.ndof = self.dofpernode * self.support.GetNumberOfNodes()

        # FE: Build the index vectors for the for coo matrix format
        self.PrintDebug("Building Connectivity matrix")
        self.edofMat = np.zeros((self.support.GetNumberOfElements(), self.nodesPerElement*self.dofpernode), dtype=np.int_)
        self.PrintDebug("Building Connectivity matrix 2")
        for i in  range(self.support.GetNumberOfElements()):
            coon = self.support.GetConnectivityForElement(i)
            self.edofMat[i, :] = np.array([(coon*self.dofpernode+y) for y in range(self.dofpernode)]).ravel('F')

        self.PrintDebug("Building Connectivity matrix Done")

        self.iK = None
        self.jK = None

        self.fixedValues = np.zeros((self.ndof, 1), dtype=np.double)

        self.PrintDebug("Treating Dirichlet 1/4")

        self.fixed = np.zeros(self.ndof, dtype=np.bool)
        if self.dirichlet_bcs is not None :
            self.dirichlet_bcs.tighten()
            self.dirichlet_bcs.eliminate_double()
            indexs = self.support.GetMonoIndexOfNode(self.dirichlet_bcs.nodes)
            indexs *= self.dofpernode
            indexs += self.dirichlet_bcs.dofs
            self.fixed[indexs] = True
            self.fixedValues[self.fixed.T,0:] = self.dirichlet_bcs.vals

        self.free = np.ones(self.ndof, dtype=np.bool)
        self.free[self.fixed] = False

        # Solution and RHS vectors
        self.f = np.zeros((self.ndof, 1), dtype=np.double)
        self.u = np.zeros((self.ndof, 1), dtype=np.double)
        self.PrintDebug("Treating Dirichlet Done")
        self.PrintDebug("Treating Neumann")

            # Set load
        if  self.neumann_bcs is not  None:

            self.neumann_bcs.tighten()
            if self.neumann_bcs.cpt:
              self.support.GenerateFullConnectivity()
              z = np.zeros((self.support.GetNumberOfNodes(),))
              z[self.support.GetMonoIndexOfNode(self.neumann_bcs.nodes)] +=  1.;
              eff = np.clip((np.sum(z[self.support.GenerateFullConnectivity()],axis=1) ),0, 1)

              MassMatrix = self.BuildMassMatrix(eff)
              self.f[self.support.GetMonoIndexOfNode(self.neumann_bcs.nodes)*self.dofpernode + self.neumann_bcs.dofs] += self.neumann_bcs.vals
              self.f[:,0] = MassMatrix*self.f[:,0]

        if  self.neumann_nodal is not  None:
            self.neumann_nodal.tighten()

            nodal_f = np.zeros((self.ndof, 1), dtype=np.double)

            nodal_f[self.support.GetMonoIndexOfNode(self.neumann_nodal.nodes)*self.dofpernode + self.neumann_nodal.dofs] += self.neumann_nodal.vals
            self.f[:,0] += nodal_f[:,0]
        self.PrintDebug("Treating Neumann Done")

        self.eed = np.zeros(self.support.GetNumberOfElements())

    def AssemblyMatrix(self,Op, Eeff = None):
        if Eeff is None:
            Eeff = np.ones(self.support.GetNumberOfElements())
            sM = ((Op.flatten()[np.newaxis]).T * Eeff.ravel()).flatten(order='F')
            self.GenerateIJs()
            M = coo_matrix((sM, (self.iK, self.jK)), shape=(self.ndof, self.ndof),dtype=float)
        else:
            self.PrintDebug(" Eeff is known")
            bool_Eeff = (Eeff>=self.minthreshold)
            nEeff = Eeff[bool_Eeff]
            sM = ((Op.flatten()[np.newaxis]).T * nEeff.ravel()).flatten(order='F')
            one = np.ones((self.nodesPerElement*self.dofpernode, 1), dtype=np.int_)
            local_iK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()
            one.shape = (1,self.nodesPerElement*self.dofpernode)
            local_jK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()
            M = coo_matrix((sM, (local_iK, local_jK)), shape=(self.ndof, self.ndof),dtype=float).tocsr()
        return M.tocsr()

    def BuildTangentMatrix(self, Eeff = None):
        self.PrintDebug("BuildTangentMatrix")
        res = self.AssemblyMatrix(self.KE, Eeff)
        self.PrintDebug("BuildTangentMatrix Done")
        return res

    def BuildMassMatrix(self, Eeff = None):
        self.PrintDebug("BuildMassMatrix")
        res = self.AssemblyMatrix(self.ME, Eeff)
        self.PrintDebug("BuildMassMatrix Done")
        return res


    def Solve(self, Eeff=None):

        # hack to integrate complex boundary condition in the mesh
        if hasattr(self,"mecaPhysics") and (self.mecaPhysics is not None):
            from BasicTools.FE.UnstructuredFeaSym import UnstructuredFeaSym
            prob = UnstructuredFeaSym()
            prob.physics.append(self.mecaPhysics)
            prob.SetMesh(self.support)
            self.mecaPhysics.ComputeDofNumberingFromConnectivity(self.support)
            prob.ComputeDofNumbering(self.support)
            k,f = prob.GetLinearProblem(computeK=False)
            if f is not None:
                f.shape = (3,len(f)//3)
                f = f.ravel(order='F')
                self.f[:,0] += f
            self.mecaPhysics = None

        self.PrintDebug("Construction of the tangent matrix")

        K = self.BuildTangentMatrix(Eeff)

        zerosdof = np.where(K.diagonal()== 0 )[0]

        self.PrintVerbose("Number of active nodes : " + str(self.ndof-len(zerosdof) ) + "  of " + str(self.ndof) + "   "+ str(float(len(zerosdof)*100.)/self.ndof)+ "% of empty dofs"  )
        Kones = coo_matrix( (np.ones((len(zerosdof),) ) ,(zerosdof,zerosdof)), shape =(self.ndof, self.ndof)).tocsr()#(self.dofpernode,self.dofpernode))
        K = (K.tocsr() + Kones.tocsr()).tocsr()

        # Remove constrained dofs from matrix
        self.PrintDebug(" Delete fixed Dofs")
        [K, rhsfixed] = deleterowcol(K, self.fixed, self.fixed, self.fixedValues)

        self.PrintDebug(" Start solver (" + str(self.linearSolver) + ")")
        rhs = self.f[self.free, 0]-rhsfixed[self.free, 0]

        self.u = np.zeros((self.ndof, 1), dtype=np.double)

        if K.nnz > 0 :
          if self.linearSolver == "Direct":
            self.u[self.free, 0] = linalg.spsolve(K, rhs)
          elif self.linearSolver == "DirectDense":
            self.u[self.free, 0] = denselinalg.solve(K.toarray(), rhs,sym_pos=True,overwrite_a=True)
          elif self.linearSolver == "CG":
            # Preconditioning
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            # normalization of the rhs term ( to treat correctly the tol of CG) (Please read the documentaion of numpy.linalg.cg)
            norm = np.linalg.norm(rhs)
            #def PrintRes(x):
            #   self.PrintDebug(np.linalg.norm(K.dot(x)-rhs/norm))
            #    , callback=PrintRes
            self.PrintDebug("tol " + str(self.tol))
            res = linalg.cg(K.tocsc(copy=False), rhs/norm, x0 = self.u[self.free, 0]/norm , M = M, tol = self.tol , atol=self.tol)
            if res[1] > 0:
                raise
            self.u[self.free, 0] = res[0]*norm

          elif self.linearSolver == "LGMRES":
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            norm = np.linalg.norm(rhs)
            res = linalg.lgmres(K, rhs/norm, x0 = self.u[self.free, 0]/norm , M = M , tol = self.tol, atol=self.tol )
            self.u[self.free, 0] = res[0]*norm
          elif self.linearSolver == "AMG":# pragma: no cover
            try:
                import pyamg
            except:#pragma: no cover
                raise Exception('AMG module not installed')#pragma: no cover
            K = K.tocsr()
            ml = pyamg.ruge_stuben_solver(K)
            res = ml.solve(rhs,x0 = self.u[self.free, 0] , tol=1e-12,accel='cg')
            self.u[self.free, 0] = res
          elif self.linearSolver == "cholesky":# pragma: no cover
            from sksparse.cholmod import cholesky
            self.u[self.free, 0] = cholesky(K)(rhs )

          else :
            from BasicTools.Linalg.LinearSolver import LinearProblem
            linSol = LinearProblem()
            linSol.tol = self.tol
            linSol.SetAlgo(self.linearSolver)
            linSol.SetOp(K)
            self.u[self.free, 0] = linSol.Solve(rhs)
          #else :

            #print("'"+self.linearSolver + "' is not a valid linear solver")#pragma: no cover
            #print('Please set a type of linear solver')#pragma: no cover
            #raise Exception()#pragma: no cover

        self.PrintDebug('Post Process')
        self.u = self.u + self.fixedValues

        # Compute element elastic energy density
        u_reshaped = self.u[self.edofMat]
        u_reshaped.shape = (self.support.GetNumberOfElements(), self.nodesPerElement*self.dofpernode)
        Ku_reshaped = np.dot(u_reshaped, self.KE)
        np.einsum('ij,ij->i', Ku_reshaped, u_reshaped, out=self.eed)
        # we divide by the volume of one element to get the energy density
        self.eed /= np.prod(self.support.GetSpacing())
        self.PrintDebug('Post Process Done')

    def GenerateIJs(self):
        # lazy generations of the IJ for the case when no density is given
        self.PrintDebug("GenerateIJs")
        if self.iK is None:
                nodesPerElement = 2**self.support.GetDimensionality()
                ones = np.ones((nodesPerElement*self.dofpernode, 1), dtype=np.int_)
                self.iK = np.kron(self.edofMat, ones).flatten()
                ones.shape = (1, nodesPerElement*self.dofpernode)
                self.jK = np.kron(self.edofMat, ones).flatten()
        self.PrintDebug("GenerateIJs Done")

    def element_elastic_energy(self, Eeff= None, OnlyOnInterface = False):

        if Eeff is None:
            return self.eed
        else:
            nEeff = np.copy(Eeff)
            bool_Eeff = Eeff<self.minthreshold
            nEeff[bool_Eeff] = 0.0
            return nEeff*self.eed


    def nodal_elastic_energy(self, Eeff=None, OnlyOnInterface = False):

        if Eeff is None:
            Eeff = np.ones(self.support.GetNumberOfElements())

        return node_averaged_element_field(self.element_elastic_energy(Eeff,OnlyOnInterface=OnlyOnInterface),self.support)

    def Write(self):

        if self.writer is not None:
            self.writer.Write(self.support,
                PointFields     = [self.u, self.f],
                PointFieldsNames= ["u", "f"]
                )

def element_averaged_node_field(node_field,support):
    nnodes = support.GetDimensions()
    ndims = support.GetDimensionality()
    if ndims == 3:
        result = np.zeros((nnodes[0]-1,nnodes[1]-1,nnodes[2]-1 ))

        field  = node_field.view()
        field.shape = tuple(x for x in support.GetDimensions())

        result += field[0:-1, 0:-1,0:-1]
        result += field[0:-1, 1:  ,0:-1]
        result += field[1:  , 0:-1,0:-1]
        result += field[1:  , 1:  ,0:-1]

        result += field[0:-1, 0:-1,1:  ]
        result += field[0:-1, 1:  ,1:  ]
        result += field[1:  , 0:-1,1:  ]
        result += field[1:  , 1:  ,1:  ]

        return result/8
    elif ndims == 2:
        result = np.zeros((nnodes[0]-1,nnodes[1]-1))

        field  = node_field.view()
        field.shape = tuple(x for x in support.GetDimensions())

        result += field[0:-1, 0:-1]
        result += field[0:-1, 1:  ]
        result += field[1:  , 0:-1]
        result += field[1:  , 1:  ]

        return result/4
    else :
        raise(Exception("only implemented for dim = 3 or 2"))

def node_averaged_element_field(element_field,support):
    nnodes = support.GetDimensions()
    ndims = support.GetDimensionality()
    if ndims == 3:
        result = np.zeros((nnodes[0],nnodes[1],nnodes[2] ))
        w = np.zeros((nnodes[0],nnodes[1],nnodes[2] ))

        field  = element_field.view()
        field.shape = tuple(x-1 for x in support.GetDimensions())

        result[0:-1, 0:-1,0:-1] += field
        result[0:-1, 1:  ,0:-1] += field
        result[1:  , 0:-1,0:-1] += field
        result[1:  , 1:  ,0:-1] += field

        result[0:-1, 0:-1,1:  ] += field
        result[0:-1, 1:  ,1:  ] += field
        result[1:  , 0:-1,1:  ] += field
        result[1:  , 1:  ,1:  ] += field

        w[0:-1, 0:-1,0:-1] += 1
        w[0:-1, 1:  ,0:-1] += 1
        w[1:  , 0:-1,0:-1] += 1
        w[1:  , 1:  ,0:-1] += 1

        w[0:-1, 0:-1,1:  ] += 1
        w[0:-1, 1:  ,1:  ] += 1
        w[1:  , 0:-1,1:  ] += 1
        w[1:  , 1:  ,1:  ] += 1

        return result/w
    else:
        result = np.zeros((nnodes[0],nnodes[1] ))
        w = np.zeros((nnodes[0],nnodes[1] ))

        field  = element_field.view()
        field.shape = tuple(x-1 for x in support.GetDimensions())

        result[0:-1, 0:-1] += field
        result[0:-1, 1:  ] += field
        result[1:  , 0:-1] += field
        result[1:  , 1:  ] += field

        w[0:-1, 0:-1] += 1
        w[0:-1, 1:  ] += 1
        w[1:  , 0:-1] += 1
        w[1:  , 1:  ] += 1

        return result/w

def CheckIntegrityThermal3D():
    print('----------------------------  Thermal3D -------------------------------------------------')

    import BasicTools.Containers.ConstantRectilinearMesh as CRM
    import BasicTools.IO.XdmfWriter  as XdmfWriter
    import time
    from BasicTools.Helpers.Tests import TestTempDir

    myMesh = CRM.ConstantRectilinearMesh()
    nx = 11; ny = 12; nz = 13;


    myMesh.SetDimensions([nx,ny,nz]);
    myMesh.SetSpacing([0.1, 0.1, 10./(nz-1)]);
    myMesh.SetOrigin([0, 0, 0]);
    print(myMesh)

    # thermal problem
    #dirichlet at plane z =0
    dirichlet_bcs = BundaryCondition()

    for x in range(nx):
        for y in range(ny):
            for z in [0]:
                dirichlet_bcs.append([x,y,z],0 , 0 )


    # Homogenous body flux
    neumann_bcs = BundaryCondition()
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                neumann_bcs.append([x,y,z],0 , 1 )

    neumann_bcs.append([0,0,0],0 , 1 )
    neumann_bcs.append([1,1,0],0 , 1 )
    neumann_bcs.eliminate_double()

    starttime = time.time()
    myElement = ElementaryMatrix()
    myElement.physics = "thermal"
    myElement.geoFactor = myMesh.GetSpacing()

    myProblem = Fea()
    myProblem.BuildProblem(myMesh, dofpernode = 1,
        KOperator = myElement.GetTangetMatrix(),
        MOperator = myElement.GetMassMatrix(),
        dirichlet_bcs = dirichlet_bcs,
        neumann_bcs = neumann_bcs)

    print("Time for Fea definition : " + str(time.time() -starttime))


    myProblem.writer = XdmfWriter.XdmfWriter(TestTempDir.GetTempPath()+"Test_constantRectilinearThermal.xdmf")
    myProblem.writer.SetBinary()
    myProblem.writer.SetTemporal()
    myProblem.writer.Open()

    myProblem.linearSolver = 'Direct'
    myProblem.Solve()
    myProblem.Write()
    myProblem.linearSolver = 'DirectDense'
    myProblem.Solve()
    myProblem.Write()
    myProblem.linearSolver = 'LGMRES'
    myProblem.Solve()
    myProblem.element_elastic_energy()
    myProblem.Write()
    try:# pragma: no cover
       import pyamg
       myProblem.linearSolver = 'AMG'
       myProblem.Solve()
       myProblem.Write()
    except Exception as inst: # pragma: no cover
       print(inst.args[0])# pragma: no cover

    myProblem.writer.Close()


    path = TestTempDir.GetTempPath() +'TestThermal3D.xmf'
    XdmfWriter.WriteMeshToXdmf(path,myMesh,
                    [myProblem.u, myProblem.f],
                    PointFieldsNames= ['Themperature','q'],
                    GridFieldsNames=[])
    print('DONE')
    print(max(myProblem.u))

    if abs(max(myProblem.u)-50.) > 1e-4:
        raise Exception()# pragma: no cover
    return("ok")


def CheckIntegrityDep3D():

    import time
    import BasicTools.Containers.ConstantRectilinearMesh as CRM
    import BasicTools.IO.XdmfWriter  as XdmfWriter

    from BasicTools.Helpers.Tests import TestTempDir

    print('------------------------------- Dep3D ----------------------------------------------')

    nx = 11; ny = 12; nz = 13;
    myMesh = CRM.ConstantRectilinearMesh()
    myMesh.SetDimensions([nx,ny,nz]);
    myMesh.SetSpacing([0.1, 0.1, 10./(nz-1)]);
    myMesh.SetOrigin([0, 0, 0]);

    # block all the faces rith
    dirichlet_bcs = BundaryCondition()

    for x in range(nx):
        for y in range(ny):
            for coor in range(3):
                for z in [0]:
                    dirichlet_bcs.append([x,y,z],coor , 0 )
                for z in [nz-1]:
                    dirichlet_bcs.append([x,y,z],coor , 1 )
#    dirichlet_bcs =( [(x, y, z, coor, 0) for x in range(nx) for y in range(ny) for z in [0]  for coor in range(3)] +
#                    [(x, y, z, coor, 1) for x in range(nx) for y in range(ny) for z in [nz-1]  for coor in range(3)])


    neumann_bcs = BundaryCondition()
    neumann_bcs.append([0,0,0],0,0) # ([x,y],coor, value)

    neumann_nodal = BundaryCondition()
    neumann_nodal.append([0,0,0],0,0) # ([x,y],coor, value)

    starttime = time.time()

    myProblem = Fea()
    myProblem.BuildProblem(myMesh, dofpernode = 3,
        dirichlet_bcs = dirichlet_bcs,
        neumann_bcs = neumann_bcs,
        neumann_nodal = neumann_nodal)

    print("Time for Fea definition : " + str(time.time() -starttime))

    starttime = time.time()
    densities  = np.ones(myMesh.GetNumberOfElements());
    myProblem.Solve(densities)
    print("Time for Fea solve : " + str(time.time() -starttime))

    XdmfWriter.WriteMeshToXdmf(TestTempDir.GetTempPath() +'TestDep3D.xmf',myMesh,
                    [myProblem.u, myProblem.f, myProblem.nodal_elastic_energy() ],
                    [densities, myProblem.element_elastic_energy() ] ,
                    [],
                    PointFieldsNames= ['Dep','F','NEnergie'],
                    CellFieldsNames=['densities','EEnergie'],
                    GridFieldsNames=[])

    print(max(myProblem.u))

    if abs(max(myProblem.u)-1.00215295) > 1e-5:
        print(TestTempDir.GetTempPath())
        raise   Exception("The value must be 1.00215295")# pragma: no cover

def CheckIntegrityThermal2D():
    import BasicTools.Containers.ConstantRectilinearMesh as CRM
    import BasicTools.IO.XdmfWriter  as XdmfWriter
    import time
    from BasicTools.Helpers.Tests import TestTempDir

    print('----------------------------  Thermal2D -------------------------------------------------')

    myMesh = CRM.ConstantRectilinearMesh(2)
    nx = 11; ny = 11;

    myMesh.SetDimensions([nx,ny]);
    myMesh.SetSpacing([0.1, 0.1]);
    myMesh.SetOrigin([0, 0]);
    print(myMesh)

    dirichlet_bcs = BundaryCondition(dim = 2)
    neumann_bcs = BundaryCondition(dim = 2)

    # thermal problem
    #dirichlet at plane z =0
    # Homogenous body flux
    for x in range(nx):
        for y in [0]:
            dirichlet_bcs.append([x,y],0,0)
        for y in range(ny):
            neumann_bcs.append([x,y],0,1)


    starttime = time.time()
    myElement = ElementaryMatrix(dim=2,physics="thermal")
    myElement.geoFactor = myMesh.GetSpacing()

    myProblem = Fea()
    myProblem.BuildProblem(myMesh, dofpernode = 1,
        KOperator = myElement.GetTangetMatrix(),
        MOperator = myElement.GetMassMatrix(),
        dirichlet_bcs = dirichlet_bcs,
        neumann_bcs = neumann_bcs)

    print("Time for Fea definition : " + str(time.time() -starttime))

    myProblem.writer = XdmfWriter.XdmfWriter(TestTempDir.GetTempPath()+"TestProblemWriterThermal2D.xdmf")
    myProblem.writer.SetBinary()
    myProblem.writer.SetTemporal()
    myProblem.writer.Open()

    myProblem.linearSolver = 'Direct'
    myProblem.Solve()
    myProblem.Write()
    myProblem.linearSolver = 'DirectDense'
    myProblem.Solve()
    myProblem.Write()
    myProblem.linearSolver = 'LGMRES'
    myProblem.Solve()
    myProblem.element_elastic_energy()
    myProblem.Write()
    try:# pragma: no cover
       import pyamg
       myProblem.linearSolver = 'AMG'
       myProblem.Solve()
       myProblem.Write()
    except Exception as inst: # pragma: no cover
        print(inst.args[0])# pragma: no cover

    myProblem.writer.Close()


    path = TestTempDir.GetTempPath() +'TestThermal2D.xmf'
    XdmfWriter.WriteMeshToXdmf(path,myMesh,
                    [myProblem.u, myProblem.f],
                    PointFieldsNames= ['Themperature','q'],
                    GridFieldsNames=[])
    print('DONE')
    print(max(myProblem.u))

    if abs(max(myProblem.u)-.5) > 1e-5:
        raise Exception()# pragma: no cover

    return 'OK'
    #dirichlet_bcs =( [(0, y, z, cor) for y in range(ny) for z in range(nz) for cor in range(3)] )
    #neumann_bcs = ([(nx-1, y, z, 2) for y in range(ny) for z in range(nz) ])
    #Fea(myMesh, dofpernode = 1, Operator= None, dirichlet_bcs = dirichlet_bcs, neumann_bcs =neumann_bcs):

def CheckIntegrityDep2D():

    import BasicTools.Containers.ConstantRectilinearMesh as CRM
    import BasicTools.IO.XdmfWriter  as XdmfWriter
    import time
    from BasicTools.Helpers.Tests import TestTempDir

    print('----------------------- 2D dep ------------------------------------------------------')
    myMesh = CRM.ConstantRectilinearMesh(2)
    nx = 11; ny = 12;

    myMesh.SetDimensions([nx,ny]);
    myMesh.SetSpacing([1./(nx-1), 1./(ny-1)]);
    myMesh.SetOrigin([0, 0]);
    print(myMesh)
    # block all the faces rith

    dirichlet_bcs = BundaryCondition(dim=2)
    for x in range(nx):
        for y in [0]:
            for coor in range(2):
                dirichlet_bcs.append([x,y], coor, 0 )
        for y in [ny-1]:
            for coor in range(2):
                dirichlet_bcs.append([x,y], coor, 1 )

    neumann_bcs = BundaryCondition(dim=2)
    neumann_bcs.append([0,0],0,0) # ([x,y],coor, value)

    neumann_nodal = BundaryCondition(dim=2)
    neumann_nodal.append([0,0],0,0) # ([x,y],coor, value)

    starttime = time.time()

    myProblem = Fea()
    myProblem.BuildProblem(myMesh, dofpernode = 2,
        dirichlet_bcs = dirichlet_bcs,
        neumann_bcs = neumann_bcs,
        neumann_nodal = neumann_nodal)

    print("Time for Fea definition : " + str(time.time() -starttime))

    starttime = time.time()
    densities  = np.ones(myMesh.GetNumberOfElements());
    myProblem.Solve(densities)

    print("Time for Fea solve : " + str(time.time() -starttime))
    print(myProblem.u.T)

    # build mass matrix
    myProblem.BuildMassMatrix();
    myProblem.BuildMassMatrix(densities);

    XdmfWriter.WriteMeshToXdmf(TestTempDir.GetTempPath() +'TestDep2D.xmf',myMesh,
                    [myProblem.u, myProblem.f, myProblem.nodal_elastic_energy(), myProblem.fixed.astype(int) ],
                    [densities, myProblem.element_elastic_energy() ] ,
                    [],
                    PointFieldsNames= ['Dep','F','NEnergie', 'ufixed'],
                    CellFieldsNames=['densities','EEnergie'],
                    GridFieldsNames=[])

    return 'ok'

def CheckIntegrity():

    print(CheckIntegrityThermal3D())
    print(CheckIntegrityDep3D())
    print(CheckIntegrityThermal2D())
    print(CheckIntegrityDep2D())

    from BasicTools.Helpers.Tests import TestTempDir
    print(TestTempDir.GetTempPath())
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
