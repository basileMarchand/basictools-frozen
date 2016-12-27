# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np

from scipy.sparse import coo_matrix
import scipy.sparse.linalg as linalg
import scipy.linalg as denselinalg

import  scipy.sparse as sps
from OTTools.FE.Hexa8Cuboid import Hexa8Cuboid
from OTTools.FE.Quad4Rectangle import  Quad4Rectangle

import OTTools.FE.FeaBase as FeaBase
from OTTools.Helpers.BaseOutputObject import BaseOutputObject

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

    def eliminate_double(self):

        m = np.amax(self.nodes, axis=0)
        fcpt = np.zeros(m+1,dtype=np.int);

        i = 0 ;
        self.PrintDebug(self.cpt)
        fac = 100./self.cpt
        c = 0
        while (i < self.cpt):
          if i*fac >c:
              self.Print2(round(i*fac))
              c +=1

          node = self.nodes[i]
          if fcpt[node[0],node[1],node[2]] == 0:
              fcpt[node[0],node[1],node[2]] += 1
              i += 1
              continue
          fcpt[node[0],node[1],node[2]] += 1

          dofs = self.dofs[i]
          #val = self.vals[i]
          j = 0
          while(j < i-1 and i < self.cpt ):
             xyz =self.nodes[j]
             if ( xyz[0] == node[0] and xyz[1] == node[1] and xyz[2] == node[2]  and self.dofs[j] == dofs   ):
                 #(self.nodes[j] == node).all()
                 self.nodes[j] = self.nodes[self.cpt-1]
                 self.dofs[j] = self.dofs[self.cpt-1]
                 self.vals[j] = self.vals[self.cpt-1]
                 self.cpt -=1
                 continue
             j +=1

          i +=1
        self.tighten()

    def append(self, nodes, dof,val):
        if self.cpt >= self.sz:
            self.reserve(self.sz*2)
        self.nodes[self.cpt,:] = nodes
        self.dofs[self.cpt] = dof
        self.vals[self.cpt] = val
        self.cpt += 1

class Fea(FeaBase.FeaBase):

    def __init__(self, support, dofpernode = 1, dirichlet_bcs = None, neumann_bcs= None, KOperator= None, MOperator= None,  neumann_nodal= None):
        super(Fea,self).__init__()

        if support.IsConstantRectilinear() == False :
            raise Exception("Must be a ConstantRectilinear mesh type ") #pragma: no cover

        self.linearSolver = "CG"
        self.outer_v = []

        self.support = support ;

        self.dofpernode = dofpernode
        self.nodesPerElement = 2**self.support.GetDimensionality()

        self.writer = None

        self.minthreshold = 0.9e-3
        self.tol = 1.e-6

        if KOperator is not None:
            self.KE    = KOperator
            self.ME    = MOperator
        else:
            if support.GetDimensionality() == 3:
                myElem = Hexa8Cuboid()
            else:
               myElem = Quad4Rectangle()
            myElem.delta = support.GetSpacing()
            self.KE = myElem.GetIsotropDispK(1.,0.3);
            self.ME = myElem.GetIsotropDispM(1.);

        # dofs:
        self.ndof = dofpernode * support.GetNumberOfNodes()

        # FE: Build the index vectors for the for coo matrix format
        self.PrintDebug("Building Connectivity matrix")
        self.edofMat = np.zeros((support.GetNumberOfElements(), self.nodesPerElement*dofpernode), dtype=np.int_)
        self.PrintDebug("Building Connectivity matrix 2")
        for i in  xrange(support.GetNumberOfElements()):
            coon = support.GetConnectivityForElement(i)
            self.edofMat[i, :] = np.array([(coon*dofpernode+y) for y in xrange(dofpernode)]).flatten('F')
        self.PrintDebug("Building Connectivity matrix Done")

        self.iK = None
        self.jK = None

        self.fixedValues = np.zeros((self.ndof, 1), dtype=np.double)

        self.PrintDebug("Treating Dirichlet 1/4")

        self.fixed = np.zeros(self.ndof, dtype=np.bool)
        if dirichlet_bcs is not None :
            dirichlet_bcs.tighten()
            indexs = support.GetMonoIndexOfNode(dirichlet_bcs.nodes)
            indexs *= dofpernode
            indexs += dirichlet_bcs.dofs
            self.fixed[indexs] = True
            self.fixedValues[self.fixed.T,0:] = dirichlet_bcs.vals

        self.free = np.ones(self.ndof, dtype=np.bool)
        self.free[self.fixed] = False

        # Solution and RHS vectors
        self.f = np.zeros((self.ndof, 1), dtype=np.double)
        self.u = np.zeros((self.ndof, 1), dtype=np.double)
        self.PrintDebug("Treating Dirichlet Done")
        self.PrintDebug("Treating Neumann")

            # Set load
        if  neumann_bcs is not  None:

            neumann_bcs.tighten()
            self.support.GenerateFullConnectivity()
            z = np.zeros((self.support.GetNumberOfNodes(),))
            z[support.GetMonoIndexOfNode(neumann_bcs.nodes)] +=  1.;
            eff = np.clip((np.sum(z[self.support.connectivity],axis=1) ),0, 1)

            MassMatrix = self.BuildMassMatrix(eff)
            self.f[support.GetMonoIndexOfNode(neumann_bcs.nodes)*dofpernode + neumann_bcs.dofs] += neumann_bcs.vals
            self.f[:,0] = MassMatrix*self.f[:,0]

        if  neumann_nodal is not  None:
            neumann_nodal.tighten()

            nodal_f = np.zeros((self.ndof, 1), dtype=np.double)

            nodal_f[support.GetMonoIndexOfNode(neumann_nodal.nodes)*dofpernode + neumann_nodal.dofs] += neumann_nodal.vals
            self.f[:,0] += nodal_f[:,0]
        self.PrintDebug("Treating Neumann Done")

        self.eed = np.zeros(support.GetNumberOfElements())

    def BuildMassMatrix(self, Eeff = None):

        self.PrintDebug("BuildMassMatrix")
        if Eeff is None:
            self.PrintDebug(" Eeff is None")
            Eeff = np.ones(self.support.GetNumberOfElements())
            sM = ((self.ME.flatten()[np.newaxis]).T * Eeff.ravel()).flatten(order='F')
            self.GenerateIJs()
            self.PrintDebug("Asm")
            M = coo_matrix((sM, (self.iK, self.jK)), shape=(self.ndof, self.ndof),dtype=float)
            self.PrintDebug("BuildMassMatrix Done")
            return  M.tocsr()#(self.dofpernode,self.dofpernode))
        else:
            self.PrintDebug(" Eeff is known")
            bool_Eeff = (Eeff>=self.minthreshold)
            nEeff = Eeff[bool_Eeff]
            sM = ((self.ME.flatten()[np.newaxis]).T * nEeff.ravel()).flatten(order='F')

            one = np.ones((self.nodesPerElement*self.dofpernode, 1), dtype=np.int_)
            local_iK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()
            one.shape = (1,self.nodesPerElement*self.dofpernode)
            local_jK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()
            self.PrintDebug("Asm")
            M = coo_matrix((sM, (local_iK, local_jK)), shape=(self.ndof, self.ndof),dtype=float).tocsr()

            self.PrintDebug("BuildMassMatrix Done")
            return M


    def Solve(self, Eeff=None):

        self.PrintDebug("Construction of the tangent matrix")
        if Eeff is None:
            Eeff = np.ones(self.support.GetNumberOfElements())
            sK = ((self.KE.flatten()[np.newaxis]).T * Eeff.ravel()).flatten(order='F')
            self.GenerateIJs()
            K = coo_matrix((sK, (self.iK, self.jK)), shape=(self.ndof, self.ndof)).tocsr()#(self.dofpernode,self.dofpernode))
        else:
            bool_Eeff = Eeff>self.minthreshold
            nEeff = Eeff[bool_Eeff]
            sK = ((self.KE.flatten()[np.newaxis]).T * nEeff.ravel()).flatten(order='F')

            one = np.ones((self.nodesPerElement*self.dofpernode, 1), dtype=np.int_)
            local_iK = np.kron(self.edofMat[bool_Eeff,:],one).flatten()
            one.shape = (1,self.nodesPerElement*self.dofpernode)
            local_jK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()

            K = coo_matrix((sK, (local_iK, local_jK)), shape=(self.ndof, self.ndof)).tocsr()#(self.dofpernode,self.dofpernode))

        zerosdof = np.where(K.diagonal()== 0 )[0]

        self.PrintVerbose("Number of active nodes : " + str(self.ndof-len(zerosdof) ) + "  of " + str(self.ndof) + "   "+ str(float(len(zerosdof)*100.)/self.ndof)+ "% of empty dofs"  )
        Kones = coo_matrix( (np.ones((len(zerosdof),) ) ,(zerosdof,zerosdof)), shape =(self.ndof, self.ndof)).tocsr()#(self.dofpernode,self.dofpernode))
        K = (K.tocsr() + Kones.tocsr()).tocsr()

        # Remove constrained dofs from matrix
        self.PrintDebug(" Delete fixed Dofs")
        [K, rhsfixed] = FeaBase.deleterowcol(K, self.fixed, self.fixed, self.fixedValues)

        self.PrintDebug(" Start solver")
        rhs = self.f[self.free, 0]-rhsfixed[self.free, 0]

        self.u = np.zeros((self.ndof, 1), dtype=np.double)

        if K.nnz > 0 :
          if self.linearSolver == "Direct":
            self.u[self.free, 0] = linalg.spsolve(K, self.f[self.free, 0]-rhsfixed[self.free, 0])
          elif self.linearSolver == "DirectDense":
            self.u[self.free, 0] = denselinalg.solve(K.toarray(), self.f[self.free, 0]-rhsfixed[self.free, 0],sym_pos=True,overwrite_a=True)
          elif self.linearSolver == "CG":
            # Preconditioning
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            # normalization of the rhs term ( to treat correctly the tol of CG) (Please read the documentaion of numpy.linalg.cg)
            norm = np.linalg.norm(rhs)
            res = linalg.cg(K, rhs/norm, x0 = self.u[self.free, 0]/norm , M = M, tol = self.tol)
            self.u[self.free, 0] = res[0]*norm
          elif self.linearSolver == "LGMRES":
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            norm = np.linalg.norm(rhs)
            res = linalg.lgmres(K, rhs/norm, x0 = self.u[self.free, 0]/norm , M = M  )
            self.u[self.free, 0] = res[0]*norm
          elif self.linearSolver == "AMG":
            try:
                import pyamg
            except:#pragma: no cover
                raise Exception('AMG module not installed')#pragma: no cover
            K = K.tocsr()
            ml = pyamg.ruge_stuben_solver(K)
            res = ml.solve(rhs,x0 = self.u[self.free, 0] , tol=1e-12,accel='cg')
            self.u[self.free, 0] = res
          else :
            print("'"+self.linearSolver + "' is not a valid linear solver")#pragma: no cover
            print('Please set a type of linear solver')#pragma: no cover
            raise Exception()#pragma: no cover

        self.PrintDebug('Post Process')
        self.u = self.u + self.fixedValues

        # Compute element elastic energy density
        u_reshaped = self.u[self.edofMat]
        u_reshaped.shape = (self.support.GetNumberOfElements(), self.nodesPerElement*self.dofpernode)
        Ku_reshaped = np.dot(u_reshaped, self.KE)
        np.einsum('ij,ij->i', Ku_reshaped, u_reshaped, out=self.eed)

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
            Eeff = np.ones(self.support.GetNumberOfElements())

        bool_Eeff = Eeff<self.minthreshold
        nEeff = Eeff.copy()
        nEeff[bool_Eeff] = 0.0
        nEeff *= self.eed

        return nEeff.ravel()

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

    from OTTools.FE.Hexa8Cuboid import Hexa8Cuboid
    import OTTools.FE.ConstantRectilinearMesh as CRM
    import OTTools.IO.XdmfWriter  as XdmfWriter
    import time
    from OTTools.Helpers.Tests import TestTempDir

    myMesh = CRM.ConstantRectilinearMesh()
    nx = 11; ny = 12; nz = 13;


    myMesh.SetDimensions([nx,ny,nz]);
    myMesh.SetSpacing([0.1, 0.1, 10/(nz-1)]);
    myMesh.SetOrigin([0, 0, 0]);
    print(myMesh)

    # thermal problem
    #dirichlet at plane z =0
    dirichlet_bcs = BundaryCondition()

    for x in xrange(nx):
        for y in xrange(ny):
            for z in [0]:
                dirichlet_bcs.append([x,y,z],0 , 0 )
                

    # Homogenous body flux
    neumann_bcs = BundaryCondition()
    for x in xrange(nx):
        for y in xrange(ny):
            for z in xrange(nz):
                neumann_bcs.append([x,y,z],0 , 1 )

    neumann_bcs.append([0,0,0],0 , 1 )
    neumann_bcs.append([1,1,0],0 , 1 )
    neumann_bcs.eliminate_double()

    starttime = time.time()
    myElement = Hexa8Cuboid()
    myElement.delta =  myMesh.GetSpacing()

    myProblem = Fea(myMesh, dofpernode = 1,
        KOperator = myElement.GetIsotropLaplaceK(1.),
        MOperator = myElement.GetIsotropLaplaceM(1.),
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
    try:
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
    import OTTools.FE.ConstantRectilinearMesh as CRM
    import OTTools.IO.XdmfWriter  as XdmfWriter

    from OTTools.Helpers.Tests import TestTempDir

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

    myProblem = Fea(myMesh, dofpernode = 3,
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
        raise   # pragma: no cover

def CheckIntegrityThermal2D():
    import OTTools.FE.ConstantRectilinearMesh as CRM
    from OTTools.FE.Quad4Rectangle import Quad4Rectangle
    import OTTools.IO.XdmfWriter  as XdmfWriter
    import time
    from OTTools.Helpers.Tests import TestTempDir

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
    myElement = Quad4Rectangle()
    myElement.delta =  myMesh.GetSpacing()

    myProblem = Fea(myMesh, dofpernode = 1,
        KOperator = myElement.GetIsotropLaplaceK(1.),
        MOperator = myElement.GetIsotropLaplaceM(1.),
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
    try:
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

    import OTTools.FE.ConstantRectilinearMesh as CRM
    import OTTools.IO.XdmfWriter  as XdmfWriter
    import time
    from OTTools.Helpers.Tests import TestTempDir

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

    myProblem = Fea(myMesh, dofpernode = 2,
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
    myProblem.BuildMassMatrix(densities);

    fixed = np.zeros((myProblem.ndof, 1), dtype=np.double)

    for i in myProblem.fixed:
        fixed[i] =  1;

    XdmfWriter.WriteMeshToXdmf(TestTempDir.GetTempPath() +'TestDep2D.xmf',myMesh,
                    [myProblem.u, myProblem.f, myProblem.nodal_elastic_energy(), fixed ],
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
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
