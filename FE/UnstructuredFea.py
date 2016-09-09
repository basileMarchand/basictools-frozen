# -*- coding: utf-8 -*-



from __future__ import division
import numpy as np

from scipy.sparse import coo_matrix
import scipy.sparse.linalg as linalg
import scipy.linalg as denselinalg
import  scipy.sparse as sps

from OTTools.FE.FeaBase import FeaBase as FeaBase

import OTTools.FE.ElementBuildier  as ElementBuildier
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
        self.nodes.resize( (size,self.dim))
        self.dofs.resize( ( size,))
        self.vals.resize( ( size,1))
        self.sz = size

    def tighten(self):
        self.reserve(self.cpt)

    def append(self, nodes, dof,val):
        if self.cpt >= self.sz:
            self.reserve(self.sz*2)
        self.nodes[self.cpt,:] = nodes
        self.dofs[self.cpt] = dof
        self.vals[self.cpt] = val
        self.cpt += 1


class Fea(FeaBase):
    def __init__(self, support, dofpernode = 1, dirichlet_bcs = None, neumann_bcs= None,  neumann_nodal= None):
        super(Fea,self).__init__()

        self.linearSolver = "CG"
        #self.outer_v = []
        self.support = support ;

        self.dofpernode = dofpernode

        self.writer = None

        self.minthreshold = 1.e-3

        self.tol = 1.e-6

        # dofs:
        self.ndof = dofpernode * support.GetNumberOfNodes()

        self.fixedValues = np.zeros((self.ndof, 1), dtype=np.double)

        self.PrintDebug("Treating Dirichlet 1/4")

        self.fixed = np.zeros(self.ndof, dtype=np.bool)
        if dirichlet_bcs is not None :
            dirichlet_bcs.tighten()
            print(dirichlet_bcs.nodes)
            indexs = support.GetMonoIndexOfNode(dirichlet_bcs.nodes)
            indexs *= dofpernode
            indexs += dirichlet_bcs.dofs
            self.fixed[indexs] = True
            self.fixedValues[self.fixed.T,0:] = dirichlet_bcs.vals
        self.PrintDebug("Treating Dirichlet 2/4")

        self.free = np.ones(self.ndof, dtype=np.bool)
        self.free[self.fixed] = False

        self.PrintDebug("Treating Dirichlet 3/4")
        self.f = np.zeros((self.ndof, 1), dtype=np.double)
        self.PrintDebug("Treating Dirichlet 4/4")
        self.u = np.zeros((self.ndof, 1), dtype=np.double)
        self.PrintDebug("Treating Neumann")

#        if (neumann_bcs is not None) or  (neumann_nodal is not None):
#            if  neumann_bcs is not  None:
#                MassMatrix = self.BuildMassMatrix()
#                self.f[support.GetMonoIndexOfNode(neumann_bcs.nodes)*dofpernode + neumann_bcs.dofs] += neumann_bcs.vals
#                #for n in neumann_bcs:
#                #    self.f[dofpernode * support.GetMonoIndexOfNode(n[0:ndims]) + n[ndims]] += n[ndims+1]
#                self.PrintDebug(np.linalg.norm(self.f[:,0]))
#                self.PrintDebug(MassMatrix)
#                self.f[:,0] = MassMatrix*self.f[:,0]
#                self.PrintDebug(np.linalg.norm(self.f[:,0]))
#
#            if  neumann_nodal is not  None:
#                nodal_f = np.zeros((self.ndof, 1), dtype=np.double)
#
#                #self.PrintDebug(neumann_nodal.nodes)
#                #self.PrintDebug(neumann_nodal.dofs)
#                nodal_f[support.GetMonoIndexOfNode(neumann_nodal.nodes)*dofpernode + neumann_nodal.dofs] += neumann_nodal.vals
#                #for n in neumann_nodal:
#                #    nodal_f[dofpernode * support.GetMonoIndexOfNode(n[0:ndims]) + n[ndims]] += n[ndims+1]
#
#                #dd = MassMatrix.diagonal()
#                #self.f[:,0] += np.multiply(dd,nodal_f[:,0])
#                self.f[:,0] += nodal_f[:,0]
#
#        self.PrintDebug("Treating Neumann Done")
        # Element elastic energy density
        # (equivalent to the elastic energy for an element fully in the structure)
        self.eed = np.zeros(support.GetNumberOfElements())


#    def BuildMassMatrix(self, Eeff = None):

#        self.PrintDebug("BuildMassMatrix")
#        if Eeff is None:
#            self.PrintDebug(" Eeff is None")
#
#            Eeff = np.ones(self.support.GetNumberOfElements())
#
#            sM = ((self.ME.flatten()[np.newaxis]).T * Eeff.ravel()).flatten(order='F')
#            #print(sM[0:10])
#            #sM = np.tile(self.ME.flatten(),self.support.GetNumberOfElements())
#            #print(sM[0:10])
#            self.GenerateIJs()
#            self.PrintDebug("Asm")
#            M = coo_matrix((sM, (self.iK, self.jK)), shape=(self.ndof, self.ndof),dtype=float)
#            self.PrintDebug("BuildMassMatrix Done")
#            return  M.tocsr()#(self.dofpernode,self.dofpernode))
#        else:
#            bool_Eeff = Eeff>self.minthreshold
#            nEeff = Eeff[bool_Eeff]
#            sM = ((self.ME.flatten()[np.newaxis]).T * nEeff.ravel()).flatten(order='F')

#            nodesPerElement = 2**self.support.GetDimensionality()

#            one = np.ones((nodesPerElement*self.dofpernode, 1), dtype=np.int_)
#            local_iK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()
#            one.shape = (1,nodesPerElement*self.dofpernode)
#            local_jK = np.kron(self.edofMat[bool_Eeff,:], one).flatten()
#            self.PrintDebug("Asm")
#            M = coo_matrix((sM, (local_iK, local_jK)), shape=(self.ndof, self.ndof),dtype=float).tocsr()#(self.dofpernode,self.dofpernode))

#            self.PrintDebug("BuildMassMatrix Done")
#            return M


    def IntegrateRight(self,weak, tag=None):
        pass

    def IntegrateLeft(self,weak, tag=None):
        #the filter send a type

        #first pass to count the number of element to integrate
        #number of elements to treat
        nelems = 0
        #number of values
        nvalues = 0
        pos = self.support.GetPosOfNodes()

        if tag is None:
            for ntype, elem in self.support.elements.iteritems():
               nelems +=elem.GetNumberOfElements()
               nvalues += elem.GetNumberOfNodesPerElement()*self.dofpernode

        else:
           for ntype, elem in self.support.elements.iteritems():
              if elem.tags.has_key(tag):
                  nelems += elem.GetNumberOfElements()
                  nvalues += elem.GetNumberOfElements()*((elem.GetNumberOfNodesPerElement()*self.dofpernode)**2 )

           if nvalues == 0:
               self.Print("Nothing to integrate")
               return

           sK = np.empty((nvalues,))
           iK = np.empty((nvalues,))
           jK = np.empty((nvalues,))

           cpt = 0
           for ntype, elem in self.support.elements.iteritems():
               if elem.tags.has_key(tag):
                 t = elem.GetTag(tag)
                 E = ElementBuildier.GetElementFromName(ntype)
                 conectivities = elem.connectivity[t.id,:]
                 nelemintag = len(t)
                 nodesPerElement = elem.GetNumberOfNodesPerElement()
                 local_nvals = (nodesPerElement*self.dofpernode)**2
                 self.Print("Building matrix for elements of type : " + ntype)
                 for i in xrange(nelemintag):
                     coon = conectivities[i,:]
                     coords = pos[coon,:]
                     edofMat= np.array([(coon*self.dofpernode+y) for y in xrange(self.dofpernode)]).flatten('F')
                     local_iK = np.kron(edofMat, np.ones((nodesPerElement*self.dofpernode, 1), dtype=np.int_)).flatten()
                     local_jK = np.kron(edofMat, np.ones((1,nodesPerElement*self.dofpernode), dtype=np.int_)).flatten()
                     iK[cpt:cpt+local_nvals] = local_iK
                     jK[cpt:cpt+local_nvals] = local_jK

                     vals = weak(E, coords).flatten()
                     #self.Print("vals")
                     #self.Print(vals)
                     sK[cpt+0:cpt+local_nvals] =  vals
                     cpt +=local_nvals

        K = coo_matrix((sK, (iK,jK)), shape=(self.ndof, self.ndof)).tocsr()#(self.dofpernode,self.dofpernode))
        return K

    def Solve(self, CleanZerosLines = True):

        print(self.K.diagonal())
        if CleanZerosLines and len(np.where( self.K.diagonal() == 0 )[0]) > 0:
            zerosdof = np.where(self.K.diagonal()== 0 )[0]
            self.PrintDebug("Number of active nodes : " + str(self.ndof-len(zerosdof) ) + "  of " + str(self.ndof) + "   "+ str(float(len(zerosdof)*100.)/self.ndof)+ "% of empty nodes"  )
            Kones = coo_matrix( (np.ones((len(zerosdof),) ) ,(zerosdof,zerosdof)), shape =(self.ndof, self.ndof)).tocsr()#(self.dofpernode,self.dofpernode))
            K = (self.K.tocsr() + Kones.tocsr())

            self.PrintDebug(" Delete fixed Dofs")
            [K, rhsfixed] = deleterowcol(K, self.fixed, self.fixed, self.fixedValues)
        else:
            self.PrintDebug(" Delete fixed Dofs")
            [K, rhsfixed] = deleterowcol(self.K, self.fixed, self.fixed, self.fixedValues)

        self.PrintDebug(" Start solver")
        rhs = self.f[self.free, 0]-rhsfixed[self.free, 0]
        self.u = np.zeros((self.ndof, 1), dtype=np.double)

        if K.nnz > 0 :
          if self.linearSolver == "Direct":
            self.u[self.free, 0] = linalg.spsolve(K, self.f[self.free, 0]-rhsfixed[self.free, 0])
          elif self.linearSolver == "DirectDense":
            self.u[self.free, 0] = denselinalg.solve(K.toarray(), self.f[self.free, 0]-rhsfixed[self.free, 0],sym_pos=True,overwrite_a=True)
          elif self.linearSolver == "CG":
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            norm = np.linalg.norm(rhs)
            res = linalg.cg(K, rhs/norm, x0 = self.u[self.free, 0]/norm , M = M, tol = self.tol)
            self.u[self.free, 0] = res[0]*norm
          elif self.linearSolver == "LGMRES":
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            res = linalg.lgmres(K, rhs, x0 = self.u[self.free, 0] , M = M  )
            #, outer_k=1, store_outer_Av= True,outer_v = self.outer_v
            #print(self.outer_v)
            self.u[self.free, 0] = res[0]
          elif self.linearSolver == "AMG":
            try:
                import pyamg
            except:#pragma: no cover
                raise Exception('AMG module not installed')#pragma: no cover
            K = K.tocsr()
            ml = pyamg.ruge_stuben_solver(K)
            res = ml.solve(rhs,x0 = self.u[self.free, 0] , tol=1e-6,accel='cg')
            self.u[self.free, 0] = res
          else :
            print("'"+self.linearSolver + "' is not a valid linear solver")#pragma: no cover
            print('Please set a type of linear solver')#pragma: no cover
            raise Exception()#pragma: no cover


        self.PrintDebug('Post Process')
        self.u = self.u + self.fixedValues
        self.PrintDebug('Post Process Done')



    def Write(self):

        if self.writer is not None:
          if self.last_element_elastic_energy is not None:
            self.writer.Write(self.support,
                PointFields     = [self.u, self.f],
                PointFieldsNames= ["u", "f"],
                CellFields     = [self.last_element_elastic_energy],
#, self.last_element_Eeff
                CellFieldsNames= ["elastic_enery"]
                #, "Eeff",
                )
          else:
            self.writer.Write(self.support,
                PointFields     = [self.u, self.f],
                PointFieldsNames= ["u", "f"]
                )


def deleterowcol(A, delrow, delcol, fixedValues ):
    # Assumes that matrix is in symmetric csc form !
    m = A.shape[0]

    rhs = A*fixedValues
    #keep = np.delete (np.arange(0, m), delrow)
    A = A[np.logical_not(delrow) , :]
    #keep = np.delete (np.arange(0, m), delcol)
    A = A[:, np.logical_not(delcol)]

    return [A, rhs]



def CheckIntegrityold():

    import OTTools.IO.GmshReader as GR
    import OTTools.TestData as test
    import OTTools.IO.XdmfWriter  as XdmfWriter
    from OTTools.Helpers.Tests import TestTempDir
    import time

    myMesh = GR.ReadGmsh( test.GetTestDataPath() + 'mesh1.msh')
    print(myMesh)


    def MyWeak(Elem,pos):
        return Elem.GetIsotropDispK(1,0.5,pos)

    myProblem = Fea(myMesh, dofpernode = 2)
    myProblem.SetGlobalDebugMode()
    myProblem.K = myProblem.IntegrateLeft(MyWeak,tag="GeoTag10")


    return
    # thermal problem
    #dirichlet at plane z =0
    dirichlet_bcs =( [(x, y, z, 0) for x in range(nx) for y in range(ny) for z in [0] ]  )

    # Homogenous body flux
    neumann_bcs = ([(x, y, z, 0, 1) for x in range(nx)  for y in range(ny) for z in range(nz) ])

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
       myProblem.linearSolver = 'AMG'
       myProblem.Solve()
       myProblem.Write()
    except Exception as inst: # pragma: no cover
       print(inst.args[0])# pragma: no cover

    myProblem.writer.Close()


    path = TestTempDir.GetTempPath() +'TestThermal.xmf'
    XdmfWriter.WriteMeshToXdmf(path,myMesh,
                    [myProblem.u, myProblem.f],
                    PointFieldsNames= ['Themperature','q'],
                    GridFieldsNames=[])
    print('DONE')
    print(max(myProblem.u))

    if abs(max(myProblem.u)-50.) > 1e-5:
        raise Exception()# pragma: no cover

    print('-----------------------------------------------------------------------------')

    # block all the faces rith

    dirichlet_bcs =( [(x, y, z, coor, 0) for x in range(nx) for y in range(ny) for z in [0]  for coor in range(3)] +
                    [(x, y, z, coor, 1) for x in range(nx) for y in range(ny) for z in [nz-1]  for coor in range(3)])

    neumann_bcs = ((0,0,0,0,0),)
    neumann_nodal = ((0,0,0,0,0),)
    starttime = time.time()

    #myElement = Hexa8Cuboid()
    #myElement.delta =  myMesh.GetSpacing()

    myProblem = Fea(myMesh, dofpernode = 3,
        dirichlet_bcs = dirichlet_bcs,
        neumann_bcs = neumann_bcs,
        neumann_nodal = neumann_nodal)

    print("Time for Fea definition : " + str(time.time() -starttime))

    starttime = time.time()
    densities  = np.ones(myMesh.GetNumberOfElements());
    myProblem.Solve(densities)
    print("Time for Fea solve : " + str(time.time() -starttime))

    XdmfWriter.WriteMeshToXdmf(TestTempDir.GetTempPath() +'TestDep.xmf',myMesh,
                    [myProblem.u, myProblem.f, myProblem.nodal_elastic_energy() ],
                    [densities, myProblem.element_elastic_energy() ] ,
                    [],
                    PointFieldsNames= ['Dep','F','NEnergie'],
                    CellFieldsNames=['densities','EEnergie'],
                    GridFieldsNames=[])

    #print(max(myProblem.u))

    if abs(max(myProblem.u)-1.0128810548) > 1e-5:
        raise   # pragma: no cover

    print('----------------------------  2D -------------------------------------------------')
    del nz
    myMesh = CRM.ConstantRectilinearMesh(2)
    nx = 11; ny = 11;


    myMesh.SetDimensions([nx,ny]);
    myMesh.SetSpacing([0.1, 0.1]);
    myMesh.SetOrigin([0, 0]);
    print(myMesh)

    # thermal problem
    #dirichlet at plane z =0
    dirichlet_bcs =( [(x, y,  0) for x in range(nx) for y in [0] ]  )

    # Homogenous body flux
    neumann_bcs = ([(x, y, 0, 1) for x in range(nx)  for y in range(ny) ])

    starttime = time.time()
    myElement = Quad4Rectangle()
    myElement.delta =  myMesh.GetSpacing()

    myProblem = Fea(myMesh, dofpernode = 1,
        KOperator = myElement.GetIsotropLaplaceK(1.),
        MOperator = myElement.GetIsotropLaplaceM(1.),
        dirichlet_bcs = dirichlet_bcs,
        neumann_bcs = neumann_bcs)

    print("Time for Fea definition : " + str(time.time() -starttime))


    myProblem.writer = XdmfWriter.XdmfWriter(TestTempDir.GetTempPath()+"Test_constantRectilinearThermal2D.xdmf")
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
    myMesh.SetSpacing([0.1, 1./(ny-1)]);
    myMesh.SetOrigin([0, 0]);
    print(myMesh)
    # block all the faces rith
    dirichlet_bcs = BundaryCondition(dim=2);

    for x in range(nx):
        for y in [0]:
            for coor in range(2):
               dirichlet_bcs.append([x,y],coor , 0 )
        for y in [ny-1]:
            for coor in range(2):
               dirichlet_bcs.append([x,y],coor , 1 )

    print(dirichlet_bcs)
    #dirichlet_bcs =( [(x, y, coor, 0) for x in range(nx) for y in [0]  for coor in range(2)] +
    #                [(x, y, coor, 1) for x in range(nx) for y in  [ny-1]  for coor in range(2)])
    #print(dirichlet_bcs )

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
    print()
    print(myProblem.u.T)

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

    #print(max(myProblem.u))
    import os
    os.system('cmd /C C:\Users\D584808\Apps\ParaView-5.0.1-Qt4-OpenGL2-Windows-64bit\\bin\paraview.exe ' + TestTempDir.GetTempPath() +'TestDep2D.xmf')


    return 'ok'
    #if abs(max(myProblem.u)-1.0128810548) > 1e-5:
    #    raise   # pragma: no cover

def CheckIntegrity():
    #try:
        #print(CheckIntegrityThermal3D())
        #print(CheckIntegrityDep3D())
        #print(CheckIntegrityThermal2D())
        print(CheckIntegrityDep2D())
    #except:# pragma: no cover
    #    return "Not ok"# pragma: no cover
    #return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())#pragma: no cover
