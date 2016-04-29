# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np

from scipy.sparse import coo_matrix
import scipy.sparse.linalg as linalg
import scipy.linalg as denselinalg

import  scipy.sparse as sps
from OTTools.FE.Hexa8Cuboid import Hexa8Cuboid
from OTTools.Helpers.BaseOutputObject import BaseOutputObject as BOO


class Fea(BOO):
    
    def __init__(self, theMeshObj, dofpernode = 1, dirichlet_bcs = [], neumann_bcs= [], KOperator= None, MOperator= None,  neumann_nodal= None):
        super(Fea,self).__init__()
        
        self.linearSolver = "CG"
        self.outer_v = []
        
        self.theMeshObj = theMeshObj ;
        
        self.dofpernode = dofpernode
        
        self.writer = None
        
        if KOperator is not None:
            self.KE    = KOperator            
            self.ME    = MOperator            
        else:
            myElem = Hexa8Cuboid()
            myElem.delta = theMeshObj.GetSpacing()
            self.KE = myElem.GetIsotropDispK(1.,0.3);
            self.ME = myElem.GetIsotropDispM(1,);
        
        # dofs:
        self.ndof = dofpernode * theMeshObj.GetNumberOfNodes()
        

        # FE: Build the index vectors for the for coo matrix format
        self.edofMat = np.zeros((theMeshObj.GetNumberOfElements(), 8*dofpernode), dtype=np.int_)
        for i in  range(theMeshObj.GetNumberOfElements()):
            coon = theMeshObj.GetConnectivityForElement(i)
            self.edofMat[i, :] = np.array([(x*dofpernode+y) for x in coon for y in xrange(dofpernode)])
#            self.edofMat[i, :] = np.array([coon[0]*3+0, coon[0]*3+1, coon[0]*3+2, 
#                                           coon[1]*3+0, coon[1]*3+1, coon[1]*3+2,
#                                           coon[2]*3+0, coon[2]*3+1, coon[2]*3+2,
#                                           coon[3]*3+0, coon[3]*3+1, coon[3]*3+2,
#                                           coon[4]*3+0, coon[4]*3+1, coon[4]*3+2,
#                                           coon[5]*3+0, coon[5]*3+1, coon[5]*3+2,
#                                           coon[6]*3+0, coon[6]*3+1, coon[6]*3+2,
#                                           coon[7]*3+0, coon[7]*3+1, coon[7]*3+2])
#            
        # Construct the index pointers for the coo format
        self.iK = np.kron(self.edofMat, np.ones((8*dofpernode, 1), dtype=np.int_)).flatten()
        self.jK = np.kron(self.edofMat, np.ones((1, 8*dofpernode), dtype=np.int_)).flatten()

        # BC's and support
        self.dofs = np.arange(self.ndof)
        self.fixed = np.array([dofpernode * theMeshObj.GetMonoIndexOfNode(d[0:3])  + d[3] for d in dirichlet_bcs])
        
        self.fixedValues = np.zeros((self.ndof, 1), dtype=np.double)
        
        if len(dirichlet_bcs) and len(dirichlet_bcs[0])>4:
            vals = np.array([[ d[4] for d in dirichlet_bcs]]).T
            self.fixedValues[self.fixed.T,0:] = vals

        self.free = np.setdiff1d(self.dofs, self.fixed)
        # Solution and RHS vectors
        self.f = np.zeros((self.ndof, 1), dtype=np.double)
        self.u = np.zeros((self.ndof, 1), dtype=np.double)


        if (neumann_bcs is not None) or  (neumann_nodal is not None):
            MassMatrix = self.BuildMassMatrix()    
            # Set load
            if  neumann_bcs is not  None:
                for n in neumann_bcs:
                    self.f[dofpernode * theMeshObj.GetMonoIndexOfNode(n[0:3]) + n[3]] = n[4]
                
                self.f[:,0] = MassMatrix*self.f[:,0]
        
            if  neumann_nodal is not  None:
                nodal_f = np.zeros((self.ndof, 1), dtype=np.double)
                for n in neumann_nodal:
                    nodal_f[dofpernode * theMeshObj.GetMonoIndexOfNode(n[0:3]) + n[3]] = n[4]
        
                dd = MassMatrix.diagonal()
                self.f[:,0] += np.multiply(dd,nodal_f[:,0])
                
        # Element elastic energy density
        # (equivalent to the elastic energy for an element fully in the structure)
        self.eed = np.zeros(theMeshObj.GetNumberOfElements())
        
        self.last_element_elastic_energy = None
        self.last_element_Eeff = None
        
    def BuildMassMatrix(self, Eeff = None):
        if Eeff is None:
            Eeff = np.ones(self.theMeshObj.GetNumberOfElements())
            
        sM = ((self.ME.flatten()[np.newaxis]).T * Eeff.ravel()).flatten(order='F')
        M = coo_matrix((sM, (self.iK, self.jK)), shape=(self.ndof, self.ndof)).tocsr()
        return M
    

    def Solve(self, Eeff=None):
#        starttime = time.time()
        if Eeff is None:
            Eeff = np.ones(self.theMeshObj.GetNumberOfElements())
#        print("Time for Eeff definition : " + str(time.time() -starttime))
#        starttime = time.time()
        
        # Setup and solve FE problem
        
        self.PrintDebug("Construction of the tangent matrix")
        sK = ((self.KE.flatten()[np.newaxis]).T * Eeff.ravel()).flatten(order='F')
        #K = coo_matrix((sK, (self.iK, self.jK)), shape=(self.ndof, self.ndof)).tocsc()
        K = coo_matrix((sK, (self.iK, self.jK)), shape=(self.ndof, self.ndof)).tocsr()
        zerosdof = np.where(K.diagonal()== 0 )[0]
        Kones = coo_matrix( (np.ones((len(zerosdof),) ) ,(zerosdof,zerosdof)), shape =(self.ndof, self.ndof)).tocsr()
        K = (K + Kones).tocsr()

        # Remove constrained dofs from matrix
        self.PrintDebug(" Delete fixed Dofs")
        [K, rhsfixed] = deleterowcol(K, self.fixed, self.fixed, self.fixedValues)
        
        self.PrintDebug(" Start solver")
        rhs = self.f[self.free, 0]-rhsfixed[self.free, 0]
        if self.linearSolver == "Direct":
            self.u[self.free, 0] = linalg.spsolve(K, self.f[self.free, 0]-rhsfixed[self.free, 0])            
        elif self.linearSolver == "DirectDense":
            self.u[self.free, 0] = denselinalg.solve(K.toarray(), self.f[self.free, 0]-rhsfixed[self.free, 0],sym_pos=True,overwrite_a=True)            
        elif self.linearSolver == "CG":
            M = sps.dia_matrix((1./K.diagonal(),0), shape=K.shape)
            res = linalg.cg(K, rhs, x0 = self.u[self.free, 0] , M = M)#,tol=  1e-4) 
            self.u[self.free, 0] = res[0]
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
        u_reshaped.shape = (self.theMeshObj.GetNumberOfElements(), 8*self.dofpernode)
        Ku_reshaped = np.dot(u_reshaped, self.KE)
        np.einsum('ij,ij->i', Ku_reshaped, u_reshaped, out=self.eed)
        self.PrintDebug('Post Process Done')
        
    def element_elastic_energy(self, Eeff= None):
        
        if Eeff is None:
            Eeff = np.ones(self.theMeshObj.GetNumberOfElements())
            
        result = np.ravel(Eeff) * self.eed
        
        self.last_element_elastic_energy = result
        self.last_element_Eeff = Eeff
        
        return result

    def nodal_elastic_energy(self, Eeff=None):
        
        if Eeff is None:
            Eeff = np.ones(self.theMeshObj.GetNumberOfElements())
            
        return node_averaged_element_field(self.element_elastic_energy(Eeff),self.theMeshObj)

    def Write(self):
        
        if self.writer is not None:
          if self.last_element_elastic_energy is not None:
            self.writer.Write(self.theMeshObj,
                PointFields     = [self.u, self.f],
                PointFieldsNames= ["u", "f"],
                CellFields     = [self.last_element_elastic_energy, self.last_element_Eeff],
                CellFieldsNames= ["elastic_enery", "Eeff"]
                )
          else:
            self.writer.Write(self.theMeshObj,
                PointFields     = [self.u, self.f],
                PointFieldsNames= ["u", "f"]
                )
                
                
def deleterowcol(A, delrow, delcol, fixedValues ):
    # Assumes that matrix is in symmetric csc form !
    m = A.shape[0]
   
    rhs = A*fixedValues
    keep = np.delete (np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete (np.arange(0, m), delcol)
    A = A[:, keep]
    
    return [A, rhs]

def node_averaged_element_field(element_field,support):
    nnodes = support.GetDimensions()
    
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
    
    

    
def CheckIntegrity():
    
    from OTTools.FE.Hexa8Cuboid import Hexa8Cuboid
    import OTTools.FE.ConstantRectilinearMesh as CRM
    import OTTools.IO.XdmfWriter  as XdmfWriter 
    import time
    from OTTools.Helpers.Tests import TestTempDir

    myMesh = CRM.ConstantRectilinearMesh()
    nx = 11; ny = 11; nz = 101;
    
    
    myMesh.SetDimensions([nx,ny,nz]);
    myMesh.SetSpacing([0.1, 0.1, 0.1]);
    myMesh.SetOrigin([0, 0, 0]);
    print(myMesh)
    
    # thermal problem 
    #dirichlet at plane z =0
    dirichlet_bcs =( [(x, y, z, 0) for x in range(ny) for y in range(ny) for z in [0] ]  )

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
        
    
    #print('-----------------------------------------------------------------------------')
    
    # block all the faces rith
    
    dirichlet_bcs =( [(x, y, z, coor, 0) for x in range(ny) for y in range(ny) for z in [0]  for coor in range(3)] +
                    [(x, y, z, coor, 1) for x in range(ny) for y in range(ny) for z in [nz-1]  for coor in range(3)])

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
        
    return 'OK'
    #dirichlet_bcs =( [(0, y, z, cor) for y in range(ny) for z in range(nz) for cor in range(3)] )
    #neumann_bcs = ([(nx-1, y, z, 2) for y in range(ny) for z in range(nz) ])
    #Fea(myMesh, dofpernode = 1, Operator= None, dirichlet_bcs = dirichlet_bcs, neumann_bcs =neumann_bcs):

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover