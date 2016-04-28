# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import MeshBase
from OTTools.FE.Hexa8Cuboid import Hexa8Cuboid

import numpy as np
from scipy.sparse import coo_matrix

class ConstantRectilinearMesh(MeshBase):
    
    def IsConstantRectilinear(self): 
        return True
        
    def __init__(self):
        super(ConstantRectilinearMesh,self).__init__()
        #Number of nodes
        self.__dimensions = np.array([2, 2, 2]);
        self.__origin = np.array([0, 0, 0])
        self.__spacing = np.array([1, 1, 1])
        self.elemTags = {}
        
    def GetNamesOfCellTags(self):
        return self.elemTags.keys()
    
    def GetSubSuperMesh(self,_newDimensions):
        newDimensions = np.array(_newDimensions)
        ## to generate meshes with more or less elements in each directions
        ## return the mesh and two operators: from old to new and from new to old
        
            
        newSpac = (self.__dimensions-1)*self.__spacing/(newDimensions-1)
        
        res = type(self)()
        res.SetSpacing(newSpac)
        res.SetDimensions(newDimensions)
        res.SetOrigin(self.GetOrigin())
        
        return res
        
    def GetNodeTrasfertMatrix(self, destination):
        # newVector   = oldToNew * oldVector
        oldToNewVals = np.zeros((destination.GetNumberOfNodes(),8))
        oldToNewIK = np.zeros((destination.GetNumberOfNodes(),8), dtype=np.int_)
        oldToNewJK = np.zeros((destination.GetNumberOfNodes(),8), dtype=np.int_)
        
        for i in  xrange(destination.GetNumberOfNodes()):
            
            pos= destination.GetPosOfNode(i)
            el = self.GetElementAtPos(pos)
            coon = self.GetConnectivityForElement(el)
            xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)
            oldToNewVals[i,:] = xiChiEta
            oldToNewIK[i,:] = i
            oldToNewJK[i,:] = coon
            
        oldToNew =  coo_matrix((oldToNewVals.flatten(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destination.GetNumberOfNodes(), self.GetNumberOfNodes())).tocsc()

        return oldToNew
        
    def GetElementTrasfertMatrix(self, destination):
        
        nps = 10
        nps3 = nps**3
        oldToNewVals = np.zeros((destination.GetNumberOfNodes(),nps3))
        oldToNewIK = np.zeros((destination.GetNumberOfNodes(),nps3), dtype=np.int_)
        oldToNewJK = np.zeros((destination.GetNumberOfNodes(),nps3), dtype=np.int_)
                
        for i in  xrange(destination.GetNumberOfElements()):
            coon = destination.GetConnectivityForElement(i)
            
            n0pos = destination.GetPosOfNode(coon[0])
            cpt =0
            for cx in range(0,nps):
                for cy in range(0,nps):
                    for cz in range(0,nps):
                        pos = n0pos + destination.GetSpacing()*([cx+0.5,cy+0.5,cz+0.5])/nps
                        el = self.GetElementAtPos(pos)
                        oldToNewVals[i,cpt] += 1./nps3
                        oldToNewIK[i,cpt] += i
                        oldToNewJK[i,cpt] += el
                        cpt +=1
                        
            #pos = destination.GetPosOfNode(coon[0]) + destination.GetSpacing()/2
            #el = self.GetElementAtPos(pos)
            #oldToNewVals[i,:] = 1
            #oldToNewIK[i,:] = i
            #oldToNewJK[i,:] = el
        oldToNew =  coo_matrix((oldToNewVals.flatten(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destination.GetNumberOfElements(), self.GetNumberOfElements())).tocsc()

        return oldToNew            
    def SetDimensions(self,data):
        self.__dimensions = np.array(data);
   
    def GetDimensions(self):
        return self.__dimensions ;
        
    def SetSpacing(self,data):
        self.__spacing = np.array(data, "float");
        
    def GetSpacing(self):
        return self.__spacing;

    def GetOrigin(self):
        return self.__origin;        
        
    def SetOrigin(self,data):
        self.__origin = np.array(data);        
      
    def GetNumberOfNodes(self):
        return self.__dimensions[0]*self.__dimensions[1]*self.__dimensions[2]
        
    def GetNumberOfElements(self):
        res = 1;
        if self.__dimensions[0] >= 1:
            res = res * (self.__dimensions[0]-1)

        if self.__dimensions[1] >= 1:
            res = res * (self.__dimensions[1]-1)

        if self.__dimensions[0] >= 1:
            res = res * (self.__dimensions[2]-1)
            
        return res
        
    def GetDimensionality(self):
        return 3
            
    def GetMultiIndexOfNode(self,index):
        index = int(index)
        planesize = self.__dimensions[1] *self.__dimensions[2]
        nx = index / planesize
        resyz = index - nx*(planesize)
        ny = resyz /self.__dimensions[2]
        nz =  resyz - ny*self.__dimensions[2]
        return np.array([nx,ny,nz],dtype= np.int_)
        
    def GetMonoIndexOfNode(self,indexs):
        planesize = self.__dimensions[1] *self.__dimensions[2]
        return planesize*indexs[0]+indexs[1]*self.__dimensions[2] +indexs[2]
        
     
    def GetMultiIndexOfElement(self,index):
        index = int(index)
        planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)
        nx = index / planesize
        resyz = index - nx*(planesize)
        ny = resyz /(self.__dimensions[2]-1)
        nz =  resyz - ny*(self.__dimensions[2]-1)
        return np.array([nx,ny,nz])        
        
    def GetMonoIndexOfElement(self,indexs):
        planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)
        return indexs[0]*planesize+indexs[1]*(self.__dimensions[2]-1) +indexs[2]
       
    def GetPosOfNode(self,index):
        [nx,ny,nz] = self.GetMultiIndexOfNode(index)
        return np.multiply([nx,ny,nz],self.__spacing)+self.__origin
        
    def GetPosOfNodes(self):
        np.meshgrid()
        x = np.arange(self.__dimensions[0])*self.__spacing[0]+self.__origin[0]
        y = np.arange(self.__dimensions[1])*self.__spacing[1]+self.__origin[1]
        z = np.arange(self.__dimensions[2])*self.__spacing[2]+self.__origin[2]
        #x = np.arange(self.__origin[0],(self.__dimensions[0])*self.__spacing[0]+self.__origin[0],self.__spacing[0])
        #y = np.arange(self.__origin[1],(self.__dimensions[1])*self.__spacing[1]+self.__origin[1],self.__spacing[1])
        #z = np.arange(self.__origin[2],(self.__dimensions[2])*self.__spacing[2]+self.__origin[2],self.__spacing[2])
        #print(x.shape)
        #print(y.shape)
        #print(z.shape)
        xv, yv, zv = np.meshgrid(x, y,z,indexing='ij')
        return np.array([xv.ravel(),yv.ravel(),zv.ravel()]).T
        
    def GetElementAtPos(self,pos):
        pos = pos-self.__origin
        pos /= self.__spacing
        #print("pos" + str(pos))
        #elemindex = np.ceil(pos).astype(int)
        elemindex = pos.astype(int)
        ## border effect correction 
        
        elemindex = np.minimum(elemindex, self.__dimensions-2)
        elemindex = np.maximum(elemindex, [0,0,0])
        
        return self.GetMonoIndexOfElement(elemindex)
    
    def GetElementShapeFunctionsAtPos(self,el, pos):
        
        
        coon = self.GetConnectivityForElement(el)
        p0 = self.GetPosOfNode(coon[0])
        n0 = (pos-p0)*2./self.__spacing - 1.
        myElem = Hexa8Cuboid()
        myElem.delta = self.__spacing
        return myElem.GetShapeFunc(n0)
        
    def GetValueAtPos(self,field,pos):
        
            el = self.GetElementAtPos(pos)
            coon = self.GetConnectivityForElement(el)
            xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)
            return field[coon].dot(xiChiEta)
        
    def GetConnectivityForElement(self, index):
        exyz = self.GetMultiIndexOfElement(index)
        
        n0 = exyz[0]*self.__dimensions[1]*self.__dimensions[2] +exyz[1]*self.__dimensions[2] + exyz[2]
        n1 = n0 + self.__dimensions[1]*self.__dimensions[2]
        n2 = n1 + self.__dimensions[2]
        n3 = n0 + self.__dimensions[2]
        
        n4 = n0 + 1
        n5 = n1 + 1
        n6 = n2 + 1
        n7 = n3 + 1
        
        ##u0,u1,u3
        return np.array([n0, n1, n2, n3, n4, n5, n6, n7])                
    
    def __str__(self):
        res = ''
        res = res + ("Number of Nodes    : "+str(self.GetNumberOfNodes())) + "\n"
        res = res + ("Number of Elements : "+str(self.GetNumberOfElements())) + "\n"
        res = res + ("dimensions         : "+str(self.__dimensions ))         + "\n"
        res = res + ("origin             : "+str(self.__origin)) + "\n"
        res = res + ("spacing            : "+str(self.__spacing)) + "\n"
        return res
        

def CheckIntegrity():
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2,2]);
    myMesh.SetSpacing([1, 1, 1]);
    #myMesh.SetOrigin([-2.5,-1.2,-1.5]);
    
    print(myMesh)
    print(myMesh.GetConnectivityForElement(3))
    
    np.set_printoptions(threshold='nan')
    newmesh = myMesh.GetSubSuperMesh([3,3,3])
    
    oldtonew = myMesh.GetNodeTrasfertMatrix(newmesh)
    print(newmesh)
    print(oldtonew.todense())
    print(newmesh.GetNodeTrasfertMatrix(myMesh).todense())
    
    print(myMesh.GetElementTrasfertMatrix(newmesh).todense())
    print(newmesh.GetElementTrasfertMatrix(myMesh).todense())
    print(newmesh.GetDimensionality())
    print(myMesh.GetValueAtPos(np.array([1,2,3,4,5,6,7,8]),[0.5,0.5,0.5]))
    return "OK"
    
if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover 