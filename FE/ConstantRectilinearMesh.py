# -*- coding: utf-8 -*-

from BasicTools.FE.MeshBase import MeshBase
from BasicTools.FE.Hexa8Cuboid import Hexa8Cuboid
from BasicTools.FE.Quad4Rectangle import Quad4Rectangle

import numpy as np
from scipy.sparse import coo_matrix

class ConstantRectilinearMesh(MeshBase):

    def IsConstantRectilinear(self):
        return True

    def __init__(self,dim = 3):
        super(ConstantRectilinearMesh,self).__init__()
        #Number of nodes
        self.__dimensions = np.ones((dim,),dtype=int)*2;
        self.__origin = np.zeros((dim,) )
        self.__spacing = np.ones((dim,))
        self.elemTags = {}
        self.connectivity = None
        self.nodes = None

    def GetNamesOfElemTags(self):
        return self.elemTags.keys()

    def GetSubSuperMesh(self,_newDimensions):
        newDimensions = np.array(_newDimensions)
        ## to generate meshes with more or less elements in each directions
        ## return the mesh
        newSpac = (self.__dimensions-1)*self.__spacing/(newDimensions-1)

        res = type(self)(dim=self.GetDimensionality())
        res.SetSpacing(newSpac)
        res.SetDimensions(newDimensions)
        res.SetOrigin(self.GetOrigin())

        return res

    def GetNodeTrasfertMatrix(self, destination):
        # newVector   = oldToNew * oldVector
        nbNodes = 2**self.GetDimensionality()
        oldToNewVals = np.zeros((destination.GetNumberOfNodes(),nbNodes))
        oldToNewIK = np.zeros((destination.GetNumberOfNodes(),nbNodes), dtype=np.int_)
        oldToNewJK = np.zeros((destination.GetNumberOfNodes(),nbNodes), dtype=np.int_)

        for i in  xrange(destination.GetNumberOfNodes()):

            pos= destination.GetPosOfNode(i)
            el = self.GetElementAtPos(pos)
            coon = self.GetConnectivityForElement(el)
            xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)
            oldToNewVals[i,:] = xiChiEta
            oldToNewIK[i,:] = i
            oldToNewJK[i,:] = coon

        oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destination.GetNumberOfNodes(), self.GetNumberOfNodes())).tocsc()

        return oldToNew

    def GetElementTrasfertMatrix(self, destination):

        nps = 10
        nps3 = nps**self.GetDimensionality()
        oldToNewVals = np.zeros((destination.GetNumberOfNodes(),nps3))
        oldToNewIK = np.zeros((destination.GetNumberOfNodes(),nps3), dtype=np.int_)
        oldToNewJK = np.zeros((destination.GetNumberOfNodes(),nps3), dtype=np.int_)

        for i in  xrange(destination.GetNumberOfElements()):
            coon = destination.GetConnectivityForElement(i)

            n0pos = destination.GetPosOfNode(coon[0])
            cpt =0
            for cx in range(0,nps):
                for cy in range(0,nps):
                    if self.GetDimensionality() == 3:
                        for cz in range(0,nps):
                            pos = n0pos + destination.GetSpacing()*([cx+0.5,cy+0.5,cz+0.5])/nps
                            el = self.GetElementAtPos(pos)
                            oldToNewVals[i,cpt] += 1./nps3
                            oldToNewIK[i,cpt] += i
                            oldToNewJK[i,cpt] += el
                            cpt +=1
                    else:
                        pos = n0pos + destination.GetSpacing()*([cx+0.5,cy+0.5])/nps
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
        oldToNew =  coo_matrix((oldToNewVals.ravel(), (oldToNewIK.flatten(), oldToNewJK.flatten())), shape=(destination.GetNumberOfElements(), self.GetNumberOfElements())).tocsc()

        return oldToNew

    def SetDimensions(self,data):
        self.__dimensions = np.array(data,"int");
        self.nodes = None

    def GetDimensions(self):
        return self.__dimensions ;

    def SetSpacing(self,data):
        self.__spacing = np.array(data, "float");
        self.nodes = None

    def GetSpacing(self):
        return self.__spacing;

    def SetOrigin(self,data):
        self.__origin = np.array(data);
        self.nodes = None

    def GetOrigin(self):
        return self.__origin;


    def GetNumberOfNodes(self):
        return np.prod(self.__dimensions)

    def GetNumberOfElements(self):
        res = 1;
        if self.__dimensions[0] >= 1:
            res = res * (self.__dimensions[0]-1)

        if self.__dimensions[1] >= 1:
            res = res * (self.__dimensions[1]-1)

        if self.GetDimensionality() == 2:
            return res

        if self.__dimensions[2] >= 1:
            res = res * (self.__dimensions[2]-1)

        return res

    def GetDimensionality(self):
        return len(self.__dimensions)

    def GetMultiIndexOfNode(self,index):
        index = int(index)
        if self.GetDimensionality() == 3:
            planesize = self.__dimensions[1] *self.__dimensions[2]
            nx = index / planesize
            resyz = index - nx*(planesize)
            ny = resyz /self.__dimensions[2]
            nz =  resyz - ny*self.__dimensions[2]
            return np.array([nx,ny,nz],dtype= np.int_)
        else:
            nx = index / self.__dimensions[1]
            ny = index - nx*(self.__dimensions[1])
            return np.array([nx,ny],dtype= np.int_)

    def GetMonoIndexOfNode(self,_indexs):
        _indexs= np.asarray(_indexs)
        if len(_indexs.shape) == 1:
            indexs = _indexs[np.newaxis]
        else:
            indexs = _indexs

        if self.GetDimensionality() == 3:

            planesize = self.__dimensions[1] *self.__dimensions[2]
            res = planesize*indexs[:,0]
            res += indexs[:,1]*self.__dimensions[2]
            res += indexs[:,2]
            return res
        else:
            planesize = self.__dimensions[1]
            return planesize*indexs[:,0]+indexs[:,1]

    def GetMultiIndexOfElement(self,index):
        index = int(index)
        if self.GetDimensionality() == 3:
            planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)

            res = np.empty(3,dtype=int)
            res[0] = index / planesize
            resyz = index - res[0]*(planesize)
            res[1] = resyz /(self.__dimensions[2]-1)
            res[2] =  resyz - res[1]*(self.__dimensions[2]-1)
            return res

        else:
            planesize = (self.__dimensions[1]-1)
            nx = index / planesize
            ny = index - nx*(planesize)
            return np.array([nx,ny])

    def GetMonoIndexOfElement(self,indexs):
        if self.GetDimensionality() == 3:
            planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)
            return indexs[0]*planesize+indexs[1]*(self.__dimensions[2]-1) +indexs[2]
        else :
            planesize = (self.__dimensions[1]-1)
            return indexs[0]*planesize+indexs[1]

    def GetPosOfNode(self,index):
        if self.GetDimensionality() == 3 :
            [nx,ny,nz] = self.GetMultiIndexOfNode(index)
            return np.multiply([nx,ny,nz],self.__spacing)+self.__origin
        else:
            [nx,ny] = self.GetMultiIndexOfNode(index)
            return np.multiply([nx,ny],self.__spacing)+self.__origin

    def GetPosOfNodes(self):

        x = np.arange(self.__dimensions[0])*self.__spacing[0]+self.__origin[0]
        y = np.arange(self.__dimensions[1])*self.__spacing[1]+self.__origin[1]
        if self.GetDimensionality() == 2:
          xv, yv = np.meshgrid(x, y,indexing='ij')
          return np.array([xv.ravel(),yv.ravel()]).T

        z = np.arange(self.__dimensions[2])*self.__spacing[2]+self.__origin[2]
        xv, yv, zv = np.meshgrid(x, y,z,indexing='ij')

        if self.nodes is None:
            self.nodes = np.array([xv.ravel(),yv.ravel(),zv.ravel()]).T

        return self.nodes

    def GetElementAtPos(self,pos):
        pos = pos-self.__origin
        pos /= self.__spacing
        elemindex = pos.astype(int)
        elemindex = np.minimum(elemindex, self.__dimensions-2)
        elemindex = np.maximum(elemindex, 0.)
        return self.GetMonoIndexOfElement(elemindex)

    def GetElementShapeFunctionsAtPos(self,el, pos):
        coon = self.GetConnectivityForElement(el)
        p0 = self.GetPosOfNode(coon[0])
        n0 = (pos-p0)*2./self.__spacing - 1.

        if self.GetDimensionality() == 3:
            myElem = Hexa8Cuboid()
        else:
            myElem = Quad4Rectangle()

        myElem.delta = self.__spacing
        return myElem.GetShapeFunc(n0)

    def GetValueAtPos(self,field,pos):
        el = self.GetElementAtPos(pos)
        coon = self.GetConnectivityForElement(el)
        xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)
        return field[coon].dot(xiChiEta)

    def GetElementShapeFunctionsDerAtPos(self, el,pos):
        coon = self.GetConnectivityForElement(el)
        p0 = self.GetPosOfNode(coon[0])
        n0 = (pos-p0)*2./self.__spacing - 1.
        if self.GetDimensionality() ==3:
            myElem = Hexa8Cuboid()
        else:
            myElem = Quad4Rectangle()
        myElem.delta = self.__spacing
        return myElem.ShapeFuncDer(n0)

    def GetDefValueAtPos(self,field,pos):
        el = self.GetElementAtPos(pos)
        coon = self.GetConnectivityForElement(el)
        xiChiEta = self.GetElementShapeFunctionsDerAtPos(el,pos)
        return field[coon].dot(xiChiEta.T)

    def GetConnectivityForElement(self, index):
        exyz = self.GetMultiIndexOfElement(index)

        if self.GetDimensionality() == 3:
            res = np.empty(8,dtype=int)
            #n0
            res[0] = exyz[0]*self.__dimensions[1]*self.__dimensions[2] +exyz[1]*self.__dimensions[2] + exyz[2]
            #n1
            res[1]= res[0] + self.__dimensions[1]*self.__dimensions[2]
            res[2] = res[1] + self.__dimensions[2]
            res[3] = res[0] + self.__dimensions[2]

            res[4:8] = res[0:4] + 1
            #n5 = n1 + 1
            #n6 = n2 + 1
            #n7 = n3 + 1

            ##u0,u1,u3
            return res
            #np.array([n0, n1, n2, n3, n4, n5, n6, n7])
        else:
            n0 = exyz[0]*self.__dimensions[1] +exyz[1]
            n1 = n0 + self.__dimensions[1]
            n2 = n1 + 1
            n3 = n0 + 1
            ##u0,u1,u3
            return np.array([n0, n1, n2, n3])

    def GenerateFullConnectivity(self):
        if(self.connectivity is None):
            self.connectivity = np.empty((self.GetNumberOfElements(),2**self.GetDimensionality() ), dtype=np.int)
            for i in xrange(self.GetNumberOfElements()):
                self.connectivity[i,:] = self.GetConnectivityForElement(i)
        from BasicTools.FE.UnstructuredMesh import ElementsContainer as ElementsContainer
        import BasicTools.FE.ElementNames
        self.elements = {}
        self.elements[BasicTools.FE.ElementNames.Hexaedron_8 ] = ElementsContainer(BasicTools.FE.ElementNames.Hexaedron_8)
        self.elements[BasicTools.FE.ElementNames.Hexaedron_8 ].connectivity = self.connectivity
        self.elements[BasicTools.FE.ElementNames.Hexaedron_8 ].tags =  self.elemTags
        return self.connectivity


    def __str__(self):
        res = ''
        res  = "ConstantRectilinearMesh \n"
        res = res + ("  Number of Nodes    : "+str(self.GetNumberOfNodes())) + "\n"
        res = res + ("  Number of Elements : "+str(self.GetNumberOfElements())) + "\n"
        res = res + ("  dimensions         : "+str(self.__dimensions ))         + "\n"
        res = res + ("  origin             : "+str(self.__origin)) + "\n"
        res = res + ("  spacing            : "+str(self.__spacing)) + "\n"
        return res


def CheckIntegrity():
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2,2]);
    myMesh.SetSpacing([1, 1, 1]);
    #myMesh.SetOrigin([-2.5,-1.2,-1.5]);

    print(myMesh)
    print(myMesh.IsConstantRectilinear())
    print(myMesh.GetNamesOfElemTags())
    print(myMesh.GetDimensions())
    print(myMesh.GetMonoIndexOfNode(np.array([0,0,0]) ))
    print(myMesh.GetMonoIndexOfNode(np.array([[0,0,0],[1,1,1]]) ))
    print(myMesh.GetPosOfNodes())

    print(myMesh.GetConnectivityForElement(0))

    print(myMesh.GetElementShapeFunctionsDerAtPos(0,[0.5,0.5,0.5] ))

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

    print("-----------------2D const rectilinear mesh------------------------")
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2]);
    myMesh.SetSpacing([1, 1]);
    myMesh.SetOrigin([0.,0.]);

    print(myMesh)

    print(myMesh.GetMonoIndexOfNode(np.array([0,0]) ))
    print(myMesh.GetMonoIndexOfNode(np.array([[0,0],[1,1]]) ))
    print(myMesh.GetPosOfNodes())

    print(myMesh.GetConnectivityForElement(0))
    print(myMesh.GetElementShapeFunctionsDerAtPos(0,[0.5,0.5] ))



    np.set_printoptions(threshold='nan')
    newmesh = myMesh.GetSubSuperMesh([3,3])

    oldtonew = myMesh.GetNodeTrasfertMatrix(newmesh)
    print(newmesh)
    print(oldtonew.todense())
    print(newmesh.GetNodeTrasfertMatrix(myMesh).todense())

    print(myMesh.GetElementTrasfertMatrix(newmesh).todense())
    print(newmesh.GetElementTrasfertMatrix(myMesh).todense())
    print(newmesh.GetDimensionality())
    print(myMesh.GetValueAtPos(np.array([1,2,3,4]),[0.5,0.5]))
    print(myMesh.GetDefValueAtPos(np.array([1,2,3,4]),[0.5,0.5] ))
    return "OK"

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover
