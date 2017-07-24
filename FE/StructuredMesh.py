# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

from BasicTools.FE.MeshBase import MeshBase

class StructuredMesh(MeshBase):

    def IsStructured(self):
        return True

    def __init__(self,dim = 3):
        super(StructuredMesh,self).__init__()

        if dim == 3:
            self.nodes = np.array( [[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1.1,0,0],
                              [1.2,0,1],
                              [1.3,1,0],
                              [1.4,1,1]])
        else :
            self.nodes = np.array( [[0,0,0],
                              [0,1,0],
                              [1.1,0,0],
                              [1.2,1,0]])

        self.__dimensions = np.ones((dim,),dtype=int)*2;
        self.elemTags = {}

    def GetNamesOfElemTags(self):
        return self.elemTags.keys()

    def SetDimensions(self,data):
        self.__dimensions = np.array(data,"int");

    def GetDimensions(self):
        return self.__dimensions ;

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
            nx = index // planesize
            resyz = index - nx*(planesize)
            ny = resyz //self.__dimensions[2]
            nz =  resyz - ny*self.__dimensions[2]
            return np.array([nx,ny,nz],dtype= np.int_)
        else:
            nx = index // self.__dimensions[1]
            ny = index - nx*(self.__dimensions[1])
            return np.array([nx,ny],dtype= np.int_)

    def GetMonoIndexOfNode(self,indexs):
        if self.GetDimensionality() == 3:
            planesize = self.__dimensions[1] *self.__dimensions[2]
            return planesize*indexs[0]+indexs[1]*self.__dimensions[2] +indexs[2]
        else:
            planesize = self.__dimensions[1]
            return planesize*indexs[0]+indexs[1]

    def GetMultiIndexOfElement(self,index):
        index = int(index)
        if self.GetDimensionality() == 3:
            planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)
            nx = index // planesize
            resyz = index - nx*(planesize)
            ny = resyz //(self.__dimensions[2]-1)
            nz =  resyz - ny*(self.__dimensions[2]-1)
            return np.array([nx,ny,nz])
        else:
            planesize = (self.__dimensions[1]-1)
            nx = index // planesize
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
        return self.nodes[index,:]

    def GetPosOfNodes(self):
        return self.nodes

    def GetConnectivityForElement(self, index):
        exyz = self.GetMultiIndexOfElement(index)

        if self.GetDimensionality() == 3:
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
        else:
            n0 = exyz[0]*self.__dimensions[1] +exyz[1]
            n1 = n0 + self.__dimensions[1]
            n2 = n1 + 1
            n3 = n0 + 1
            ##u0,u1,u3
            return np.array([n0, n1, n2, n3])

    def __str__(self):
        res = ''
        res  = "StructuredMesh \n"
        res = res + ("  Number of Nodes    : "+str(self.GetNumberOfNodes())) + "\n"
        res = res + ("  Number of Elements : "+str(self.GetNumberOfElements())) + "\n"
        res = res + ("  dimensions         : "+str(self.__dimensions ))         + "\n"
        return res


def CheckIntegrity():
    print("-----------------3D const rectilinear mesh------------------------")
    myMesh = StructuredMesh()
    myMesh.SetDimensions([2,2,2]);
    myMesh.nodes = np.array( [[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1.1,0,0],
                              [1.2,0,1],
                              [1.3,1,0],
                              [1.4,1,1]])

    print(myMesh)
    print(myMesh.GetConnectivityForElement(0))
    print(myMesh.GetMultiIndexOfNode(0))
    print(myMesh.GetMonoIndexOfNode([1,0,0]))
    print(myMesh.GetMonoIndexOfElement([0,0,0]))
    print(myMesh.GetPosOfNode(2))
    print("-----------------2D const rectilinear mesh------------------------")
    myMesh = StructuredMesh(2)
    myMesh.SetDimensions([2,2]);
    myMesh.nodes = np.array( [[0,0,0],
                              [0,1,0],
                              [1.1,0,0],
                              [1.2,1,0]])

    print(myMesh)
    print(myMesh.GetConnectivityForElement(0))
    print(myMesh.GetMultiIndexOfNode(0))
    print(myMesh.GetMonoIndexOfNode([1,0]))
    print(myMesh.GetMonoIndexOfElement([0,0]))
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
