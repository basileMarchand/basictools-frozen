# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.FE.Hexa8Cuboid import Hexa8Cuboid
from BasicTools.FE.Quad4Rectangle import Quad4Rectangle
from BasicTools.Containers.MeshBase import Tags
from BasicTools.Containers.UnstructuredMesh import ElementsContainer as ElementsContainer
from BasicTools.Containers.UnstructuredMesh import AllElements as AllElements
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject


class ConstantRectilinearElementContainer(BaseOutputObject):
    def __init__(self,caller):
        super(ConstantRectilinearElementContainer,self).__init__(None)
        self.caller = caller
        self.tags = self.caller.elemTags
        self._connectivity = None

        if caller.GetDimensionality() == 3:
            self.elementType = ElementNames.Hexaedron_8
        elif caller.GetDimensionality() == 2 :
            self.elementType = ElementNames.Quadrangle_4
        else:
             raise(Exception("cant build a mesh of this dimensionality"))

    @property
    def connectivity(self):
        if(self._connectivity is None):
            self._connectivity = np.empty((self.caller.GetNumberOfElements(),2**self.caller.GetDimensionality() ), dtype=np.int)
            for i in range(self.caller.GetNumberOfElements()):
                self._connectivity[i,:] = self.caller.GetConnectivityForElement(i)
        return self._connectivity

    def GetNumberOfElements(self):
        return self.caller.GetNumberOfElements()

    def GetNumberOfNodesPerElement(self):
        return 2*self.caller.GetDimensionality()

class ConstantRectilinearMesh(MeshBase):

    def IsConstantRectilinear(self):
        return True

    def __init__(self,dim = 3):
        super(ConstantRectilinearMesh,self).__init__()
        #Number of nodes
        self.__dimensions = np.ones((dim,),dtype=int)*2;
        self.__origin = np.zeros((dim,) )
        self.__spacing = np.ones((dim,))
        self.elemTags = Tags()
        self.nodes = None
        self.elements = AllElements()
        structElements = ConstantRectilinearElementContainer(self)
        self.elements[structElements.elementType] = structElements

    def GetNamesOfElemTags(self):
        return self.elemTags.keys()

    def GetElementsInTag(self,tagname):
        return self.elemTags[tagname].GetIds()

    def SetDimensions(self,data):
        self.__dimensions = np.array(data,int);
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

    def GetNumberOfElements(self,dim=None):
        if dim is None:
            dim = self.GetDimensionality()
        res = 1;
        if self.__dimensions[0] >= 1:
            res = res * (self.__dimensions[0]-1)

        if self.__dimensions[1] >= 1:
            res = res * (self.__dimensions[1]-1)

        if self.GetDimensionality() == 2 :
            if dim == 2:
                return res
            else:
                return 0

        if self.__dimensions[2] >= 1:
            if dim == 3:
                return  res * (self.__dimensions[2]-1)
            else:
                return 0

    def GetDimensionality(self):
        return len(self.__dimensions)

    def GetMultiIndexOfNode(self,index):
        index = int(index)
        if self.GetDimensionality() == 3:
            planesize = self.__dimensions[1] *self.__dimensions[2]
            nx = index // planesize
            resyz = index - nx*(planesize)
            ny = resyz // self.__dimensions[2]
            nz =  resyz - ny*self.__dimensions[2]
            return np.array([nx,ny,nz],dtype= np.int_)
        else:
            nx = index // self.__dimensions[1]
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
            res[0] = index // planesize
            resyz = index - res[0]*(planesize)
            res[1] = resyz //(self.__dimensions[2]-1)
            res[2] =  resyz - res[1]*(self.__dimensions[2]-1)
            return res

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
        if self.nodes is not None:
            return self.nodes[index,:]

        if self.GetDimensionality() == 3 :
            [nx,ny,nz] = self.GetMultiIndexOfNode(index)
            return np.multiply([nx,ny,nz],self.__spacing)+self.__origin
        else:
            [nx,ny] = self.GetMultiIndexOfNode(index)
            return np.multiply([nx,ny],self.__spacing)+self.__origin

    def GetPosOfNodes(self):
        """
        Space coordinates for all nodes in the mesh.

        Returns
        -------
        numpy.array
            A 2-dimensional array, the first axis corresponds to the node
            index, the second axis corresponds to space dimension index.
        """
        if self.nodes is None:
            x = np.arange(self.__dimensions[0])*self.__spacing[0]+self.__origin[0]
            y = np.arange(self.__dimensions[1])*self.__spacing[1]+self.__origin[1]
            if self.GetDimensionality() == 2:
              xv, yv = np.meshgrid(x, y,indexing='ij')
              return np.array([xv.ravel(),yv.ravel()]).T

            z = np.arange(self.__dimensions[2])*self.__spacing[2]+self.__origin[2]
            xv, yv, zv = np.meshgrid(x, y,z,indexing='ij')

            self.nodes = np.empty((self.GetNumberOfNodes(),3))
            self.nodes[:,0] = xv.ravel()
            self.nodes[:,1] = yv.ravel()
            self.nodes[:,2] = zv.ravel()

        return self.nodes

    def GetNodalIndicesOfBorder(self,border=0):

        dim =  np.maximum(self.__dimensions-border*2,0)

        if np.any(dim <= 1):
            raise(Exception("Cube to small "))

        def GetMonoIndexOfIndexTensorProduct2D(a,b):
            x,y = np.meshgrid(a,b,indexing='ij')
            faceindexs = (np.hstack((x.flatten()[:,np.newaxis],y.flatten()[:,np.newaxis])))

            face = self.GetMonoIndexOfNode(faceindexs)
            return face

        def GetMonoIndexOfIndexTensorProduct3D(a,b,c):
            x,y,z = np.meshgrid(a,b,c,indexing='ij')
            faceindexs = (np.hstack((x.flatten()[:,np.newaxis],y.flatten()[:,np.newaxis],z.flatten()[:,np.newaxis])))
            face = self.GetMonoIndexOfNode(faceindexs)
            return face

        d2 =  np.maximum(dim-2,0)
        # first and last
        f = border
        l = np.maximum(self.__dimensions-border,f)
        cpt = 0

        if self.GetDimensionality() == 3:
            #the faces, the edges, the corners
            res = np.empty(dim[0]*dim[1]*2+
                           dim[1]*d2[2]*2+
                           d2[0]*d2[2]*2,dtype=np.int)


            face = GetMonoIndexOfIndexTensorProduct3D(range(f,l[0]),range(f,l[1]),[f, l[2]-1])
            res[cpt:cpt+face.size] = face

            cpt += face.size
            face = GetMonoIndexOfIndexTensorProduct3D([f, l[0]-1],range(f,l[1]),range(f+1,l[2]-1))
            res[cpt:cpt+face.size] = face
            cpt += face.size
            face = GetMonoIndexOfIndexTensorProduct3D(range(f+1,l[0]-1),[f,l[1]-1],range(f+1,l[2]-1))
            res[cpt:cpt+face.size] = face
            cpt += face.size
        else:
            #the faces, the edges, the corners
            res = np.empty(dim[0]*2+
                           d2[1]*2,dtype=np.int)


            face = GetMonoIndexOfIndexTensorProduct2D(range(f,l[0]),[f, l[1]-1])
            res[cpt:cpt+face.size] = face

            cpt += face.size
            face = GetMonoIndexOfIndexTensorProduct2D([f, l[0]-1],range(f+1,l[1]-1))
            res[cpt:cpt+face.size] = face
            cpt += face.size
        return res

    def GetClosestPointToPos(self,pos, MultiIndex=False):
        pos = (pos-self.__origin)-self.__spacing/2.
        pos /= self.__spacing
        elemindex = pos.astype(int)
        elemindex = np.minimum(elemindex, self.__dimensions-1)
        elemindex = np.maximum(elemindex, 0)
        if MultiIndex :
            return elemindex

        return self.GetMonoIndexOfNode(elemindex)


    def GetElementAtPos(self,pos, MultiIndex=False):
        pos = pos-self.__origin
        pos /= self.__spacing
        elemindex = pos.astype(int)
        elemindex = np.minimum(elemindex, self.__dimensions-2)
        elemindex = np.maximum(elemindex, 0)
        if MultiIndex :
            return elemindex

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
        import BasicTools.Containers.ElementNames as ElementNames

        if self.GetDimensionality() == 3:
            return self.elements[ElementNames.Hexaedron_8 ].connectivity
        elif self.GetDimensionality() == 2:
            return self.elements[ElementNames.Quadrangle_4].connectivity
        raise(Exception("Invalid dimensionality : " + str(self.GetDimensionality())) )

    def __str__(self):
        res = ''
        res  = "ConstantRectilinearMesh \n"
        res = res + "  Number of Nodes    : "+str(self.GetNumberOfNodes()) + "\n"
        res = res + "  Number of Elements : "+str(self.GetNumberOfElements()) + "\n"
        res = res + "  dimensions         : "+str(self.__dimensions )         + "\n"
        res = res + "  origin             : "+str(self.__origin) + "\n"
        res = res + "  spacing            : "+str(self.__spacing) + "\n"
        return res


def CheckIntegrity():
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2,2]);
    myMesh.SetSpacing([1, 1, 1]);
    #myMesh.SetOrigin([-2.5,-1.2,-1.5]);

    print(myMesh)
    print((myMesh.IsConstantRectilinear()))
    print((myMesh.GetNamesOfElemTags()))
    print((myMesh.GetDimensions()))
    print((myMesh.GetMonoIndexOfNode(np.array([0,0,0]) )))
    print((myMesh.GetMonoIndexOfNode(np.array([[0,0,0],[1,1,1]]) )))
    print((myMesh.GetPosOfNodes()))

    print((myMesh.GetConnectivityForElement(0)))

    print((myMesh.GetElementShapeFunctionsDerAtPos(0,[0.5,0.5,0.5] )))

    print(myMesh.GenerateFullConnectivity())

    np.set_printoptions(threshold=np.nan)

    print((myMesh.GetValueAtPos(np.array([1,2,3,4,5,6,7,8]),[0.5,0.5,0.5])))


    res = (myMesh.GetNodalIndicesOfBorder(0))
    print(res)
    print(res.size)


    print("-----------------2D const rectilinear mesh------------------------")
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([3,3]);
    myMesh.SetSpacing([1, 1]);
    myMesh.SetOrigin([0.,0.]);

    print(myMesh)

    print((myMesh.GetMonoIndexOfNode(np.array([0,0]) )))
    print((myMesh.GetMonoIndexOfNode(np.array([[0,0],[1,1]]) )))
    print((myMesh.GetPosOfNodes()))

    print((myMesh.GetConnectivityForElement(0)))
    print((myMesh.GetElementShapeFunctionsDerAtPos(0,[0.5,0.5] )))



    np.set_printoptions(threshold=np.nan)


    print((myMesh.GetValueAtPos(np.array([1,2,3,4,5,6,7,8,9]),[0.5,0.5])))
    print((myMesh.GetDefValueAtPos(np.array([1,2,3,4,5,6,7,8,9]),[0.5,0.5] )))


    res = (myMesh.GetNodalIndicesOfBorder(0))
    print(res)
    if np.any(np.sort(res) != [0, 1, 2, 3, 5, 6, 7, 8 ]):
        return "Not Ok on 'GetNodalIndicesOfBorder(0)'"

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover