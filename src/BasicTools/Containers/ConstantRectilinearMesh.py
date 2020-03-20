# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np


from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tags
from BasicTools.Containers.UnstructuredMesh import ElementsContainer as ElementsContainer
from BasicTools.Containers.UnstructuredMesh import AllElements as AllElements
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Containers.UnstructuredMeshModificationTools import ComputeSkin

class ConstantRectilinearElementContainer(BaseOutputObject):
    def __init__(self,__dimensions):
        super(ConstantRectilinearElementContainer,self).__init__(None)
        #self.caller = caller
        self.__dimensions = None
        self.SetDimensions(__dimensions)
        self.tags = Tags()
        self._connectivity = None
        self.mutable = False
        self.space = None
    @property
    def connectivity(self):
        if(self._connectivity is None):
            self._connectivity = np.empty((self.GetNumberOfElements(),self.GetNumberOfNodesPerElement() ), dtype=np.int)
            for i in range(self.GetNumberOfElements()):
                self._connectivity[i,:] = self.GetConnectivityForElement(i)
        return self._connectivity

    def SetDimensions(self,data):
        if self.__dimensions is None:
            self.__dimensions = np.array(data,int);
        else:
            if len(self.__dimensions) != len(data):
                raise(Exception("Cant change the dimensionality after creation "))
            else:
                self.__dimensions = np.array(data,int);

        self.nodesPerElement = 2**len(self.__dimensions)

        if len(self.__dimensions)  == 3:
            self.elementType = ElementNames.Hexaedron_8
            from BasicTools.FE.Spaces.HexaSpaces import Hexa_P1_Lagrange
            self.space = Hexa_P1_Lagrange()
        elif len(self.__dimensions) == 2 :
            self.elementType = ElementNames.Quadrangle_4
            from BasicTools.FE.Spaces.QuadSpaces import Quad_P1_Lagrange
            self.space = Quad_P1_Lagrange()
        else:
             raise(Exception("cant build a mesh of this dimensionality"))

    def GetDimensionality(self):
        return len(self.__dimensions)

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

    def GetNumberOfElements(self,dim=None):
        if dim is None:
            dim = len(self.__dimensions)

        res = 1;

        if self.__dimensions[0] >= 1:
            res = res * (self.__dimensions[0]-1)

        if self.__dimensions[1] >= 1:
            res = res * (self.__dimensions[1]-1)

        if dim == 2:
            return res

        if self.__dimensions[2] >= 1:
            res = res * (self.__dimensions[2]-1)

        if dim == 3:
            return  res

    def GetNumberOfNodesPerElement(self):
        return 2**len(self.__dimensions)

    def GetTag(self, tagName):
        """
        return the tag based is a name
        if the tag does not exist a new tag is created
        """
        return self.tags.CreateTag(tagName,False)
    def __str__(self):
        res  = "    ConstantRectilinearElementContainer, "
        res += "  Type : ({},{}), ".format(self.elementType,self.GetNumberOfElements())
        res += "  Tags : " + " ".join([ ("("+x.name+":"+str(len(x)) +")") for x in self.tags]) + "\n"
        return res

class ConstantRectilinearMesh(MeshBase):

    def IsConstantRectilinear(self):
        return True

    def __init__(self,dim = 3):
        super(ConstantRectilinearMesh,self).__init__()
        #Number of nodes
        self.__dimensions = np.ones((dim,),dtype=int)*2;
        self.__origin = np.zeros((dim,) )
        self.__spacing = np.ones((dim,))
        self.nodes = None
        self.elements = AllElements()
        self.structElements = ConstantRectilinearElementContainer(self.__dimensions)
        self.elements[self.structElements.elementType] = self.structElements

    def GetNamesOfElemTagsBulk(self):
        return [ tag.name for tag in self.structElements.tags]

    def GetElementsInTagBulk(self,tagname):
        return self.structElements.tags[tagname].GetIds()

    def SetDimensions(self,data):
        self.__dimensions = np.array(data,int);
        self.structElements.SetDimensions(self.__dimensions)
        self.nodes = None

    def GetDimensions(self):
        return np.array(self.__dimensions)

    def SetSpacing(self,data):
        self.__spacing = np.array(data, "float");
        self.nodes = None

    def GetSpacing(self):
        return self.__spacing;

    def GetdV(self):
        """Get the volume of one element."""
        return np.prod(self.GetSpacing())

    def SetOrigin(self,data):
        self.__origin = np.array(data);
        self.nodes = None

    def GetOrigin(self):
        return self.__origin;

    @property
    def boundingMin(self):
        return self.GetOrigin()

    @property
    def boundingMax(self):
        return self.GetOrigin() + (self.GetDimensions()-1)*self.GetSpacing()

    def GetNumberOfNodes(self):
        return np.prod(self.__dimensions)

    def GetNumberOfElements(self,dim=None):
        return self.structElements.GetNumberOfElements()

    def GetMultiIndexOfElement(self,index):
        return self.structElements.GetMultiIndexOfElement(index)

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
              self.nodes = np.empty((self.GetNumberOfNodes(),2))
              self.nodes[:,0] = xv.ravel()
              self.nodes[:,1] = yv.ravel()
              return self.nodes

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

        return self.structElements.space.GetShapeFunc(n0)

    def GetValueAtPos(self,field,pos):
        el = self.GetElementAtPos(pos)
        coon = self.GetConnectivityForElement(el)
        xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)
        return field[coon].dot(xiChiEta)

    def GetValueAtPosMultipleField(self,fields,pos):
        """
        All fields must have same mask (i.e. NaN values outside the support)
        """

        alPos = [pos]

        dim = len(self.__dimensions)
        for i in range(2**dim):
            tempPos = 1.*pos
            res0 = ""
            for j in range(dim-len(bin(i)[2:])):
              res0 += "0"
            res0 += bin(i)[2:]
            for j,b in enumerate(res0):
                if b=="0":
                    tempPos[j] -= self.__spacing[j]
                else:
                    tempPos[j] += self.__spacing[j]
            alPos.append(tempPos)

        acceptableSolution = False

        count = 0
        while acceptableSolution == False:
            pos = alPos[count]
            res = []
            el = self.GetElementAtPos(pos)
            coon = self.GetConnectivityForElement(el)
            xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)

            locFfield = fields[0][coon]
            nans = np.argwhere(np.isnan(locFfield))
            notNans = np.argwhere(~np.isnan(locFfield))

            if notNans.shape[0] != 0:
                acceptableSolution = True
                for i in range(fields.shape[0]):
                    locFfield = fields[i][coon]
                    locFfield[nans] = np.mean(locFfield[notNans])
                    res.append(locFfield.dot(xiChiEta))
            count += 1
        return res

    def GetElementShapeFunctionsDerAtPos(self, el,pos):
        coon = self.GetConnectivityForElement(el)
        p0 = self.GetPosOfNode(coon[0])
        n0 = (pos-p0)*2./self.__spacing - 1.
        f = 2./self.__spacing
        return f[:,np.newaxis]*self.structElements.space.GetShapeFuncDer(n0)

    def GetDefValueAtPos(self,field,pos):
        el = self.GetElementAtPos(pos)
        coon = self.GetConnectivityForElement(el)
        xiChiEta = self.GetElementShapeFunctionsDerAtPos(el,pos)
        return field[coon].dot(xiChiEta.T)

    def GetConnectivityForElement(self, index):
        return self.structElements.GetConnectivityForElement(index)

    def GenerateFullConnectivity(self):
        return self.structElements.connectivity

    def ComputeGlobalOffset(self):
        """
        Recompute the Global Offset,
        This is necessary for some operation.
        Recomendation : Call it after changing the topology
        """

        cpt = 0
        for type, data in self.elements.items():
            data.globaloffset = cpt
            n = data.GetNumberOfElements()
            cpt = cpt + n

    def __str__(self):
        res = ''
        res  = "ConstantRectilinearMesh \n"
        res = res + "  Number of Nodes    : "+str(self.GetNumberOfNodes()) + "\n"
        res += "    Tags : " + " ".join( ["("+x.name+":"+str(len(x))+")" for x in  self.nodesTags ]) + "\n"

        res = res + "  Number of Elements : "+str(self.GetNumberOfElements()) + "\n"
        res = res + "  dimensions         : "+str(self.__dimensions )         + "\n"
        res = res + "  origin             : "+str(self.__origin) + "\n"
        res = res + "  spacing            : "+str(self.__spacing) + "\n"
        for name,data in self.elements.items():
            res += str(data)
        res += "\n"
        res += "  Node Tags          : " + str(self.nodesTags) + "\n"
        res += "  Cell Tags          : " + str([x for x in self.GetNamesOfElemTags()])+ "\n"
        if len(self.nodeFields.keys()):
            res += "  nodeFields         : " + str(list(self.nodeFields.keys())) + "\n"
        if len(self.elemFields.keys()):
            res += "  elemFields         : " + str(list(self.elemFields.keys())) + "\n"
        return res


def CheckIntegrity():
    import sys

    # Error checking tests
    try:
        # not implemented for dim = 1 this line must fail
        myMesh = ConstantRectilinearMesh( dim=1)
        raise("Error detecting bad argument") # pragma: no cover
    except:
        pass

    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([1,1,1]);
    myMesh.SetSpacing([1, 1, 1]);


    try:
        # not implemented for dim = 1 this line must fail
        myMesh.GetNodalIndicesOfBorder()
        raise("Error detecting bad mesh props")# pragma: no cover
    except:
        pass


    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2,2]);
    myMesh.SetSpacing([1, 1, 1]);
    #myMesh.SetOrigin([-2.5,-1.2,-1.5]);

    print(myMesh)
    print(myMesh.elements[ElementNames.Hexaedron_8].GetNumberOfElements())
    print(myMesh.elements[ElementNames.Hexaedron_8].GetNumberOfNodesPerElement())
    print((myMesh.IsConstantRectilinear()))
    print((myMesh.GetNamesOfElemTags()))
    print((myMesh.GetDimensions()))
    print((myMesh.GetMonoIndexOfNode(np.array([0,0,0]) )))
    print((myMesh.GetMonoIndexOfNode(np.array([[0,0,0],[1,1,1]]) )))
    print((myMesh.GetPosOfNodes()))
    print(myMesh.GetClosestPointToPos([0,0.5,1.5]))
    print(myMesh.GetClosestPointToPos([0,0.5,1.5],MultiIndex=True))

    print(myMesh.GetElementAtPos([0,0.5,1.5]))
    print(myMesh.GetElementAtPos([0,0.5,1.5],MultiIndex=True))

    myMesh.elements[ElementNames.Hexaedron_8].tags.CreateTag("TestTag",False).SetIds([0,1])

    if len(myMesh.GetElementsInTagBulk("TestTag")) != 2 :
        raise(Exception("Tag system not working corretly") )# pragma: no cover



    print((myMesh.GetConnectivityForElement(0)))

    print((myMesh.GetElementShapeFunctionsDerAtPos(0,[0.5,0.5,0.5] )))

    print(myMesh.GenerateFullConnectivity())

    np.set_printoptions(threshold=sys.maxsize)

    print((myMesh.GetValueAtPos(np.array([1,2,3,4,5,6,7,8]),[0.5,0.5,0.5])))


    res = (myMesh.GetNodalIndicesOfBorder(0))
    print(res)
    print(res.size)


    print("-----------------2D const rectilinear mesh------------------------")
    myMesh = ConstantRectilinearMesh(dim=2)
    myMesh.SetDimensions([3,3]);
    myMesh.SetSpacing([1, 1]);
    myMesh.SetOrigin([0.,0.]);

    print(myMesh)

    print((myMesh.GetMonoIndexOfNode(np.array([0,0]) )))
    print((myMesh.GetMonoIndexOfNode(np.array([[0,0],[1,1]]) )))
    print((myMesh.GetPosOfNodes()))

    print((myMesh.GetConnectivityForElement(0)))
    print((myMesh.GetElementShapeFunctionsDerAtPos(0,[0.5,0.5] )))



    np.set_printoptions(threshold=sys.maxsize)


    print((myMesh.GetValueAtPos(np.array([1,2,3,4,5,6,7,8,9]),[0.5,0.5])))
    print((myMesh.GetDefValueAtPos(np.array([1,2,3,4,5,6,7,8,9]),[0.5,0.5] )))


    res = (myMesh.GetNodalIndicesOfBorder(0))
    print(res)
    if np.any(np.sort(res) != [0, 1, 2, 3, 5, 6, 7, 8 ]): # pragma: no cover
        return "Not Ok on 'GetNodalIndicesOfBorder(0)'"

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
