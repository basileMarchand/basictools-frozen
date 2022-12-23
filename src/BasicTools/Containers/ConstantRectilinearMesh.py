# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType
from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tags
from BasicTools.Containers.UnstructuredMesh import AllElements as AllElements
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject, froze_it
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1

@froze_it
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
        self.originalIds = np.empty((0,),dtype=PBasicIndexType)
        self.originalOffset = 0

    @property
    def connectivity(self):
        if(self._connectivity is None):
            self._connectivity = self.GetConnectivityForElements(np.arange(self.GetNumberOfElements()))
            self._connectivity.flags.writeable = False
        return self._connectivity

    def SetDimensions(self,data):
        if self.__dimensions is None:
            self.__dimensions = np.array(data,dtype=PBasicIndexType)
        else:
            if len(self.__dimensions) != len(data):
                raise(Exception("Cant change the dimensionality after creation "))
            else:
                self.__dimensions = np.array(data,dtype=PBasicIndexType)

        self.nodesPerElement = 2**len(self.__dimensions)

        if len(self.__dimensions)  == 3:
            self.elementType = ElementNames.Hexaedron_8
            self.space = LagrangeSpaceP1[ElementNames.Hexaedron_8]
        elif len(self.__dimensions) == 2 :
            self.elementType = ElementNames.Quadrangle_4
            self.space =  LagrangeSpaceP1[ElementNames.Quadrangle_4]
        else:
             raise(Exception("cant build a mesh of this dimensionality"))
        self.space.Create()
        self.originalIds = np.arange(self.GetNumberOfElements(),dtype=PBasicIndexType)


    def GetDimensionality(self):
        return len(self.__dimensions)

    def GetConnectivityForElements(self, indices):

        exyz = self.GetMultiIndexOfElements(np.asarray(indices))

        if self.GetDimensionality() == 3:
            res = np.empty((exyz.shape[0],8),dtype=PBasicIndexType)
            #n0
            res[:,0] = exyz[:,0]*self.__dimensions[1]*self.__dimensions[2] +exyz[:,1]*self.__dimensions[2] + exyz[:,2]
            #n1
            res[:,1]= res[:,0] + self.__dimensions[1]*self.__dimensions[2]
            res[:,2] = res[:,1] + self.__dimensions[2]
            res[:,3] = res[:,0] + self.__dimensions[2]

            res[:,4:8] = res[:,0:4] + 1
            return res
        else:
            res = np.empty((exyz.shape[0],4),dtype=PBasicIndexType)
            res[:,0] = exyz[:,0]*self.__dimensions[1] +exyz[:,1]
            res[:,1] = res[:,0] + self.__dimensions[1]
            res[:,2] = res[:,1] + 1
            res[:,3] = res[:,0] + 1
            return res

    def GetConnectivityForElement(self, index):
        return self.GetConnectivityForElements([index])[0,:]

    def GetMultiIndexOfElements(self,indices):
        indices = np.asarray(indices,dtype=PBasicIndexType)
        if self.GetDimensionality() == 3:
            planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)

            res = np.empty((len(indices),3),dtype=PBasicIndexType)
            res[:,0] = indices // planesize
            resyz = indices - res[:,0]*(planesize)
            res[:,1] = resyz //(self.__dimensions[2]-1)
            res[:,2] =  resyz - res[:,1]*(self.__dimensions[2]-1)
            return res

        else:
            res = np.empty((len(indices),2),dtype=PBasicIndexType)
            planesize = (self.__dimensions[1]-1)
            res[:,0] = indices // planesize
            res[:,1] = indices - res[:,0]*(planesize)
            return res

    def GetMultiIndexOfElement(self,index):
        return self.GetMultiIndexOfElements([index])[0,:]

    def GetNumberOfElements(self):
        return np.prod((self.__dimensions-1)[self.__dimensions>=1] )


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

    def tighten(self):
        self.tags.Tighten()

@froze_it
class ConstantRectilinearMesh(MeshBase):

    def IsConstantRectilinear(self):
        return True

    def __init__(self,dim = 3):
        super(ConstantRectilinearMesh,self).__init__()
        #Number of nodes
        self.__dimensions = np.ones((dim,),dtype=PBasicIndexType)*2
        self.__origin = np.zeros((dim,) )
        self.__spacing = np.ones((dim,))
        self.nodes = None
        self.originalIDNodes = None
        self.elements = AllElements()
        self.structElements = ConstantRectilinearElementContainer(self.__dimensions)
        self.elements[self.structElements.elementType] = self.structElements

    def __copy__(self):
        res = ConstantRectilinearMesh(dim = len(self.__dimensions) )
        res._assign(self)
        res.__dimensions = self.__dimensions
        res.__origin = self.__origin
        res.__spacing = self.__spacing
        res.nodes = self.nodes
        res.originalIDNodes = self.originalIDNodes
        res.elements = self.elements
        res.structElements = self.structElements
        return res

    def GetElementsOriginalIDs(self,dim = None):
        """
        return a single list with all the originalid concatenated
        """
        res = np.empty(self.GetNumberOfElements(dim=dim),dtype=PBasicIndexType)
        cpt = 0
        from BasicTools.Containers.Filters import ElementFilter
        for name,data,ids in ElementFilter(self,dimensionality = dim):
            res[0+cpt:len(ids)+cpt] = data.originalIds[ids]
            cpt += len(ids)
        return res

    def GetNamesOfElemTagsBulk(self):
        return [ tag.name for tag in self.structElements.tags]

    def GetElementsInTagBulk(self,tagname):
        return self.structElements.tags[tagname].GetIds()

    def SetDimensions(self,data):
        self.__dimensions = np.array(data,int)
        self.structElements.SetDimensions(self.__dimensions)
        self.nodes = None
        self.originalIDNodes = None

    def GetDimensions(self):
        return np.array(self.__dimensions)

    def SetSpacing(self,data):
        self.__spacing = np.array(data, float)
        self.nodes = None

    def GetSpacing(self):
        return self.__spacing

    def GetdV(self):
        """Get the volume of one element."""
        return np.prod(self.GetSpacing())

    def SetOrigin(self,data):
        self.__origin = np.array(data)
        self.nodes = None

    def GetOrigin(self):
        return self.__origin

    @property
    def boundingMin(self):
        return self.GetOrigin()

    @property
    def boundingMax(self):
        return self.GetOrigin() + (self.GetDimensions()-1)*self.GetSpacing()

    def GetNumberOfNodes(self):
        return np.prod(self.__dimensions)

    def GetNumberOfElements(self,dim=None):
        """
        Compute and return the total number of elements in the mesh
        """

        numberOfElements = 0
        if dim == None:
            for elemname, data in self.elements.items():
                numberOfElements += data.GetNumberOfElements()
        else:
            for elemname, data in self.elements.items():
                if ElementNames.dimension[elemname] == dim:
                    numberOfElements += data.GetNumberOfElements()

        return numberOfElements

    def GetMultiIndexOfElements(self,indices):
        return self.structElements.GetMultiIndexOfElements(indices)

    def GetMultiIndexOfElement(self,index):
        return self.structElements.GetMultiIndexOfElement(index)

    def GetDimensionality(self):
        return len(self.__dimensions)

    def GetPointsDimensionality(self):
        return len(self.__dimensions)

    def GetMultiIndexOfNodes(self,indices):
        indices = np.asarray(indices,dtype=PBasicIndexType)

        if self.GetDimensionality() == 3:
            planesize = self.__dimensions[1] *self.__dimensions[2]
            res = np.empty((len(indices),3),dtype=PBasicIndexType)
            res[:,0] = indices // planesize
            resyz = indices - res[:,0]*(planesize)
            res[:,1] = resyz // self.__dimensions[2]
            res[:,2] =  resyz - res[:,1]*self.__dimensions[2]
            return res
        else:
            res = np.empty((len(indices),2),dtype=PBasicIndexType)
            res[:,0] = indices // self.__dimensions[1]
            res[:,1] = indices - res[:,0]*(self.__dimensions[1])
            return res

    def GetMultiIndexOfNode(self,index):
        return self.GetMultiIndexOfNodes([index])[0,:]

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

    def GetMonoIndexOfElements(self,indices):
        indices = np.asarray(indices,dtype=PBasicIndexType)
        if self.GetDimensionality() == 3:
            planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)
            return indices[:,0]*planesize+indices[:,1]*(self.__dimensions[2]-1) +indices[:,2]
        else :
            planesize = (self.__dimensions[1]-1)
            return indices[:,0]*planesize+indices[:,1]

    def GetMonoIndexOfElement(self,index):
        return self.GetMonoIndexOfElements([index])[0]

    def GetPosOfNode(self,index):
        if self.nodes is not None:
            return self.nodes[index,:]

        if self.GetDimensionality() == 3 :
            nxnynz = self.GetMultiIndexOfNode(index)
            return np.multiply(nxnynz,self.__spacing)+self.__origin
        else:
            nxny = self.GetMultiIndexOfNode(index)
            return np.multiply(nxny,self.__spacing)+self.__origin

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
              self.nodes = np.empty((self.GetNumberOfNodes(),2),dtype=PBasicFloatType)
              self.nodes[:,0] = xv.ravel()
              self.nodes[:,1] = yv.ravel()

              self.originalIDNodes = np.arange(self.GetNumberOfNodes(),dtype=PBasicIndexType)
              return self.nodes

            z = np.arange(self.__dimensions[2])*self.__spacing[2]+self.__origin[2]
            xv, yv, zv = np.meshgrid(x, y,z,indexing='ij')

            self.nodes = np.empty((self.GetNumberOfNodes(),3),dtype=PBasicFloatType)
            self.nodes[:,0] = xv.ravel()
            self.nodes[:,1] = yv.ravel()
            self.nodes[:,2] = zv.ravel()

            self.originalIDNodes = np.arange(self.GetNumberOfNodes(),dtype=PBasicIndexType)

        return self.nodes


    def GetElementsInTag(self,tagname,useOriginalId=False) :
        """
        return a list with the ids of the elements in a tag (only for the structElements)
        """
        if tagname in self.structElements.tags:
            return self.structElements.tags[tagname].GetIds()
        return np.zeros((0,),dtype=PBasicIndexType)

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
                           d2[0]*d2[2]*2,dtype=PBasicIndexType)


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
                           d2[1]*2,dtype=PBasicIndexType)


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
        pos = np.asarray(pos)
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

    # Error checking tests
    try:
        # not implemented for dim = 1 this line must fail
        myMesh = ConstantRectilinearMesh( dim=2)
        myMesh.SetDimensions([1,2,3])
        raise("Error detecting bad argument") # pragma: no cover
    except:
        pass


    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([1,1,1])
    myMesh.SetSpacing([1, 1, 1])

    try:
        # not implemented for dim = 1 this line must fail
        myMesh.GetNodalIndicesOfBorder()
        raise("Error detecting bad mesh props")# pragma: no cover
    except:
        pass


    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2,2])
    myMesh.SetSpacing([1, 1, 1])
    myMesh.ComputeGlobalOffset()
    myMesh.nodeFields["simpleNF"] = np.arange(myMesh.GetNumberOfNodes())
    myMesh.elemFields["simpleEF"] = np.arange(myMesh.GetNumberOfElements())
    myMesh.structElements.tags.CreateTag("First Element").SetIds([1])
    #myMesh.SetOrigin([-2.5,-1.2,-1.5])


    print(myMesh)
    if myMesh.GetNumberOfElements() != 1:
        raise Exception("Wrong number of elements")# pragma: no cover

    myMesh.structElements.tighten()

    import copy
    myNewMesh = copy.copy(myMesh)

    if myMesh.GetNumberOfElements() != myNewMesh.GetNumberOfElements():
        raise Exception("Error int the copy") # pragma: no cover

    print(myMesh.GetElementsOriginalIDs())
    print(myMesh.GetNamesOfElemTagsBulk())
    print(myMesh.GetdV())
    print(myMesh.boundingMin)
    print(myMesh.boundingMax)
    if myMesh.GetPointsDimensionality() != 3:
        raise Exception("Wrong dim of points") # pragma: no cover

    myMesh.GetElementsInTag("")
    myMesh.GetElementsInTag("First Element")


    myMesh.GetValueAtPosMultipleField(np.zeros((3,myMesh.GetNumberOfNodes())), [0,0,0])

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
    myMesh.SetDimensions([3,3])
    myMesh.SetSpacing([1, 1])
    myMesh.SetOrigin([0.,0.])

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


    print(myMesh.GetMultiIndexOfElement(0))
    print(myMesh.GetMultiIndexOfElement(0))
    print(myMesh.GetMultiIndexOfElements([0]))
    print(myMesh.GetMultiIndexOfElements(np.array([0,1])))
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
