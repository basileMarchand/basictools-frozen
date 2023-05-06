# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Optional, List
import numpy as np

from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType, ArrayLike
from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tags, Tag
from BasicTools.Containers.UnstructuredMesh import AllElements as AllElements
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject, froze_it
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1

@froze_it
class ConstantRectilinearElementContainer(BaseOutputObject):
    """Element container of a topologically regular array of elements.
    This class can hold regular element in 2D and 3D. Pixels and Voxel.
    The grid is defined by the dimensions (number of nodes in each direction)
    """

    def __init__(self, __dimensions: ArrayLike):
        """init this class

        Parameters
        ----------
        __dimensions : ArrayLike
            the number of nodes in each direction.
            must contains only 2 or 3  (for 2D or 3D grids)
        """
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
    def connectivity(self)-> np.ndarray:
        """Generate and retrieve the connectivity of the elements

        Returns
        -------
        np.ndarray
            Connectivity of the elements size (self.GetNumberOfElements(), 4 or 8)
        """
        if(self._connectivity is None):
            self._connectivity = self.GetConnectivityForElements(np.arange(self.GetNumberOfElements()))
            self._connectivity.flags.writeable = False
        return self._connectivity

    def SetDimensions(self, data: ArrayLike):
        """Set the number of points for the grid in each dimension

        Parameters
        ----------
        data : ArrayLike
            the number of points in each dimension for this grid
        """
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

    def GetDimensionality(self) -> int:
        """Return the dimensionality (2 for 2D, 3 for 3D) of this container

        Returns
        -------
        int
            the dimensionality of this containers
        """
        return len(self.__dimensions)

    def GetConnectivityForElements(self, indices:ArrayLike) -> np.ndarray:
        """Return the connectivity for the element listed in the indices array.
        If the connectivity is used many times, a simple call of self.connectivity will generate the
        connectivity for all the elements and keep it for later use.

        Parameters
        ----------
        indices : ArrayLike
            the indices of the elements to generate the connectivity

        Returns
        -------
        np.ndarray
            the connectivity of the selected elements size = (len(indices), 4 or 8 )
        """

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

    def GetConnectivityForElement(self, index: int )-> np.ndarray:
        """Get the connectivity matrix for one element
        please see documentation for GetConnectivityForElements

        Parameters
        ----------
        index : int
            element index

        Returns
        -------
        np.ndarray
            connectivity matrix
        """
        return self.GetConnectivityForElements([index])[0,:]

    def GetMultiIndexOfElements(self, indices:ArrayLike )-> np.ndarray:
        """Return the multi-index ijk for element with indices

        Parameters
        ----------
        indices : ArrayLike
            the indices of the elements to treat

        Returns
        -------
        np.ndarray
            an array with the ijk indices for every element in indices
            size (nb element, 2 (ij) or 3 (ijk) )
        """
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

    def GetMultiIndexOfElement(self, index : int ) -> np.ndarray:
        """Return the multi-index ijk for an element
        please see documentation for GetMultiIndexOfElements

        Parameters
        ----------
        index : int
            the index of the element

        Returns
        -------
        np.ndarray
            ij or ijk for the element
        """

        return self.GetMultiIndexOfElements([index])[0,:]

    def GetNumberOfElements(self)-> int :
        """Return the number of element in this container

        Returns
        -------
        int
            the number of elements
        """
        return np.prod((self.__dimensions-1)[self.__dimensions>=1] )


    def GetNumberOfNodesPerElement(self)-> int:
        """Return the number of node per element

        Returns
        -------
        int
            number of nodes
        """

        return 2**len(self.__dimensions)

    def GetTag(self, tagName: str) -> Tag:
        """Return the Tag by the name
        if the tag does not exist a new tag is created

        Parameters
        ----------
        tagName : str
            The name of the tag

        Returns
        -------
        Tag
            an instance of type Tag
        """
        return self.tags.CreateTag(tagName,False)

    def __str__(self)->str:
        """return a string of a summary of this container

        Returns
        -------
        str
            The description
        """
        res  = "    ConstantRectilinearElementContainer, "
        res += "  Type : ({},{}), ".format(self.elementType,self.GetNumberOfElements())
        res += "  Tags : " + " ".join([ ("("+x.name+":"+str(len(x)) +")") for x in self.tags]) + "\n"
        return res

    def tighten(self) -> None:
        """Tighten all the tags (free unused memory)
        Call Tag.Tighten on every tag
        """
        self.tags.Tighten()

@froze_it
class ConstantRectilinearMesh(MeshBase):
    """Topologically and geometrically regular array of data.
    This class can hold regular data in 2D and 3D. The positions
    of the nodes is generated using the origin, spacing and dimension
    of the mesh

    """

    def IsConstantRectilinear(self) -> bool:
        return True

    def __init__(self,dim = 3):
        """initialization of this class, the dim (dimensionality of the mesh)
        can't be changed after

        Parameters
        ----------
        dim : int, optional
            the dimensionality of the mesh (2 or 3), by default 3
        """
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

    def __copy__(self) -> ConstantRectilinearMesh:
        """Copy operator, the internal data is not copied.

        Returns
        -------
        ConstantRectilinearMesh
            a new instance of ConstantRectilinearMesh pointing to the same internal data
        """
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

    def GetElementsOriginalIDs(self, dim :Optional[int]= None) -> np.ndarray:
        """return a single list with all the originalid concatenated

        Parameters
        ----------
        dim : int, optional
            if dim != none generate the original ids only for the element of
            dimensionality dim, by default None

        Returns
        -------
        np.ndarray
            the original ids
        """
        res = np.empty(self.GetNumberOfElements(dim=dim),dtype=PBasicIndexType)
        cpt = 0
        from BasicTools.Containers.Filters import ElementFilter
        for name,data,ids in ElementFilter(self,dimensionality = dim):
            res[0+cpt:len(ids)+cpt] = data.originalIds[ids]
            cpt += len(ids)
        return res

    def GetNamesOfElemTagsBulk(self)-> List[str] :
        """Return a list of the tag names of the element of dimensionality == dim

        Returns
        -------
        List[str]
            The names of the tags
        """
        return [ tag.name for tag in self.structElements.tags]

    def GetElementsInTagBulk(self, tagname:str)->np.ndarray:
        """return the element in the tag of name tagname

        Parameters
        ----------
        tagname : str
            the name of the tag

        Returns
        -------
        np.ndarray
            the ids of the elements
        """
        return self.structElements.tags[tagname].GetIds()

    def SetDimensions(self, data:ArrayLike) -> None:
        """Set the number of nodes and element in each direction

        Parameters
        ----------
        data : ArrayLike
            Number of node in each direction
        """
        self.__dimensions = np.array(data,int)
        self.structElements.SetDimensions(self.__dimensions)
        self.nodes = None
        self.originalIDNodes = None

    def GetDimensions(self) -> np.ndarray:
        """Return the number of nodes in each direction

        Returns
        -------
        np.ndarray
            the number of nodes in each direction
        """
        return np.array(self.__dimensions)

    def SetSpacing(self, data:ArrayLike) -> None:
        """Set the spacing. The length of the element in each direction

        Parameters
        ----------
        data : ArrayLike
            size in each coordinate
        """
        self.__spacing = np.array(data, float)
        self.nodes = None

    def GetSpacing(self) -> np.ndarray:
        """Return the spacing. The length of the elements in each direction

        Returns
        -------
        np.ndarray
            size in each direction
        """
        return self.__spacing

    def GetdV(self) -> PBasicFloatType:
        """Get the volume of one element.

        Returns
        -------
        PBasicFloatType
            The volume of one element (this is the product of the spacing)
        """
        return np.prod(self.GetSpacing())

    def SetOrigin(self, data:ArrayLike) -> None:
        """Set the origin of this mesh.

        Parameters
        ----------
        data : ArrayLike
            the coordinates of the first point in this mesh
        """
        self.__origin = np.array(data)
        self.nodes = None

    def GetOrigin(self) -> np.ndarray:
        """Return the origin. the coordinate of the first node of the mesh

        Returns
        -------
        np.ndarray
            the coordinate position of the first point
        """
        return self.__origin

    @property
    def boundingMin(self) -> np.ndarray:
        """the bounding box minimum of the mesh

        Returns
        -------
        np.ndarray
            the coordinate of the lower corner of the bounding box
        """
        return self.GetOrigin()

    @property
    def boundingMax(self)-> np.ndarray:
        """The bounding box maximum of the mesh

        Returns
        -------
        np.ndarray
            the coordinates of the higher corner of the bounding box
        """
        return self.GetOrigin() + (self.GetDimensions()-1)*self.GetSpacing()

    def GetNumberOfNodes(self)-> PBasicIndexType:
        """Return the number of nodes in the mes

        Returns
        -------
        PBasicIndexType
            The number of nodes in the mesh
        """
        return np.prod(self.__dimensions)

    def GetNumberOfElements(self, dim: int=None)-> PBasicIndexType:
        """Compute and return the total number of elements in the mesh

        Parameters
        ----------
        dim : int, optional
            the user can filter by the dimensionality, by default None

        Returns
        -------
        PBasicIndexType
            number of element in the mesh
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

    def GetMultiIndexOfElements(self, indices:ArrayLike) -> np.ndarray:
        """Please see ConstantRectilinearElementContainer.GetMultiIndexOfElements Documentation

        Parameters
        ----------
        indices : ArrayLike
            indices

        Returns
        -------
        np.ndarray
            MultiInndex
        """
        return self.structElements.GetMultiIndexOfElements(indices)

    def GetMultiIndexOfElement(self, index: int ) -> np.ndarray:
        """Please see ConstantRectilinearElementContainer.GetMultiIndexOfElement Documentation

        Parameters
        ----------
        index : int
            the index of the element

        Returns
        -------
        np.ndarray
             ij or ijk for the element
        """
        return self.structElements.GetMultiIndexOfElement(index)

    def GetDimensionality(self)-> int:
        """Get the dimensionality of the mesh 2 for 2D or 3 for 3D

        Returns
        -------
        int
            The dimensionality of the mesh
        """
        return len(self.__dimensions)

    def GetPointsDimensionality(self) -> int:
        """Get the point dimensionality of the mesh 2 for 2D or 3 for 3D

        Returns
        -------
        int
            The dimensionality of the mesh
        """
        return len(self.__dimensions)

    def GetMultiIndexOfNodes(self, indices:ArrayLike) -> np.ndarray:
        """Return the multi-index ijk for nodes with indices

        Parameters
        ----------
        indices : ArrayLike
            the indices of the nodes to treat

        Returns
        -------
        np.ndarray
            an array with the ijk indices for every point in indices
            size (nb points, 2 (ij) or 3 (ijk) )
        """
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

    def GetMultiIndexOfNode(self, index: int) -> np.ndarray:
        """Return the multi-index ijk for a node
        Please see GetMultiIndexOfNodes.ConstantRectilinearMesh

        Parameters
        ----------
        index : int
            the index of the node

        Returns
        -------
        np.ndarray
            ij or ijk for the node
        """
        return self.GetMultiIndexOfNodes([index])[0,:]

    def GetMonoIndexOfNode(self, indices: ArrayLike)-> np.ndarray:
        """return the mono index for the nodes with multi index indices

        Parameters
        ----------
        indices : ArrayLike
            the ij or ijk indices of the nodes

        Returns
        -------
        np.ndarray
            the mono index for the nodes
        """
        indices= np.asarray(indices)
        if len(indices.shape) == 1:
            indexs = indices[np.newaxis]
        else:
            indexs = indices

        if self.GetDimensionality() == 3:

            planesize = self.__dimensions[1] *self.__dimensions[2]
            res = planesize*indexs[:,0]
            res += indexs[:,1]*self.__dimensions[2]
            res += indexs[:,2]
            return res
        else:
            planesize = self.__dimensions[1]
            return planesize*indexs[:,0]+indexs[:,1]

    def GetMonoIndexOfElements(self, indices: ArrayLike)-> np.ndarray:
        """return the mono index for the elements with multi index indices

        Parameters
        ----------
        indices : ArrayLike
            the ij or ijk indices of the elements

        Returns
        -------
        np.ndarray
            the mono index for the element
        """
        indices = np.asarray(indices,dtype=PBasicIndexType)
        if self.GetDimensionality() == 3:
            planesize = (self.__dimensions[1]-1) *(self.__dimensions[2]-1)
            return indices[:,0]*planesize+indices[:,1]*(self.__dimensions[2]-1) +indices[:,2]
        else :
            planesize = (self.__dimensions[1]-1)
            return indices[:,0]*planesize+indices[:,1]

    def GetMonoIndexOfElement(self, indices:int)-> np.ndarray:
        """Return the mono index for the element with multi index indices


        Parameters
        ----------
        indices : int
            the multi index of the element

        Returns
        -------
        np.ndarray
            the mono index
        """

        return self.GetMonoIndexOfElements([indices])[0]

    def GetPosOfNode(self, indices:ArrayLike) -> np.ndarray:
        """Return the position of the nodes with indices

        Parameters
        ----------
        indices : ArrayLike
            Nodes to treat

        Returns
        -------
        np.ndarray
            Position of the nodes
        """
        if self.nodes is not None:
            return self.nodes[indices,:]

        if self.GetDimensionality() == 3 :
            nxnynz = self.GetMultiIndexOfNode(indices)
            return np.multiply(nxnynz,self.__spacing)+self.__origin
        else:
            nxny = self.GetMultiIndexOfNode(indices)
            return np.multiply(nxny,self.__spacing)+self.__origin

    def GetPosOfNodes(self):
        """
        position for all nodes in the mesh.

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


    def GetElementsInTag(self, tagname:str, useOriginalId:bool=False) -> np.ndarray:
        """return a list with the ids of the elements in a tag (only for the structElements)


        Parameters
        ----------
        tagname : str
            the name of the tag
        useOriginalId : bool, optional
            if True return the original ids, by default False

        Returns
        -------
        np.ndarray
            the id or the original ids
        """
        if tagname in self.structElements.tags:
            return self.structElements.tags[tagname].GetIds()
        return np.zeros((0,),dtype=PBasicIndexType)

    def GetNodalIndicesOfBorder(self, border:int =0)-> np.ndarray:
        """Return the ids of the nodes in the border (first layer if border is 0)
        the id of the second line of nodes (second layer if border is 1)
        ...

        Parameters
        ----------
        border : int, optional
            the layer number to extract, by default 0

        Returns
        -------
        np.ndarray
            the ids of the nodes in the layer = border
        """

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

    def GetClosestPointToPos(self, pos:ArrayLike, MultiIndex:bool=False) -> np.ndarray:
        """Get the index of the closes point to pos

        Parameters
        ----------
        pos : ArrayLike
            Position of the probe
        MultiIndex : bool, optional
            if MultiIndex is True the, the multi-index is returned, by default False

        Returns
        -------
        np.ndarray
            the index of the closest node to pos
        """
        pos = (pos-self.__origin)-self.__spacing/2.
        pos /= self.__spacing
        elemindex = pos.astype(int)
        elemindex = np.minimum(elemindex, self.__dimensions-1)
        elemindex = np.maximum(elemindex, 0)
        if MultiIndex :
            return elemindex

        return self.GetMonoIndexOfNode(elemindex)


    def GetElementAtPos(self, pos: ArrayLike, MultiIndex:bool=False) -> np.ndarray:
        """Get the index of the closes element to pos

        Parameters
        ----------
        pos : ArrayLike
            Position of the probe
        MultiIndex : bool, optional
            if MultiIndex is True the, the multi-index is returned, by default False

        Returns
        -------
        np.ndarray
            the index of the closest element to pos
        """
        pos = pos-self.__origin
        pos /= self.__spacing
        elemindex = pos.astype(int)
        elemindex = np.minimum(elemindex, self.__dimensions-2)
        elemindex = np.maximum(elemindex, 0)
        if MultiIndex :
            return elemindex

        return self.GetMonoIndexOfElement(elemindex)

    def GetElementShapeFunctionsAtPos(self, el:int, pos:ArrayLike) -> np.ndarray:
        """Get the shape function values at position pos

        Parameters
        ----------
        el : int
            the id of the element to evaluate the shape function
        pos : ArrayLike
            Real space coordinate of the position to evaluate the shape functions

        Returns
        -------
        np.ndarray
            the values of the shape function at pos
        """
        coon = self.GetConnectivityForElement(el)
        p0 = self.GetPosOfNode(coon[0])
        n0 = (pos-p0)*2./self.__spacing - 1.

        return self.structElements.space.GetShapeFunc(n0)

    def GetValueAtPos(self, field:ArrayLike, pos:ArrayLike)-> PBasicFloatType:
        """Evaluate a point field (field defined at nodes) at the position pos

        Parameters
        ----------
        field : ArrayLike
            the values at every node of the mesh
        pos : ArrayLike
            the position of the probe

        Returns
        -------
        PBasicFloatType
            the value of the field at pos
        """
        el = self.GetElementAtPos(pos)
        coon = self.GetConnectivityForElement(el)
        xiChiEta = self.GetElementShapeFunctionsAtPos(el,pos)
        return field[coon].dot(xiChiEta)


    def GetConnectivityForElement(self, index:int)->np.ndarray:
        return self.structElements.GetConnectivityForElement(index)

    def GenerateFullConnectivity(self)-> np.ndarray:
        return self.structElements.connectivity

    def ComputeGlobalOffset(self)->None:
        """
        Recompute the Global Offset,
        This is necessary for some operations.
        Recommendation : Call it after changing the topology
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



    np.set_printoptions(threshold=sys.maxsize)


    print((myMesh.GetValueAtPos(np.array([1,2,3,4,5,6,7,8,9]),[0.5,0.5])))


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
