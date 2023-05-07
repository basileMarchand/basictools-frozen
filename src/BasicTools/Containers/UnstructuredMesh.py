# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Dict, Union

import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType, ArrayLike
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tags, Tag

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject, froze_it

@froze_it
class ElementsContainer(BaseOutputObject):
    """
    Class to hold a list of element of the same type

    * elementType : a string from BasicTools.Containers.ElementNames
    * connectivity : the connectivity matrix starting form 0
    * tags : the tags holder class
    * originalIds : the id or number from the previous mesh/file

    The user can use this data to find the mapping from the initial mesh/file
    to the current mesh (self).

    * self.globaloffset  : this value is calculate automatically by the mesh
    .. deprecated:: 2.0.0

    * self.cpt : an internal counter to do efficient add of elements one by one

    The user is responsible to call self.tighten() to compact the connectivity
    matrix after the population ( calls AddNewElement(...) or allocate(...))
    """
    def __init__(self,elementType):
        super(ElementsContainer,self).__init__()
        self.elementType = elementType
        self.connectivity = np.empty((0,ElementNames.numberOfNodes[elementType]),dtype=PBasicIndexType)
        self.globaloffset   = 0
        self.tags = Tags()
        self.cpt = 0

        self.originalIds = np.empty((0,),dtype=PBasicIndexType)
        self.originalOffset = 0
        self.mutable = True

    def __eq__(self, other:ElementsContainer) -> bool:
        """return True if 2 container are equal

        Parameters
        ----------
        other : ElementsContainer
            _description_

        Returns
        -------
        bool
            True if the 2 containers are equal, False ir not
        """

        if self.elementType != other.elementType:
            return False

        self.tighten()
        other.tighten()

        if not np.array_equal(self.connectivity, other.connectivity):
            return False

        if not np.array_equal(self.originalIds, other.originalIds):
            return False

        if self.tags != other.tags:
            return False

        if self.originalOffset != other.originalOffset:
            return False

        return True

    def GetNumberOfElements(self)->int :
        """return the number of elements in this container

        Returns
        -------
        int
            Number of elements
        """
        return self.cpt
        #return self.connectivity.shape[0]

    def Merge(self,other:ElementsContainer, offset:int=None):
        """
        Merge the elements from the other container into this.

        Warning!!! Non elimination of double elements is done

        if an offset is supplied the connectivity of the other container is
        shifted by the value of the offset during the merge
        """
        other.tighten()
        if other.cpt == 0:
            return

        self.Reserve(self.cpt+other.cpt)

        if offset is None:
            offset = 0

        self.connectivity[self.cpt:,:] = other.connectivity+offset
        self.originalIds[self.cpt:] = -1*np.arange(other.cpt)

        for tag in other.tags:
            self.GetTag(tag.name).AddToTag(tag.GetIds() + self.cpt)

        self.cpt += other.cpt

    def AddNewElements(self, conn: ArrayLike, originalids: ArrayLike=None) -> int:
        """Append a new elements to the connectivity

        Parameters
        ----------
        conn : ArrayLike
            connectivity of the added elements
        originalids : ArrayLike, optional
            the original id of the added element, by default -1 is used

        Returns
        -------
        int
            the total number of elements in the container
        """
        onoe = self.GetNumberOfElements()
        self.Allocate(onoe+conn.shape[0])

        self.connectivity[onoe:onoe+conn.shape[0],:] = conn

        if originalids is None:
            self.originalIds[onoe:onoe+conn.shape[0]] = -1
        else:
            self.originalIds[onoe:onoe+conn.shape[0]] = originalids

        return self.cpt

    def AddNewElement(self, conn: ArrayLike, originalid:int)-> int:
        """append a single element to the connectivity

        Parameters
        ----------
        conn : ArrayLike
            connectivity of the added element
        originalid : int
            the original id of the added element

        Returns
        -------
        int
            the total number of elements in the container
        """
        if self.cpt >= self.connectivity.shape[0]:
            self.Reserve(2*self.cpt+1)

        self.connectivity[self.cpt,:] = conn
        self.originalIds[self.cpt] = originalid
        self.cpt +=1

        return self.cpt

    def GetNumberOfNodesPerElement(self) -> int :
        """return the number of nodes per element for the elements in this container

        Returns
        -------
        int
            the number of nodes per element for the elements in this container
        """
        if  self.connectivity.shape[1]:
            return self.connectivity.shape[1]
        return ElementNames.numberOfNodes[self.elementType]

    def GetNodesIdFor(self, ids:Union[int,ArrayLike] ) -> np.ndarray:
        """return the nodes used by the list of elements

        input:
            ids :

        Parameters
        ----------
        ids : Union[int,ArrayLike]
            the id or a list of ids of elements to treat (always a local id list)

        Returns
        -------
        np.ndarray
            the ids of the nodes used by the elements ids
        """

        return np.unique(self.connectivity[ids,:])

    def GetTag(self, tagName:str)->Tag:
        """ return the tag by name.
        If the tag does not exist a new tag is created

        Parameters
        ----------
        tagName : str
            the name of the tag to retrieve

        Returns
        -------
        Tag
            a tag with the name tagName
        """
        return self.tags.CreateTag(tagName,False)

    def Reserve(self,nbElements:int):
        """Reserve the storage for nbElements

        The user is responsible to call self.tighten() to compact the connectivity
        and the originalIds.
        matrix after the population

        Parameters
        ----------
        nbElements : int
            Number of element to reserve
        """

        if nbElements != self.connectivity.shape[0]:
            self.connectivity =  np.resize(self.connectivity, (nbElements,self.GetNumberOfNodesPerElement()))
            self.originalIds =  np.resize(self.originalIds, (nbElements,))

    def Allocate(self,nbElements:int):
        """Allocate the storage for nbElements

        the user is responsible of filling the connectivity and the originalIds
        with valid values

        Parameters
        ----------
        nbElements : int
            Number of element to Allocate
        """
        self.Reserve(nbElements)
        self.cpt = nbElements

    def tighten(self)-> None:
        """Compact the storage and free non used space.
        """
        self.Reserve(self.cpt)
        self.tags.Tighten()

    def AddElementToTag(self, globalElemNumber:int, tagname:str)-> bool:

        """Add an element to a tag using a global element number, if the global number
        The user must compute the globaloffset first to make this function work

        .. deprecated:: 2.0.0
       This function is to complex and to many side effects.

        Returns
        -------
        bool
            True if the element is added to the tag (in the case globalElemNumber element is in this
             ElementContainers ).
        """
        if globalElemNumber -self.globaloffset <  self.GetNumberOfElements():
            self.tags.CreateTag(tagname,False).AddToTag(globalElemNumber-self.globaloffset)
            return True
        else:
            return False

    def __str__(self):
        res  = "    ElementsContainer, "
        res += "  Type : ({},{}), ".format(self.elementType,self.GetNumberOfElements())
        res += "  Tags : " + " ".join([ ("("+x.name+":"+str(len(x)) +")") for x in self.tags]) + "\n"
        return res

@froze_it
class AllElements(object):
    """
    Class to store a list of element containers. This class is mostly a re-implementation of dict
    but with ordered keys.
    This class is sorted by keys, in lexicographic order, so the retrieving order is stable.

    note:
        FB: the number of different types of elements is low, I don't think
        this is gonna add alot of overhead to the library
    """

    def __init__(self):
        super(AllElements,self).__init__()
        self.storage = {}

    def __eq__(self, other):
        if len(self.storage) != len(other.storage):
            return False

        for i in self:
            if i in other:
                if self[i] != other[i]:
                    return False
            else:
                return False
        return True

    def keys(self):
        return sorted(self.storage.keys())

    def values(self):
        return [ values for key,values in sorted(self.storage.items()) ]

    def __iter__(self):
        return iter(self.keys())

    def items(self):
        return  sorted(self.storage.items())

    #send basis functions calls to the storage dictionary
    def __setitem__(self, key, value):
        self.storage[key] = value

    def __len__(self):
        return len(self.storage)

    def __contains__(self, k):
        return k in self.storage

    def __getitem__(self,key):
        return self.storage[key]

    def __delitem__(self,key):
        del self.storage[key]

    def GetElementsOfType(self, typename:str) -> ElementsContainer:
        """return a ElementContainer by name

        Parameters
        ----------
        typename : str
            the name of the elements to retrieve

        Returns
        -------
        ElementsContainer
            the container
        """
        if not typename in self:
            self[typename] = ElementsContainer(typename)
        return self[typename]

    def GetTagsNames(self)-> List[str]:
        """return a list with all the tags in this container

        Returns
        -------
        List[str]
            the list of the tag names
        """
        res = set()
        for data in self.values():
            res.update(data.tags.keys() )

        return list(res)

    def __str__(self):
        res = ""
        for data in self.storage.values():
            res += str(data)
        return res

@froze_it
class UnstructuredMesh(MeshBase):
    """
    Class to store an unstructured (i.e. general) mesh:

    * self.nodes : the points positions
    * self.originalIdNodes : the ids of the previous mesh/file
    * self.elements : the list of all the element in the mesh
    * self.boundingMin/Max : the bounding box of the mesh (use ComputeBoundingBox to compute it)

    The manual construction of this class must always end with a call to the
    function PrepareForOutput or with the context manager mesh.WithModification().
    """
    def IsUnstructured(self):
        return True

    def __init__(self):
        super(UnstructuredMesh, self).__init__()
        self.nodes = np.empty((0,3),dtype=PBasicFloatType)
        self.originalIDNodes = np.empty((0,),dtype=PBasicIndexType)
        self.elements = AllElements()
        self.boundingMin = np.array([0.,0,0])
        self.boundingMax = np.array([0.,0,0])

    def __copy__(self):
        res = UnstructuredMesh()
        res._assign(self)
        res.nodes = self.nodes
        res.originalIDNodes = self.originalIDNodes
        res.elements = self.elements
        res.boundingMin = self.boundingMin
        res.boundingMax = self.boundingMax
        return res

    def __eq__(self, other):
        if not super(UnstructuredMesh,self).__eq__(other):
            return False

        if not np.array_equal(self.nodes,other.nodes):
            return False

        if not np.array_equal(self.originalIDNodes, other.originalIDNodes):
            return False

        if self.elements != other.elements:
            return False

        return True

    def ConvertDataForNativeTreatment(self):
        """Ensure all the data is compatible with the cpp treatment (continuous in C order)
        """
        self.originalIDNodes = np.asarray(self.originalIDNodes, dtype=PBasicIndexType, order="C")
        self.nodes = np.asarray(self.nodes, dtype=PBasicFloatType, order="C")
        for data in self.elements.values():
            data.connectivity = np.asarray(data.connectivity, dtype=PBasicIndexType, order="C")

    def GetNumberOfNodes(self) -> int:
        """return the total number of nodes in the mesh

        Returns
        -------
        int
            Number of node in the mesh
        """
        return self.nodes.shape[0]

    def SetNodes(self, nodes:ArrayLike, originalIDNodes:ArrayLike=None, generateOriginalIDs:bool = False):
        """Set nodes and original Ids in the correct internal format.

        Parameters
        ----------
        nodes : ArrayLike
            Nodes positions of size (number of node, space dimensionality )
        originalIDNodes : ArrayLike, optional
            original Ids, by default None
        generateOriginalIDs : bool, optional
            if True to generate original ids with numpy.arange, by default False
        """
        self.nodes = np.require(nodes, dtype=PBasicFloatType, requirements=['C','A'])
        if originalIDNodes is not None:
            self.originalIDNodes = np.require(originalIDNodes,dtype=PBasicIndexType,requirements=['C','A'])
        elif generateOriginalIDs:
            self.originalIDNodes = np.arange(self.GetNumberOfNodes(), dtype=PBasicIndexType)

    def GetPointsDimensionality(self) -> int:
        """Return the number of coordinates of the points

        Returns
        -------
        int
            number of columns in the point array
        """
        return self.nodes.shape[1]

    def GetDimensionality(self):
        """
        .. deprecated:: 2.0.0
        Please use GetPointsDimensionality()

        return the dimensionality 2/3
        """
        return self.nodes.shape[1]

    def GetNumberOfElements(self, dim : int= None) -> int:
        """Compute and return the total number of elements in the mesh

        Parameters
        ----------
        dim : int, optional
            dimensionality filter, by default None

        Returns
        -------
        int
            number of element in the mesh (filtered by dim if dim != None)
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

    def MergeElements(self, other:UnstructuredMesh, force:bool=False):
        """Merge the elements for a second mesh into this
        the nodes array must be the same (not only equal)

        the user can force the merge if needed (force variable)

        Parameters
        ----------
        other : UnstructuredMesh
            the other mesh
        force : bool, optional
            True to force the merge even if the nodes are not the same, by default False
        """
        if (self.nodes is not other.nodes) and (not force) :
            raise(RuntimeError("the two meshes does not share the same nodes fields (potentially dangerous)"))

        for name,data in other.elements.items():
            self.GetElementsOfType(name).Merge(data)

    def ComputeGlobalOffset(self) -> Dict[str,PBasicIndexType]:
        """
        Recompute the Global Offset,
        This is necessary for some operation.
        Recomendation : Call it after changing the topology

        .. deprecated:: 2.0.0
        The globaloffset variable will be deprecated
        """
        offsets = dict()
        cpt = 0
        for data in self.elements.values():
            offsets[data.elementType] = cpt
            data.globaloffset = cpt
            cpt += data.GetNumberOfElements()
        return offsets

    def ComputeBoundingBox(self)->None:
        """to recompute the bounding box (min and max)
        """
        self.boundingMin = np.amin(self.nodes, axis=0)
        self.boundingMax = np.amax(self.nodes, axis=0)

    def AddNodeToTagUsingOriginalId(self, oid:int, tagname:str)->None:
        """add a node (using the original id ) to a tag (tagname)

        Parameters
        ----------
        oid : int
            Original id node
        tagname : str
            Tag name

        Raises
        ------
        Exception
            if the original id is not found
        """
        w = np.where(self.originalIDNodes == oid)
        if len(w[0]) > 0 :
            self.GetNodalTag(tagname).AddToTag(w[0])
        else:
            raise Exception("No node with id " + str(oid)) #pragma: no cover

    def AddElementToTagUsingOriginalId(self, oid:int, tagname:str):
        """Add an element (using the originalId) to a tag (tagname)


        Parameters
        ----------
        oid : int
            the original id of the element
        tagname : str
            Tag name

        Raises
        ------
        Exception
            if the element with the original id is not found
        """
        for data in self.elements.values():
            w = np.where(data.originalIds[:data.cpt] == oid)
            if len(w[0]) > 0 :
                data.tags.CreateTag(tagname,False).AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid)) #pragma: no cover

    def AddElementsToTag(self, globalElemNumbers:ArrayLike, tagname:str):
        """Add elements (using the global element number) to a tag (tagname)
        you must compute the globaloffset first to make this function work

        Parameters
        ----------
        globalElemNumbers : ArrayLike
            the list of the global ids of the element
        tagname : str
            Tag name

        Raises
        ------
        Exception
            if some elementid are greater than the number of elements.
        """
        elementNotTreated = np.unique(globalElemNumbers)
        cpt = 0
        for data in self.elements.values():
            cpt2 = data.GetNumberOfElements() +cpt
            f = elementNotTreated < cpt2
            dataToTreat = elementNotTreated[f]
            if len(dataToTreat):
                data.tags.CreateTag(tagname,False).AddToTag(dataToTreat-data.globaloffset)
            elementNotTreated = elementNotTreated[np.logical_not(f)]
            cpt += data.GetNumberOfElements()

        if len(elementNotTreated):
            raise Exception("No element found") #pragma: no cover


    def AddElementToTag(self, globalElemNumber:int, tagname:str):
        """Add an element (using the global element number) to a tag (tagname)
        # you must compute the globaloffset first to make this function work

        Parameters
        ----------
        globalElemNumber : int
            the element number
        tagname : str
            Tag name

        Raises
        ------
        Exception
            if the the element is not found
        """
        for data in self.elements.values():
            if data.AddElementToTag(globalElemNumber,tagname):
                return
        raise Exception("No element found") #pragma: no cover

    def DeleteElemTags(self, tagNames: List[str]):
        """Delete element tags

        Parameters
        ----------
        tagNames : List[str]
            List of element tags to be deleted
        """
        #check not a string but a list like
        assert not isinstance(tagNames, str)
        for data in self.elements.values():
            data.tags.DeleteTags(tagNames)

    def GetPosOfNode(self, i:int ) -> np.ndarray:
        """Return the position of the node i

        Parameters
        ----------
        i : int
            id of the node

        Returns
        -------
        np.ndarray
            position of the point
        """
        return self.nodes[i,:]

    def GetPosOfNodes(self) -> np.ndarray:
        """return the position of all the nodes

        Returns
        -------
        np.ndarray
            The position for all the nodes in the mesh
        """
        return self.nodes

    def GetElementsOriginalIDs(self,dim:int = None)-> np.ndarray:
        """return a single list with all the originalid concatenated

        Parameters
        ----------
        dim : int, optional
            dimensionality filter, by default None

        Returns
        -------
        np.ndarray
            array with all the originalId for the selected elements
        """
        if dim is None:
            res = np.empty(self.GetNumberOfElements(),dtype=PBasicIndexType)
            cpt =0
            for name, data in self.elements.items():
                n = data.originalIds.shape[0]
                res[0+cpt:n+cpt] = data.originalIds
                cpt += n
        else:
            res = np.empty(self.GetNumberOfElements(dim=dim),dtype=PBasicIndexType)
            cpt = 0
            from BasicTools.Containers.Filters import ElementFilter
            for name, data, ids in ElementFilter(self,dimensionality = dim):
                res[0+cpt:len(ids)+cpt] = data.originalIds[ids]
                cpt += len(ids)
        return res

    def SetElementsOriginalIDs(self, originalIDs:ArrayLike):
        """Set from a single list all the element original Ids

        Parameters
        ----------
        originalIDs : ArrayLike
            the list of all the original ids in the order of the mesh
        """
        cpt = 0
        for data in self.elements.values():
            data.originalIds = originalIDs[cpt:data.GetNumberOfElements()+cpt]
            cpt += data.GetNumberOfElements()

    def GetElementsInTag(self, tagname:str, useOriginalId:bool=False) ->np.ndarray:
        """return a list with the ids of the elements in a tag

        Parameters
        ----------
        tagname : str
            the tag name
        useOriginalId : bool, optional
            to return the list of original ids and not the global numbers, by default False

        Returns
        -------
        np.ndarray
            the list off all the element in the tag named "tagname"
        """

        res = np.zeros((self.GetNumberOfElements(),),dtype=PBasicIndexType)
        cpt =0
        offsets = self.ComputeGlobalOffset()
        for name,data in self.elements.items():
            if tagname in data.tags:
                tag = data.tags[tagname].GetIds()
                if useOriginalId:
                    res[cpt:cpt+len(tag) ] = data.originalIds[tag]
                else:
                    res[cpt:cpt+len(tag) ] = offsets[name]+tag
                cpt +=  len(tag)
        return res[0:cpt]

    def PrepareForOutput(self) -> None:
        """Function to free the extra memory used after the incremental creation of a mesh
        and final treatment (offset computation).
        """
        self.nodesTags.Tighten()
        for data in self.elements.values():
            data.tags.RemoveDoubles()
            data.tighten()
        self.ComputeGlobalOffset()
        self.VerifyIntegrity()

    def __str__(self) -> str:
        """Return a string with a summary of the mesh

        Returns
        -------
        str
            summary of the mesh
        """
        res  = "UnstructuredMesh \n"
        res += "  Number Of Nodes    : {} \n".format(self.GetNumberOfNodes())
        res += "    Tags : " + " ".join( ["("+x.name+":"+str(len(x))+")" for x in  self.nodesTags ]) + "\n"

        res += "  Number Of Elements : {} \n".format(self.GetNumberOfElements())
        for data in self.elements.values():
            res += str(data)
        if len(self.nodeFields.keys()) > 0:
            res += "  nodeFields         : " + str(list(self.nodeFields.keys())) + "\n"
        if len(self.elemFields.keys()) > 0:
            res += "  elemFields         : " + str(list(self.elemFields.keys())) + "\n"
        return res

    def Clean(self)->None:
        """Remove unnecessary data :
            1) empty tags
            2) empty element containers
        """
        self.nodesTags.RemoveEmptyTags()
        for k in list(self.elements.keys()):
            if self.elements[k].GetNumberOfElements() == 0:
                del self.elements[k]
                continue
            self.elements[k].tags.RemoveEmptyTags()

    def VerifyIntegrity(self) -> None:
        """Integrity verification  of the mesh.,

        Raises
        ------
        Exception
            if the integrity of the mesh is compromised
        """
        #verification nodes an originalIdNodes are compatible
        if len(self.nodes.shape) !=2:
            raise Exception("Error in the shape of nodes")

        if self.nodes.flags['C_CONTIGUOUS'] is False:
            raise Exception("Error in the order of nodes")

        if len(self.originalIDNodes.shape) !=1:
            print(self.originalIDNodes.shape)
            raise Exception("Error in the shape of originalIDNodes")

        if self.originalIDNodes.flags['C_CONTIGUOUS'] is False:
            raise Exception("Error in the order of originalIDNodes")

        if self.originalIDNodes.shape[0] != self.nodes.shape[0]:
            print(self.originalIDNodes.shape)
            print(self.nodes.shape)
            raise Exception("nodes and originalIDNodes are incompatible")

        nbnodes = self.nodes.shape[0]

        #verification of min max in nodes tags
        for elemtype,data in self.nodesTags.items():
            ids = data.GetIds()
            if len(ids) == 0:
                continue
            if ids[0] < 0:
                raise Exception("Ids of '"+ elemtype +"' tag out of bound (<0)")
            if ids[-1] >= nbnodes:
                print(ids)
                print(nbnodes)
                raise Exception("Ids of '"+ elemtype +"' tag out of bound > nbnodes")

        # verification of elements
        for elemtype, data in self.elements.items():
            if data.connectivity.flags['C_CONTIGUOUS'] is False:
                raise Exception("Error :  connectivity not C_CONTIGUOUS")

            if len(data.connectivity.shape) != 2:
                raise Exception("Wrong shape of connetivity of elements '"+ elemtype +"'")

            if data.connectivity.shape is False:
                raise Exception("Error in the order of connecitivty")

            if data.originalIds.shape[0] != data.connectivity.shape[0]:
                print(data.originalIds.shape[0])
                print(data.connectivity.shape[0])
                raise Exception("connectivity and originalIds are incompatible '"+elemtype+"'")

            if data.originalIds.flags['C_CONTIGUOUS'] is False:
                raise Exception("Error in the order of originalIds")

            if ElementNames.numberOfNodes[elemtype] != data.connectivity.shape[1]:
                print(elemtype)
                print(ElementNames.numberOfNodes[elemtype])
                print(data.connectivity.shape[1])
                raise Exception("Incompatible number of columns of the connectivity array")

            if len(data.connectivity) and  np.amin(data.connectivity) < 0:
                raise Exception("connectivity of '"+elemtype+"' out of bound (<0)")

            if len(data.connectivity) and np.amax(data.connectivity) >= nbnodes:
                print(data.connectivity)
                print(nbnodes)
                raise Exception("connectivity of '"+elemtype+"' out of bound > nbnodes")

            nbelements = data.connectivity.shape[0]


            #verification of min max in nodes tags
            for tagName,tagData in data.tags.items():
                ids = tagData.GetIds()
                if len(ids) == 0:
                    continue
                if ids[0] < 0:
                    print(nbelements)
                    print(ids)
                    raise Exception("Ids of '"+ tagName +"' tag out of bound (<0)")
                if ids[-1] >= nbelements:
                    print(nbelements)
                    print(ids)
                    raise Exception("Ids of '"+ tagName +"' tag out of bound > nbelements")

def CheckIntegrity():
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,2,3],[0,1,0]], [[0,1,2]])
    res.ConvertDataForNativeTreatment()

    elements = res.GetElementsOfType(ElementNames.Triangle_3)

    elements = res.GetElementsOfType(ElementNames.Bar_2)
    elements.AddNewElement([1,2],1)

    elements.GetNumberOfNodesPerElement()

    print(res.IsUnstructured())

    res.ComputeGlobalOffset()
    res.AddElementToTag(1,"SecondElement")

    if res.GetNumberOfElements() != 2:
        raise Exception()
    res.ComputeGlobalOffset()

    print(res.GetDimensionality())
    res.ComputeBoundingBox()
    print(res.boundingMin)
    print(res.boundingMax)

    res.nodesTags.CreateTag('toto')
    print(res.GetNodalTag('toto'))
    print(res.GetNodalTag('toto2'))

    res.AddElementToTagUsingOriginalId(1,"bars")

    if res.GetPosOfNodes()[1,1] != 2:
        raise Exception()

    print(res.PrepareForOutput())
    print(res.GetElementsInTag("bars"))
    print(res.GetElementsInTag("bars",useOriginalId=True))

    print(res.GetElementsOfType(ElementNames.Bar_2).GetTag("toto"))
    print(res.GetPosOfNode(0))
    print(res)
    res.DeleteElemTags(["SecondElement"])
    print(res)

    resII = CreateMeshOfTriangles([[0,0,0],[1,2,3],[0,1,0]], [[0,1,2]])
    resII.AddNodeToTagUsingOriginalId(0,"First Point")
    resII.GenerateManufacturedOriginalIDs()

    if resII.GetPointsDimensionality() != 3:
        raise  #pragma: no cover

    if resII.GetElementsDimensionality() != 2:
        raise  #pragma: no cover

    resII.Clean()

    try:
        resII.MergeElements(res)
        raise #pragma: no cover
    except:
        pass

    resII.MergeElements(res,force=True)
    print("----")
    print(resII)
    print(resII.GetNumberOfElements(dim=2))
    print(resII.elements.GetTagsNames())
    del resII.elements[ElementNames.Triangle_3]
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
