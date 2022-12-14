# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Optional, List, Callable, Any, Union, Collection, Iterator, Tuple
from functools import reduce

import numpy as np

from BasicTools.NumpyDefs import ArrayLike, PBasicIndexType
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject  as BOO
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers.UnstructuredMesh import ElementsContainer, UnstructuredMesh

class Filter(BOO):
    """Base class to construct node and element filters

    Parameters
    ----------
    mesh : Optional[UnstructuredMesh], optional
        The mesh to be used, by default None
    zones : Optional[List[Callable]], optional
        a list of Callable (ImplicitGeometry for example), by default None
    tags : Optional[List[str]], optional
        The tags used to filter the elements, by default None
    zone : Optional[Callable], optional
        A Callable, by default None
    tag : Optional[str], optional
        A tag, by default None
    mask : Optional[ArrayLike], optional
        a boolean vector of the size of the object to filter, by default None
    """
    def __init__(self, mesh:Optional[UnstructuredMesh] = None, zones:Optional[List[Callable]] = None, tags:Optional[List[str]] = None, zone:Optional[Callable] = None, tag:Optional[str] = None, mask:Optional[ArrayLike] = None):

        super(Filter,self).__init__()

        self.tags = list()
        if tags is not None:
            self.SetTags(tags)
        if tag is not None:
            self.AddTag(tag)

        self.zones = list()
        if zones is not None:
            self.SetZones(zones)
        if zone is not None:
            self.AddZone(zone)

        self.mask = None
        if mask is not None:
            self.mask = mask

        self.mesh = mesh

    def IsEquivalent(self, other:Any)->bool:
        """To check if 2 element filter are equivalent ()

        Parameters
        ----------
        other : Any
            other object to check the equivalency

        Returns
        -------
        bool
            True if the two filters are equal
        """
        if id(self) == id(other):
            return True
        if isinstance(other, self.__class__):
            if sorted(self.tags) != sorted(other.tags):
                return False
            if sorted(self.zones) != sorted(other.zones):
                return False
            if not np.array_equal(self.mask, other.mask):
                return False
            return True
        else:
            return False

    def SetMesh(self, mesh:UnstructuredMesh):
        """Set the mesh

        Parameters
        ----------
        mesh : UnstructuredMesh
            The mesh to be used
        """
        self.mesh = mesh

    def SetTags(self, tagNames:List[str]):
        """Set the tag list name to treat

        Parameters
        ----------
        tagNames : List[str]
            the list of string with the tag names
        """
        self.tags = list(tagNames)

    def AddTag(self, tagName:str):
        """Add a tagname to the list of tag to treat

        Parameters
        ----------
        tagName : str
            the tag name to be added
        """
        self.tags.append(tagName)

    def SetZones(self, zonesList:List[Callable]):
        """Set the zone list to treat

        Parameters
        ----------
        zonesList : List[Callable]
            The list of zone to be
        """
        self.zones = zonesList

    def AddZone(self, zone:Callable):
        """dd a zone to the list of zone to be treated by the filter

        Parameters
        ----------
        zone : Callable
            a callable object capable of taking one argument with the
            points positions, and returning a vector of size pos.shape[0]
            with negative values for the entities to be selected by the filter
        """
        self.zones.append(zone)

    def _CheckTags_(self, tags, numberOfObjects) -> Union[np.ndarray,None]:
        """Internal function to compute the ids to be treated based on the tags

        Parameters
        ----------
        tags : _type_
            The tags container
        numberOfObjects : _type_
            the total number of object (number of points or number of element in the current element container)

        Returns
        -------
        np.ndarray or None
            np.ndarray of the the ids to be kept
            None if no tags present on to filter
        """

        if len(self.tags) == 0:
            return None

        tagList = list(tags[t] for t in self.tags if t in tags)

        if len(tagList) == 0:
            return []

        #fast path
        if len(tagList) == 1 :
            return tagList[0].GetIds()

        return reduce(np.union1d, (tag.GetIds()  for tag in tagList)  )

    def _CheckMask_(self, start:PBasicIndexType, size:PBasicIndexType) -> Union[np.ndarray,None]:
        """Internal function to compute the ids based on the mask

        Parameters
        ----------
        start : PBasicIndexType
            the start position for the current request
        size : PBasicIndexType
            the stop position for the current request

        Returns
        -------
        np.ndarray or None
            np.ndarray of the the ids to be kept
            None if no mask present on to filter
        """
        if self.mask is not None:
            return np.where(self.mask[start:start+size])[0]
        return None

    def intersect1D(self, first:Union[ArrayLike,None], second:Union[ArrayLike,None]) -> Union[np.ndarray,None]:
        """Function to generate an intersection of two vectors (like np.intersect1d)
        but with the particularity of treat the case where the inputs can be None
        None represent a non filtered inputs,

        Parameters
        ----------
        first : ArrayLike|None
            first vector of indices
        second : ArrayLike|None
            second vector of indices

        Returns
        -------
        Union[np.ndarray,None]
            The intersection of two list
        """
        if first is None:
            if second is None:
                return None
            else:
                return second
        else:
            if second is None:
                return first
            else:
                return np.intersect1d(first,second,assume_unique=True)

class NodeFilter(Filter):
    """class for node filtering zone, tag, mask
    the rest of the parameter are passed to the constructor of the base class Filter

    Parameters
    ----------
    mesh : Optional[UnstructuredMesh], optional
        The mesh to be used, by default None
    etags : Optional[List[str]], optional
       The element tags to be used to filter the nodes, by default None
    etag : Optional[str], optional
        A element tag name to be used to filter the nodes, by default None
    """
    def __init__(self, mesh:Optional[UnstructuredMesh]=None, etags:Optional[List[str]]=None, etag:Optional[str]=None, **kwargs):
        super().__init__(mesh=mesh,**kwargs)
        self.etags = list()
        if etags is not None:
            self.SetETags(etags)
        if etag is not None:
            self.AddETag(etag)

    def SetETags(self, tagNames:List[str]):
        """Set the tag list name to treat

        Parameters
        ----------
        tagNames : List[str]
            list of tag names
        """
        self.etags = list(tagNames)

    def AddETag(self, tagName:str):
        """Add a tagname to the list of tag to treat

        Parameters
        ----------
        tagName : str
            tag name to add
        """
        self.etags.append(tagName)

    def _CheckZones_(self, pos:ArrayLike, numberOfObjects:PBasicIndexType) -> Union[np.ndarray,None]:
        """Internal function to compute the ids to be treated based on the zones

        Parameters
        ----------
        pos : ArrayLike
            (n,3) size array with the positions to be treated
        numberOfObjects : PBasicIndexType
            total number of points

        Returns
        -------
        Union[np.ndarray,None]
            list of nodes to treat or None for all nodes
        """
        if len(self.zones) == 0:
            return None

        if len(self.zones) == 1 :
            return np.where(self.zones[0](pos)<=0)[0]

        return reduce(np.logical_or, (zone(pos)<=0 for zone in self.zones )  )

    def GetIdsToTreat(self,notUsed:Any=None)-> Union[np.ndarray,Collection]:
        """Get the nodes selected by this filter

        Parameters
        ----------
        notUsed : Any, optional
            Not Used, by default None

        Returns
        -------
        Union[np.ndarray,Collection]
            The filtered entities
        """
        if len(self.etags) > 0:
            ff =ElementFilter(self.mesh,tags=self.etags)
            class OP():
                def __init__(self):
                    self.set = set()

                def __call__(self,name,data,ids):
                    self.set.update(data.GetNodesIdFor(ids))
            op = OP()
            ff.ApplyOnElements(op)
            if len(op.set) > 0:
                resE = list(op.set)
            else:
                resE = None
        else:
            resE = None

        res  = self._CheckTags_(self.mesh.nodesTags,self.mesh.GetNumberOfNodes())
        res2 = self._CheckZones_(self.mesh.nodes,self.mesh.GetNumberOfNodes() )
        resM = self._CheckMask_(0,self.mesh.GetNumberOfNodes())
        res3 = self.intersect1D(res,res2)
        res3 = self.intersect1D(res3,resE)
        res3 = self.intersect1D(res3,resM)
        if res3 is None:
            return range(self.mesh.GetNumberOfNodes())
        else:
            return res3

    def ApplyOnNodes(self,op:Callable):
        """Apply the filter using an operator

        Parameters
        ----------
        op : Callable
            An instance of a callable object, the object can have
            the PreCondition function and/or the PostCondition function. Theses
            functions are called (if exist) (with the mesh as the first argument)
            before and after the main call ( op(mesh,nodes,ids) )
        """
        pc = getattr(op,"PreCondition",None)

        if callable(pc):
            pc(self.mesh)

        op(self.mesh,self.mesh.nodes,self.GetIdsToTreat() )

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)

class FilterOP(BOO):
    """Base class to construct derived Filters (Union filters for example)

    Parameters
    ----------
    mesh : Optional[UnstructuredMesh], optional
        The mesh to work on, by default None
    filters : List[Union[Filter,FilterOP]], optional
        the list of filters to use, by default None

    """
    def __init__(self,mesh:Optional[UnstructuredMesh]=None,filters:List[Union[Filter,FilterOP]]=None):

        super(FilterOP,self).__init__()

        self._mesh= None
        if filters is not None:
            self.filters = filters
        else:
            self.filters = []

        if mesh is not None:
            self.mesh = mesh
        else:
            #if no mesh provided we pick the mesh of the first object
            #(if available)
            if len(self.filters):
                self.mesh = self.filters[0].mesh
        self.withError = False

    def IsEquivalent(self, other:Any) -> bool:
        """To check if 2 element filter are equivalent ()

        Parameters
        ----------
        other : Any
            other object to check the equivalency

        Returns
        -------
        bool
            True if the two filters are equal
        """

        if id(self) == id(other):
            return True
        if isinstance(other, self.__class__):
            if self.filters != other.filters:
                return False
            return True
        else:
            return False

    def SetMesh(self, mesh:UnstructuredMesh):
        """Set the mesh, this will push the mesh to the internal filters

        Parameters
        ----------
        mesh : UnstructuredMesh
            The mesh to be used
        """
        self._mesh = mesh

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, m):
        self._mesh = m
        for f in self.filters:
            f.mesh = m

    @property
    def zoneTreatment(self, zt):
        raise Exception("Cant ask zoneTreatment to a FilterOP")

    @zoneTreatment.setter
    def zoneTreatment(self, zt):
        for f in self.filters:
            f.zoneTreatment = zt

    def Complementary(self):
        # the complementary of the complementary is the original filter
        if isinstance(self,ComplementaryObject):
            return self.filters[0]
        return ComplementaryObject(mesh=self.mesh,filters=[self])

    def GetFrozenFilter(self, mesh=None):
        if isinstance(self,FrozenFilter):
            if mesh is None:
                return self
            if mesh is self.mesh:
                return self
            else:
                raise Exception("Can't freeze a FrozenFilter with a new mesh")

        if mesh is not None:
            return FrozenFilter(mesh=mesh, filters=[self])

        if self.mesh is None:
            raise Exception("Need to set the mesh first on the filter or provide one")
        else:
            return FrozenFilter(mesh=self.mesh, filters=[self])

    def __iter__(self):
        """
        Iteration interface to ease the use of the filter

        :example:

            myFilter = DerivedElementFilterName(myMesh)
            myFilter.filters.append(myOtherFilter1)
            myFilter.filters.append(myOtherFilter2)

            for name,elements,ids in myFilter:

                print("This function is called on the union of the 2 filters")

                print("Number of element of type " + str(name)+ " is : "  + str(len(ids) )
        """
        elementsFound = False
        for name,data in self.mesh.elements.items():
            ids = self.GetIdsToTreat(data)
            if len(ids) == 0: continue
            elementsFound = True
            yield name, data, ids

        if elementsFound is False and self.withError:
            raise Exception("Error!! Zero element in the element filter : \n" + str(self))

    def ApplyOnElements(self,op):
        """
        Function to apply the filter  using an operator

        :param callable op: An instance of a callable object, the object can have
            the PreCondition function and/or the post-condition function. Theses
            functions are called (if exist) (with the mesh as the first argument)
            before and after the main call ( op(name,elements,ids) )
        """
        op.PreCondition(self.mesh)
        for name,elements,ids in self:
            op(name,elements,ids)
        op.PostCondition(self.mesh)
        return op

    def ApplyOnNodes(self,op):
        """
        Function to apply filter using an operator

        :param callable op: An instance of a callable object, the object can have
            the PreCondition function and/or the post-condition function. Theses
            functions are called (if exist) (with the mesh as the first argument)
            before and after the main call ( op(mesh,nodes,ids) )
        """

        pc = getattr(op,"PreCondition",None)

        if callable(pc):
            pc(self.mesh)

        op(self.mesh,self.mesh.nodes,self.GetIdsToTreat(None) )

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)

class UnionElementFilter(FilterOP):
    """
      Specialized class to compute the union of filter (add)
    """
    def __init__(self,mesh=None,filters=None):
        super(UnionElementFilter,self).__init__(mesh=mesh,filters=filters)

    def GetIdsToTreat(self, data):
        return reduce(np.union1d, (np.asarray(ff.GetIdsToTreat(data),dtype=PBasicIndexType) for ff in self.filters))

class IntersectionElementFilter(FilterOP):
    """
      Specialized class to compute the intersection of filters
    """
    def __init__(self,mesh=None,filters=None):
        super(IntersectionElementFilter,self).__init__(mesh=mesh,filters=filters)

    def GetIdsToTreat(self, data):
        ids = None
        for ff in self.filters:
            if ids is None:
                ids = ff.GetIdsToTreat(data)
            else:
                ids = np.intersect1d(ids,ff.GetIdsToTreat(data),assume_unique=True )
            if len(ids) == 0:
                return []
        if ids is None :
            return []
        return list(ids)

class DifferenceElementFilter(FilterOP):
    """
      Specialized class to compute the difference between two filters
    """
    def __init__(self,mesh=None,filters=None):
        if filters is not None and len(filters) != 2:
            raise Exception("Need exactly 2 filter to compute the difference")
        super(DifferenceElementFilter,self).__init__(mesh=mesh,filters=filters)

    def GetIdsToTreat(self, data):
        ids = None

        idsA = self.filters[0].GetIdsToTreat(data)
        idsB = self.filters[1].GetIdsToTreat(data)

        return np.setdiff1d(idsA, idsB)

class ComplementaryObject(FilterOP):
    def __init__(self,mesh=None,filters=None):
        super(ComplementaryObject,self).__init__(mesh=mesh,filters=filters)
        if len(self.filters) > 1 :
            raise Exception("ComplementaryObject Error!: filters must be of len = 1")

    def GetIdsToTreat(self, data):
        if len(self.filters)  > 1:
            raise Exception("ComplementaryObject Error!: filters must be of len = 1")

        f = self.filters[0]
        ids = f.GetIdsToTreat(data)
        if len(ids) == data.GetNumberOfElements(): return []
        mask = np.ones(data.GetNumberOfElements(),dtype=bool)
        mask[ids] = False
        return np.where(mask)[0]


class PartialElementFilter(FilterOP):
    """ Utility class to create a partition of a ElementFilter"""

    def __init__(self, elementFilter, partitions, partitionNumber):
        super().__init__(self,[elementFilter])
        self.partitions = partitions
        self.partitionNumber = partitionNumber

    def GetIdsToTreat(self,elements):
        res = self.filters[0].GetIdsToTreat(elements)
        return np.array_split(res,self.partitions)[self.partitionNumber]

    def Complementary(self):
        for name,data,ids  in self.filters[0].Complementary():
            ids = np.array_split(ids,self.partitions)[self.partitionNumber]
            if len(ids) == 0:
                continue
            yield name, data, ids

class IdsAsNumpyMask(FilterOP):
    def __init__(self,mesh=None,filters=None):
        super(IdsAsNumpyMask,self).__init__(mesh=mesh,filters=filters)
        if len(self.filters) > 1 :
            raise Exception("IdsToMask Error!: filters must be of len = 1")

    def __iter__(self):
        for name, data, ids  in self.filters[0]:
            mask = np.ones(data.GetNumberOfElements(),dtype=bool)
            mask[ids] = False
            yield name, data, mask

class ElementFilter(Filter):
    """Class for element filtering by dimensionality, zone, mask,elementType, and tag
    for the zones three types of treatments are possible:
      if the the center of the element is inside the zone   : self.zoneTreatment = "center"
      if all nodes of the element are inside the zone       : self.zoneTreatment = "allnodes"
      if at least one node of the element is inside the zone: self.zoneTreatment = "leastonenode"

    Parameters
    ----------
    mesh : UnstructuredMesh, optional
        The mesh to be used
    dimensionality : int, optional
        the dimensionality of the element to be included in this filter, by default None
        possible option are: -3 -2 -1 0 1 2 3 None
        the - sign is for the complementary part (-2 = all non 2D elements)
    elementType : Optional[str], optional
        the name of a element type to be included in this filter, by default None
    elementTypes : Optional[List[str]], optional
        a list of element type to be included in this filter, by default None
    zoneTreatment : str, optional
        ["center" | "allnodes" | "leastonenode"], by default "center"
    nTags : Optional[List[str]], optional
        a list of nodal tags names to use to extract elements
    nTagsTreatment : str, optional
        ["allnodes" | "leastonenode"], by default "allnodes"
    """
    def __init__(self, mesh:UnstructuredMesh=None,
                 dimensionality:Optional[int]=None,
                 elementTypes:Optional[List[str]]=None,
                 elementType:Optional[str]=None,
                 zoneTreatment:str="center",
                 nTags:Optional[List[str]] = None,
                 nTagsTreatment:str="allnodes", **kwargs):

        super(ElementFilter,self).__init__(mesh=mesh,**kwargs)
        self.dimensionality = dimensionality
        self.zoneTreatment = zoneTreatment # "center", "allnodes", "leastonenode"
        self.nTagsTreatment = nTagsTreatment

        self.elementTypes = list()

        if elementTypes is not None:
            self.SetElementTypes(elementTypes)

        if elementType is not None:
            self.AddElementType(elementType)

        self.nTags = list()
        if nTags is not None:
            self.SetNTags(nTags)

        self.withError = False

    def IsEquivalent(self, other:Any) -> bool:
        """To check if 2 element filter are equivalent ()

        Parameters
        ----------
        other : Any
            other object to check the equivalency

        Returns
        -------
        bool
            True if the two filters are equal
        """
        res = super(ElementFilter,self).IsEquivalent(other)
        if not res :
            return False

        if id(self) == id(other):
            return True
        if isinstance(other, self.__class__):
            if self.dimensionality != other.dimensionality:
                return False
            if self.zoneTreatment != other.zoneTreatment:
                return False
            if self.elementTypes != other.elementTypes:
                return False
            return True
        else:
            return False

    def __str__(self) -> str:
        res = "ElementFilter\n"
        res += "  dimensionality: "+ str(self.dimensionality) + " \n"
        res += "  tags          : "+ str(self.tags) + " \n"
        res += "  zones         : "+ str(self.zones) + " \n"
        res += "  zoneTreatment : "+ str(self.zoneTreatment) + " \n"
        res += "  elementTypes  : "+ str(self.elementTypes) + " \n"
        return res

    def SetZoneTreatment(self, zoneTreatment:str):
        """Set the way the elements are selected based on the position
        if the the center of the element is inside the zone   : self.zoneTreatment = "center"
        if all nodes of the element are inside the zone       : self.zoneTreatment = "allnodes"
        if at least one node of the element is inside the zone: self.zoneTreatment = "leastonenode"

        Parameters
        ----------
        zoneTreatment : str
            ["center" | "allnodes" | "leastonenode"]

        Raises
        ------
        Exception
            if the string is not permitted
        """
        if zoneTreatment in ["center", "allnodes", "leastonenode"]:
            self.zoneTreatment = zoneTreatment
        else:
            raise Exception(f"Zone treatment not valid ({zoneTreatment}), possible options are : center, allnodes, leastonenode")

    def SetDimensionality(self,dimensionality:int):
        """ Set the dimensionality of the elements to be treated

        Parameters
        ----------
        dim : int
            the dimensionality filter, [-3 -2 -1 0 1 2 3 or None]
            the - sign is for the complementary part (-2 = all non 2D elements).
            Set to None to reset (not to apply dimensionality as a criteria)
        """
        self.dimensionality = dimensionality

    def SetElementTypes(self, elementTypes:List[str]):
        """Set the names of the element types to be included in this filter

        Parameters
        ----------
        elementTypes : List[str]
            the list of element types
        """
        self.elementTypes = []
        self.elementTypes.extend(elementTypes)

    def SetNTags(self, nTags:List[str]):
        """Set the names of the nodal tags (nTags) to be included in this filter

        Parameters
        ----------
        elementTypes : List[str]
            the list of element types
        """
        self.nTags = []
        self.nTags.extend(nTags)

    def SetNTagsTreatment(self, nTagsTreatment:str):
        """Set the way the elements are selected based on the nTags
        if all nodes of the element are inside the nTag   : self.nTagsTreatment = "allnodes"
        if at least one node of the element is inside nTag: self.nTagsTreatment = "leastonenode"

        Parameters
        ----------
        nTagsTreatment : str
            ["allnodes" | "leastonenode"]

        Raises
        ------
        Exception
            if the nTagsTreatment string is not in ["allnodes", "leastonenode"]
        """
        if nTagsTreatment in ["allnodes", "leastonenode"]:
            self.nTagsTreatment = nTagsTreatment
        else:
            raise Exception(f"NTag treatment not valid ({nTagsTreatment}), possible options are : allnodes, leastonenode")

    def AddElementType(self, elementType:str):
        """Add an element type to be included in this filter

        Parameters
        ----------
        elementType : str
            the name of an element type
        """
        self.elementTypes.append(elementType)

    def _CheckDimensionality_(self, elements:ElementsContainer) -> Union[bool,None]:
        """Internal function check if a type of element must be included based on
        the dimensionality criteria

        Parameters
        ----------
        elements : ElementsContainer
            the incoming ElementsContainer

        Returns
        -------
        Union[bool,None]
            True if this type of elements must be included
            False if this type of elements must be excluded
            None if this the filtering by dimensionality is not active
        """
        if self.dimensionality is None:
            return None
        else:
            elementDimension = EN.dimension[elements.elementType]
            if self.dimensionality  >= 0:
                if elementDimension != self.dimensionality:
                    return False
            else:
                if elementDimension == (-self.dimensionality):
                    return False
        return True

    def _CheckElementTypes_(self, elements:ElementsContainer) -> Union[bool,None]:
        """Internal function check if a type of element must be included based on
        the elementType criteria

        Parameters
        ----------
        elements : ElementsContainer
            the incoming ElementsContainer

        Returns
        -------
        Union[bool,None]
            True if this type of elements must be included
            False if this type of elements must be excluded
            None if this the filtering by dimensionality is not active
        """

        if len(self.elementTypes) == 0:
            return None
        if elements.elementType in self.elementTypes:
            return True
        else:
            return False

    def _CheckZones_(self, elements:ElementsContainer) -> Union[np.ndarray,None]:
        """Internal function check if a type of element must be included based on
        the zone

        Parameters
        ----------
        elements : ElementsContainer
            the incoming ElementsContainer

        Returns
        -------
        Union[np.array,None]
            np.ndarray of the indices inside the selection
            None if this the filtering by zones is not active
        """

        if len(self.zones) == 0:
            return None

        if self.zoneTreatment == "center":
            from BasicTools.Containers.MeshTools import GetElementsCenters
            centers = GetElementsCenters(nodes=self.mesh.nodes,elements=elements)

        numberOfObjects = elements.GetNumberOfElements()

        res = np.zeros(numberOfObjects,dtype=bool)

        for zone in self.zones :
            if self.zoneTreatment == "center":
                res2 = zone(centers)<=0
            elif self.zoneTreatment == "allnodes":
                z = zone(self.mesh.nodes)<=0
                res2 = np.sum(z[elements.connectivity],axis=1) == elements.GetNumberOfNodesPerElement()
            elif self.zoneTreatment == "leastonenode":
                z = zone(self.mesh.nodes)<=0
                res2 = np.sum(z[elements.connectivity],axis=1) > 0
            else:#pragma: no cover
                raise Exception("zoneTreatment unknown")

            np.logical_or(res, res2 ,out=res)

        return np.where(res)[0]

    def _CheckNTags_(self, elements:ElementsContainer) -> Union[np.ndarray,None]:
        """Internal function check if a element must be included based on
        the nTags

        Parameters
        ----------
        elements : ElementsContainer
            the incoming ElementsContainer

        Returns
        -------
        Union[bool,None]
            np.ndarray of the indices inside the selection
            None if this the filtering by nTags is not active
        """

        if len(self.nTags) == 0:
            return None

        numberOfObjects = elements.GetNumberOfElements()

        res = np.zeros(numberOfObjects,dtype=bool)

        nodalMask = np.zeros(self.mesh.GetNumberOfNodes(),dtype = bool)

        for tagName in self.nTags :
            if tagName in self.mesh.nodesTags:
                nodalMask[self.mesh.nodesTags[tagName].GetIds()] = True

        if self.nTagsTreatment == "allnodes":
            res = np.sum(nodalMask[elements.connectivity],axis=1) == elements.GetNumberOfNodesPerElement()
        elif self.nTagsTreatment == "leastonenode":
            res = np.sum(nodalMask[elements.connectivity],axis=1) > 0
        else: #pragma: no cover
            raise Exception("nTagsTreatment unknown")
        return np.where(res)[0]

    def GetIdsToTreat(self,elements:ElementsContainer) -> Union[np.ndarray,Collection]:
        """Get the entities selected by this filter

        Parameters
        ----------
        elements : ElementsContainer
            Elements to treat

        Returns
        -------
        Union[np.ndarray,Collection]
            The filtered entities
        """
        if self._CheckDimensionality_(elements) == False:
            return []

        if self._CheckElementTypes_(elements) == False:
            return []

        res  = self._CheckTags_(elements.tags,elements.GetNumberOfElements())
        if res is not None and len(res) == 0:
            return res

        res_nTags  = self._CheckNTags_(elements)
        if res_nTags is not None and len(res_nTags) == 0:
            return res_nTags
        res = self.intersect1D(res, res_nTags)

        res_zones = self._CheckZones_(elements)
        # this is a hack to be able to access the evaluated zone over the used points
        # for the moment this is not part of the public API
        self.zonesField = res_zones
        res3 = self.intersect1D(res, res_zones)

        init = 0
        for name, data in self.mesh.elements.items():
            if name == elements.elementType:
                break
            init +=  data.GetNumberOfElements()

        res_mask = self._CheckMask_(init,elements.GetNumberOfElements())
        res3 = self.intersect1D(res3, res_mask)

        if res3 is None:
            return range(elements.GetNumberOfElements())
        else:
            return res3

    def Complementary(self) -> ComplementaryObject:
        """Create a filter with the complementary part

        Returns
        -------
        ComplementaryObject
            the new filter
        """

        # the complementary of the complementary is the original filter
        if isinstance(self,ComplementaryObject):
            return self.filters[0]
        return ComplementaryObject(mesh=self.mesh,filters=[self])

    def GetFrozenFilter(self, mesh:Optional[UnstructuredMesh]=None) -> FrozenFilter:
        """Generate a frozen filter. This is a filter with pre-evaluated ids.
        This class is useful when a repeated use of a filter is needed

        Parameters
        ----------
        mesh : Union[UnstructuredMesh], optional
            the mesh to pre evaluated the filter.
            If None we use the internal mesh stored in the filter, by default None

        Returns
        -------
        FrozenFilter
            A FrozenFilter pre-evaluated on the mesh

        Raises
        ------
        Exception
            If the user tries to freeze a frozen filter on a different mesh
        Exception
            if no mesh is available
        """

        if isinstance(self,FrozenFilter):
            if mesh is None:
                return self
            else:
                raise Exception("Can't freeze a FrozenFilter with a new mesh")

        if mesh is not None:
            return  FrozenFilter(mesh=mesh, filters=[self])

        if self.mesh is None:
            raise Exception("Need to set the mesh first on the filter or provide one")
        else:
            return FrozenFilter(mesh=self.mesh, filters=[self])

    def __iter__(self) -> Iterator[Tuple[str,ElementsContainer,ArrayLike]]:
        """Iterator over the selection

        Yields
        ------
        Iterator[Tuple[str,ElementsContainer,ArrayLike]]
            a tuple with three elements:
              - the elementType
              - the ElementContainer
              - a np.ndarray with the ids to treat

        Raises
        ------
        Exception
            if withError is activated raise an exception to help the debugging


        :example:

            myFilter = ElementFilter(myMesh,dimensionality=2)

            for name,elements,ids in myFilter:

                print("This function is called only for 2D elements")

                print("Number of 2D element of type " + str(name)+ " is : "  + str(len(ids) )
        """


        elementsFound = False
        for name,data in self.mesh.elements.items():
            ids = self.GetIdsToTreat(data)
            if len(ids) == 0: continue
            elementsFound = True
            yield name, data, ids

        if elementsFound  == False and self.withError:
            raise Exception("Error!! Zero element in the element filter : \n" + str(self))

    def ApplyOnElements(self, op:Callable):
        """Apply the filter using an operator

        Parameters
        ----------
        op : Callable
            An instance of a callable object, the object can have
            the PreCondition function and/or the PostCondition function. Theses
            functions are called (if exist) (with the mesh as the first argument)
            before and after the main call ( op(mesh,nodes,ids) )

        Returns
        -------
        Callable
            return op after the application of the filter
        """
        pc = getattr(op,"PreCondition",None)

        if callable(pc):
            pc(self.mesh)

        for name,elements,ids in self:
            op(name,elements,ids)

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)

        return op

class FrozenFilter(FilterOP):
    """Class to hold a pre-evaluated filter

    Parameters
    ----------
    mesh : Optional[UnstructuredMesh], optional
        The mesh to work on, by default None
    filters : Union[Filter,FilterOP], optional
        the list with only one filter to use, by default None

    """
    def __init__(self,mesh=None,filters=None):
        self.__frozenData = {}
        super(FrozenFilter,self).__init__(mesh=mesh,filters=filters)
        if len(self.filters) > 1 :
            raise Exception("ComplementaryObject Error!: filters must be of len = 1")

    @property
    def mesh(self):
        return super(FrozenFilter,self).mesh

    @mesh.setter
    def mesh(self, m):
        if m is None:
            return
        if self._mesh is not None :
            if m is self._mesh:
                return
            raise Exception("You cant set the mesh 2 times")
        self._mesh = m
        if len(self.filters) > 0:
            self.filters[0].mesh = m
            for name,data in self.filters[0].mesh.elements.items():
                self.__frozenData[name] = (data, self.filters[0].GetIdsToTreat(data))

    def IsEquivalent(self, other):
        for name,(data,ids) in self.__frozenData.items():
            ids2 = other.GetIdsToTreat(data)
            if len(ids) != len(ids2):
                return False
            if np.any(np.asarray(ids) != np.asarray(ids2)):
                return False
        return True

        return self.filters[0].IsEquivalent(other)

    def SetIdsToTreatFor(self, elements, localIds):
        self.__frozenData[elements.elementType] = (elements, localIds)

    def GetIdsToTreat(self, elements):
        return self.__frozenData[elements.elementType][1] # (data, ids)

    def __iter__(self)-> Iterator[Tuple[str,ElementsContainer,ArrayLike]]:
        """Iterator over the selection

        Yields
        ------
        Iterator[Tuple[str,ElementsContainer,ArrayLike]]
            a tuple with three elements:
              - the elementType
              - the ElementContainer
              - a np.ndarray with the ids to treat
        """
        for name,(data,ids) in self.__frozenData.items():
            if len(ids) == 0:
                continue
            yield name, data, ids

class ElementFilterBaseOperator():
    """Template for the creation of operators
    """
    def PreCondition(self,mesh):
        pass

    def __call__(self,name,data,ids):
        pass

    def PostCondition(self,mesh):
        pass

class ElementCounter(ElementFilterBaseOperator):
    """Basic ElementCounter operators

    Usage:
        numberOfElement = ElementFilter(myMesh,dimensionality=2).ApplyOnElements(ElementCounter()).cpt

    """
    def __init__(self):
        self.cpt = 0

    def PreCondition(self, mesh:UnstructuredMesh):
        """function called at the beginning to reset the counter
        """
        self.cpt = 0

    def __call__(self,name,data,ids):
        """We add the current number of element to the global counter
        """
        self.cpt += len(ids)

def ElementFilterToImplicitField(elementFilter:ElementFilter, pseudoDistance:int=2) -> np.ndarray:
    """Function to generate an iso zero level-set on the mesh to represent
    the shape of the filter. This discretized iso zero on the mesh cant always
    hold a 'perfect' representation of the filter, so a modified iso zero is
    created. An additional parameter pseudo-distance can be increased to create
    a pseudo distance. This field is created using the connectivity and not
    the real distances.

    Parameters
    ----------
    elementFilter : ElementFilter
        the element filter to process
    pseudoDistance : int, optional
        the number of element to propagate the pseudo-distance, by default 2

    Returns
    -------
    np.ndarray
        a field over the nodes with negative values inside the domain defined by the filter
    """
    def UpdateInsideOutsideNodes(mesh, phi, insideNodes=None,outsideNodes=None, iso=0.0):
        """ function to build masks (insideNodes, outsideNodes) with the information
        about if a particular nodes is connected (through an element) to the
        inside phi < iso or to the outside phi > iso
        """
        for name, data  in mesh.elements.items():
            if mesh.GetDimensionality() == EN.dimension[name]:
                phis = phi[data.connectivity]
                if insideNodes is not None:
                    elementMask = np.any(phis<iso,axis=1)
                    insideNodes[ data.connectivity[elementMask] ] = True
                if outsideNodes is not None:
                    elementMask = np.any(phis>iso,axis=1)
                    outsideNodes[ data.connectivity[elementMask] ] = True

    mesh = elementFilter.mesh
    phi = np.zeros(mesh.GetNumberOfNodes())
    insideNodes = np.zeros(mesh.GetNumberOfNodes(),dtype=bool)
    dim = 0
    for name, data, ids  in elementFilter:
        insideNodes[data.connectivity[ids,:]]  = True
        dim = max(dim,EN.dimension[name])

    outsideNodes = np.zeros(mesh.GetNumberOfNodes(),dtype=bool)
    for name, data, ids  in elementFilter.Complementary():
        outsideNodes[data.connectivity[ids,:]]  = True

    phi[insideNodes] = -1
    phi[outsideNodes] = 1
    phi[np.logical_and(insideNodes,outsideNodes)] = 0

    And = np.logical_and

    if dim == mesh.GetDimensionality():
        # correction of point values
        # if a point have only zeros or negative values on neighbors point then the point is set to inside
        # if a point have only zeros or positive values on neighbors point then the point is set to outside

        insideNodes.fill(False)
        outsideNodes.fill(False)

        UpdateInsideOutsideNodes(mesh,phi,insideNodes,outsideNodes, 0)
        mask = phi == 0
        phi[And(mask, And(np.logical_not(outsideNodes),insideNodes ))] = -1/2
        phi[And(mask, And(np.logical_not(insideNodes), outsideNodes))] = 1/2


    insideWork = True
    outsideWork = True
    for i in range(1, pseudoDistance):
        if insideWork:
            insideNodes.fill(False)
            UpdateInsideOutsideNodes(mesh, phi, insideNodes, None, float(i) )
            mask = phi == i
            finalMask = And(mask, np.logical_not(insideNodes) )
            if np.any(finalMask):
                phi[finalMask] = i+1
            else:
                insideWork = False

        if outsideWork:
            outsideNodes.fill(False)
            UpdateInsideOutsideNodes(mesh, phi, None, outsideNodes, float(-i) )
            mask = phi == -i
            finalMask = And(mask, np.logical_not(outsideNodes) )
            if np.any(finalMask):
                phi[finalMask] = -(i+1)
            else:
                outsideWork = False

        if not outsideWork and not insideWork:
            break
    return phi

def CheckIntegrity( GUI=False):
    """
    .. literalinclude:: ../../src/BasicTools/Containers/Filters.py
       :pyobject: CheckIntegrity
    """
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    nNodesX = 11
    nNodesY = 12
    nNodesZ = 13
    mesh = CreateCube(dimensions=[nNodesX,nNodesY,nNodesZ],origin=[0,0,0.], spacing=[1./(nNodesX-1),1./(nNodesY-1), 10./(nNodesZ-1)], ofTetras=True )
    print(mesh)

    class NOP():
        def __init__(self):
            self.cpt = 0

        def PreCondition(self,mesh):
            self.cpt = 0

        def __call__(self,mesh,node,ids):
            self.cpt += len(ids)

        def PostCondition(self,mesh):
            print("The counter is at {}".format(self.cpt) )

    ff = NodeFilter(mesh)
    ff.ApplyOnNodes(NOP())


    ff = NodeFilter(mesh,tags=["x0y0z0","x0y0z1",],tag="x1y0z0")
    ff.AddTag("x1y0z1")


    op = NOP()
    print(ff)
    ff.ApplyOnNodes(op)

    if op.cpt != 4:#pragma: no cover
        raise Exception("Error finding the point")

    ff.AddZone(lambda p: np.ones(p.shape[0])-1)

    op = NOP()
    print(ff)
    ff.ApplyOnNodes(op)
    if op.cpt != 4:#pragma: no cover
        raise Exception("Error finding the point")

    # example of counting the number of element in the eTag X0
    cpt = 0
    ff = ElementFilter(mesh)
    ff.AddTag("X0")
    for name, data, ids in ff:
        cpt += len(ids)

    if cpt != 264: # pragma: no cover
        raise
    print("Number of Element in tag X0 {}".format(cpt))

    cpt = 0
    ff = ElementFilter(mesh,nTags=["x1y0z1"],nTagsTreatment="leastonenode")
    for name,data,ids in ff:
        cpt += len(ids)

    if cpt != 6: # pragma: no cover
        raise
    print(f"Number of Element touching the node tag 'x1y0z1' {cpt}")

    mask = np.zeros(mesh.GetNumberOfElements(),dtype=bool)
    mask[0] = True
    ff = ElementFilter(mesh,mask=mask)
    for n,d,ids in ff:
        print(n,ids)
        if n != EN.Tetrahedron_4 or len(ids) != 1 :
            raise  # pragma: no cover
        if ids[0] != 0:
            raise  # pragma: no cover

    ff = ElementFilter(mesh, zone = lambda p: (p[:,2]-mesh.boundingMin[2]-0.001),zones=[])

    ff.GetFrozenFilter().GetFrozenFilter(mesh).ApplyOnElements(ElementFilterBaseOperator())

    ff = ElementFilter(mesh, zone = lambda p: (p[:,2]-mesh.boundingMin[2]-0.001),zones=[])
    for name,data,ids in ff:
        data.tags.CreateTag("ZZ0").SetIds(ids)

    ## to check if a filter can be used 2 times
    for name,data,ids in ff:
        data.tags.CreateTag("ZZ0O").SetIds(ids)


    # example of counting the number of element in the eTag ZZ0

    class OP(ElementCounter):
        def __init__(self):
            super(OP,self).__init__()

        def __call__(self,name,data,ids):
            super(OP,self).__call__(name,data,ids)
            print(name)

        def PostCondition(self,mesh):
            print("The counter is at {}".format(self.cpt) )

    op = OP()

    ff.SetZoneTreatment("allnodes")
    ff.ApplyOnElements(op)
    ff.zoneTreatment = "leastonenode"
    ff.ApplyOnElements(op)



    cpt = 0
    ff = ElementFilter(mesh,elementTypes=[], elementType=EN.Triangle_3)
    ff.SetElementTypes([EN.Triangle_3])
    print(ff)
    def MustFail(func,*args,**kwargs):
        try:
            func(*args,**kwargs)
        except:
            pass
        else:
            raise Exception("Error in the CheckIntegrity ")# pragma no cover

    MustFail(ff.SetZoneTreatment,"toto")
    ff.ApplyOnElements(op)

    MustFail(ComplementaryObject,filters=[ff,ff])

    ff.ApplyOnElements(op)

    f2 = ElementFilter(mesh,elementType=EN.Triangle_3)
    f2.AddTag("ZZ0")

    ff.SetDimensionality(2)

    op = f2.ApplyOnElements(OP())

    if op.cpt != (2*(nNodesX-1)*(nNodesY-1)) : # pragma: no cover
        raise Exception("Error in the number of elements in the tag = " + str(op.cpt)+ " must be " + str((2*(nNodesX-1)*(nNodesY-1))))

    IF = IntersectionElementFilter(mesh=mesh,filters=[ff,f2])
    IF.ApplyOnElements(OP())

    UF = UnionElementFilter(mesh=mesh,)
    UF.filters=[ff,f2]
    UF.ApplyOnElements(OP())
    print(UF.Complementary())
    NUF = UnionElementFilter(mesh=mesh,filters=[NodeFilter(mesh,etag="3D"),NodeFilter(mesh)])
    NUF.ApplyOnNodes(NOP())

    print("Number of Element in tag ZZ0 {}".format(op.cpt))
    print(mesh)

    if GUI: # pragma: no cover
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

    ff = ElementFilter(mesh,tag="Some")
    mesh.elements[EN.Tetrahedron_4].tags.CreateTag("Some").SetIds(list(range(4000)))


    from BasicTools.Helpers.Timer import Timer
    a = Timer("ElementFilterToImplicitField").Start()
    phi = ElementFilterToImplicitField(ff, pseudoDistance=50)
    a.Stop()
    print(a)

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()

    WriteMeshToXdmf(tempdir+"test.xdmf",mesh, PointFields=[phi],PointFieldsNames=["Phi"] )

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=False)) # pragma: no cover
