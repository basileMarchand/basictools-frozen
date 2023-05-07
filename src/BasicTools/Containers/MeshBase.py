# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from __future__ import annotations
from typing import Dict

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.NumpyDefs import PBasicIndexType
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Tags import Tag, Tags

class MeshBase(BaseOutputObject):
    """This is the base mesh of all the meshe classes in BasicTools
    """
    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = Tags()
        """Tags for the nodes """
        self.nodeFields: Dict[str,np.ndarray] = {}
        """ This is a dictionary containing fields defined over the nodes (all the nodes).
        The keys are the names of the fields
        the values are the actual data of size (nb nodes, nb of components)"""
        self.elemFields = {}
        """ This is a dictionary containing fields defined over the elements (all the elements).
        The keys are the names of the fields
        the values are the actual data of size (nb elements, nb of components)"""
        self.props = {}
        """Metadata this is just a dictionary that can be used to transport
        information with the mesh, please use the class name as key of the
        object using/setting the information
        """

    def __copy__(self):
        res = MeshBase()
        res._assign(self)
        return res

    def _assign(self, other):
        self.nodesTags = other.nodesTags
        self.nodeFields = other.nodeFields
        self.elemFields = other.elemFields
        self.props = other.props

    def __eq__(self, other):
        if other is None:
            return False

        if self.nodesTags != other.nodesTags:
            return False

        if len(self.nodeFields) != len(other.nodeFields):
            return False

        for k,v in self.nodeFields.items():
            if k in other.nodeFields:
                if not np.array_equal(v, other.nodeFields[k]):
                    return False
            else:
                return False

        if len(self.elemFields) != len(other.elemFields):
            return False

        for k,v in self.elemFields.items():
            if k in other.elemFields:
                if not np.array_equal(v, other.elemFields[k]):
                    return False
            else:
                return False

        if self.props != other.props:
            return False

        return True

    def GetElementsOfType(self, typename:str):
        """return the element container for the element name (typename)

        Parameters
        ----------
        typename : str
            the name of the elements to extract

        Returns
        -------
        ElementContainer
            The Element Container for element
        """
        return self.elements.GetElementsOfType(typename)

    def GetNamesOfElemTags(self) -> List[str]:
        """return a list containing all the element tags present in the mesh

        Returns
        -------
        List[str]
            the list of all the name of the element tags present in the mesh
        """
        res = set()
        for ntype, data in self.elements.items():
            for tag in data.tags:
                res.add(tag.name)

        return list(res)

    def CopyProperties(self, other):
        import copy
        self.props = copy.deepcopy(self.props)

    def GetNodalTag(self, tagName:str):
        """return the Tag (instance of the class) with name tagName.
        If the tag does not exist a new is created """
        return self.nodesTags.CreateTag(tagName,False)

    def GetNumberOfNodes(self) -> int :
        """Return the number of nodes

        Returns
        -------
        int
            the number of nodes

        Raises
        ------
        NotImplemented
            if not implemented in derived class
        """
        raise NotImplemented() # pragma: no cover

    def ComputeBoundingBox(self):
        pass

    def PrepareForOutput(self):
        pass    # pragma: no cover

    def GetElementsDimensionality(self) -> int:
        """Return the maximal dimension of the elements

        Returns
        -------
        int
            the max of all elements dimensionality
        """
        return int(np.max([ElementNames.dimension[elemtype] for elemtype in self.elements.keys() ]) )

    def IsConstantRectilinear(self): return False
    def IsRectilinear(self): return False
    def IsStructured(self): return False
    def IsUnstructured(self): return False

    def GenerateManufacturedOriginalIDs(self, offset:int=0) :
        """function to generate a valid originalids for the nodes and the data

        Parameters
        ----------
        offset : int, optional
            offset to use to start counting nodes and elements.
            This is useful for code needing numbering stating from 1, by default 0
        """
        self.originalIDNodes = np.arange(self.GetNumberOfNodes(),dtype=PBasicIndexType)
        self.originalIDNodes += offset

        counter = 0
        for key, value in self.elements.items():
            value.originalIds = np.arange(counter,counter+value.GetNumberOfElements(),dtype=PBasicIndexType)
            value.originalIds += offset
            counter += value.GetNumberOfElements()


    def WithModification(self) -> ClosingMeshAutomatically:
        """Context manager to release all the extra memory of the mesh after modification.
        this context manager call mesh.PrepareForOutput() automatically at the end of the block

        Example
        -------
        with mesh.WithModification():
            # mesh.PrepareForOutput() is called at the end of this block,
            # even if code in the block raises an exception

        """

        class ClosingMeshAutomatically():
            def __init__(self,mesh):
                self.mesh = mesh
            def __enter__(self):
                pass

            def __exit__(self, type, value, traceback):
                self.mesh.PrepareForOutput()
        return ClosingMeshAutomatically(self)

def CheckIntegrity():
    obj = MeshBase()
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
