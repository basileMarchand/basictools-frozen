# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.NumpyDefs import PBasicIndexType
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Tags import Tag, Tags

class MeshBase(BaseOutputObject):
    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = Tags()
        self.nodeFields = {}
        self.elemFields = {}
        """Metadata this is just a dictionary that can be used to transport
        information with the mesh, please use the class name as key of the
        object using/setting the information
        """
        self.props = {}

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

    def GetElementsOfType(self,typename):
        """
        return the element container for the typename element
        """
        return self.elements.GetElementsOfType(typename)

    def GetNamesOfElemTags(self):
        """
        return a list containing all the element tags present in the mehs
        """
        res = set()
        for ntype, data in self.elements.items():
            for tag in data.tags:
                res.add(tag.name)

        return list(res)

    def CopyProperties(self,other):
        import copy
        self.props = copy.deepcopy(self.props)


    def GetNodalTag(self, tagName):
        if tagName not in self.nodesTags:
           res = Tag(tagName)
           self.nodesTags.AddTag(res)
           return res
        else:
           return self.nodesTags[tagName]

    def GetNumberOfNodes(self):
        raise Exception()# pragma: no cover

    def ComputeBoundingBox(self):
        pass

    def PrepareForOutput(self):
        pass    # pragma: no cover

    def GetElementsDimensionality(self) -> int:
        """return the maximal dimension of the elements

        Returns
        -------
        int
            the max of all elements dimensionality
        """

        return np.max([ElementNames.dimension[elemtype] for elemtype in self.elements.keys() ])

    def IsConstantRectilinear(self): return False
    def IsRectilinear(self): return False
    def IsStructured(self): return False
    def IsUnstructured(self): return False

    def GenerateManufacturedOriginalIDs(self,offset=0):
        """
        function to generate a valid originalid data
        """
        self.originalIDNodes = np.arange(self.GetNumberOfNodes(),dtype=PBasicIndexType)
        self.originalIDNodes += offset

        counter = 0
        for key, value in self.elements.items():
           value.originalIds = np.arange(counter,counter+value.GetNumberOfElements(),dtype=PBasicIndexType)
           value.originalIds += offset
           counter += value.GetNumberOfElements()


    def WithModification(self):
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
