# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

allElements = object()
bulkElements = object() # element of dimensionality for meshdim (Tetra in 3D, Tri in 2D)
borderElements = object() # element of dimensionality for meshdim-1 (Tri in 3D, Edge in 2D)
borderborderElements = object() # element of dimensionality for meshdim-2 (Edge in 3D, Point in 2D)

class Tag(object):
    def __init__(self,tagname):
        self.name = tagname
        self._id = np.empty(0,dtype=np.int)
        self.cpt =0


    def AddToTag(self,tid):
        if type(tid).__module__ == np.__name__:
            self._id = np.resize(self._id, (self.cpt+tid.size,))
            self._id[self.cpt:] = tid

            self.cpt += tid.size
        else:
            if len(self._id) <= self.cpt:
                self._id = np.resize(self._id, (self.cpt*2+1,))

            self._id[self.cpt] = tid
            self.cpt += 1

    def __len__(self):
        return self.cpt

    def Tighten(self):
        if self._id.shape[0] != self.cpt:
            self._id = np.resize(self._id, (self.cpt,))

    def RemoveDoubles(self):
        self.Tighten()
        self.SetIds(self._id)

    def SetIds(self, ids):
        self._id = np.unique(np.asarray(ids,dtype=np.int))
        self.cpt = len(self._id)

    def SetId(self, pos, ids):
        self._id[pos] = ids

    def Allocate(self,l):
        self.cpt = l
        self.Tighten()

    def GetIds(self):
        self.Tighten()
        return self._id

    def GetIdsAsMask(self,totalNumberOfObjects=None,output=None,erase=True):
        """
        .. py:classmethod:: GetIdsAsMask(self,totalNumberOfObjects,out=None,erase=True)

        Get mask of tag items

        :param int totalNumberOfObjects: total number of Objects to allocate the output
        :param numpy.ndarray out: output array for inplace work
        :param bool erase: option to erase the output before asignment
        :return: mask items contained in tag (True,False)
        :rtype: numpy.ndarray
        """

        self.Tighten()

        if output is None:
            output = np.zeros(totalNumberOfObjects,dtype=bool)
        else:
            if erase :
                output.fill(False)

        output[self._id] = True

        return output

    def Merge(self,other=None,ids=None):
        if other is not None:
            self.SetIds(list(set().union(self.GetIds(),other.GetIds())))

        if ids is not None:
            self.SetIds(list(set().union(self.GetIds(),ids)))

    def __str__(self):
        res = ''
        res  = str(self.name) + " " + str(self.cpt) + " "
        return res

class Tags(BaseOutputObject):
    def __init__(self):
        super(Tags,self).__init__()
        self.storage = []

    def Tighten(self):
        for tag in self:
            tag.Tighten()

    def AddTag(self,item):
        if item.name in self:
            raise Exception("Cant add the tag two times!!")# pragma: no cover

        self.storage.append(item)
        return item

    def DeleteTags(self,tagNames):
        for name in tagNames:
            for i in range(len(self.storage)):
                if self.storage[i].name == name:
                    self.storage.pop(i)
                    break

    def CreateTag(self,name, errorIfAlreadyCreated=True):

        if name in self:
            if errorIfAlreadyCreated :
                raise Exception("Tag name '"+name+"' already exist")# pragma: no cover
            else:
                return self[name]
        else:
            return self.AddTag(Tag(name))

    def RenameTag(self,name,newName,noError= False):
        if name in self:
            self[name].name = newName
        else:
            if noError: return
            raise Exception("Tag '" + str(name) + "' does not exist")# pragma: no cover

    def RemoveEmptyTags(self):
        for tagname in  list(self.keys()):
            tag = self[tagname]
            if tag.cpt == 0:
                self.storage.remove(tag)

## function to act like a dict

    def keys(self):
        return [ x.name for x in self.storage ]

    def __contains__(self, k):
        for item in self.storage:
            if item.name == k:
                return True
        return False

    def __getitem__(self,key):
        for item in self.storage:
            if item.name == key:
                return item
        raise KeyError("Tag '"+ str(key) + "' not found")# pragma: no cover
        #return self.storage[k]

    def __iter__(self):
        return iter(self.storage)

    def __len__(self):
        return len(self.storage)

    def __str__(self):
        res = ''
        res  = str(list(self.keys()))
        return res

class MeshBase(BaseOutputObject):
    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = Tags()
        self.nodeFields = {}
        self.elemFields = {}
        """Metadata this is just a dictionary that can be used to tranport
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

    def IsConstantRectilinear(self): return False
    def IsRectilinear(self): return False
    def IsStructured(self): return False
    def IsUnstructured(self): return False

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
    tag = obj.GetNodalTag("toto")
    tag.Allocate(1)
    tag.SetId(0,0)
    tag.AddToTag(1)
    len(tag)

    tagII = obj.GetNodalTag("TagII")
    tagII.Allocate(2)
    tagII.SetId(0,1)
    tagII.SetId(1,2)

    tag.Merge(tagII)
    print(tag.GetIds())


    print(obj.nodesTags)
    print('toto' in obj.nodesTags)
    obj.nodesTags.RenameTag('toto','newtoto')
    obj.nodesTags.RenameTag('toto','newtoto', noError= True)

    obj.nodesTags.DeleteTags(["newtoto"])

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
