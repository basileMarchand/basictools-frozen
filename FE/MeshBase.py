# -*- coding: utf-8 -*-

import numpy as np
from OTTools.Helpers.BaseOutputObject import BaseOutputObject

class Tag():
    def __init__(self,tagname):
        self.name = tagname
        self.id = np.empty((0,1),dtype=np.int)
        self.cpt =0

    def AddToTag(self,tid):
        if len(self.id) <= self.cpt:
            self.id.resize((self.cpt*2+1,1))

        self.id[self.cpt] = tid
        self.cpt += 1

    def __len__(self):
        return self.cpt

    def tighten(self):
        self.id = np.resize(self.id, (self.cpt,))

    def SetIds(self, ids):
        self.id = ids
        self.cpt = len(ids)

class Tags(BaseOutputObject):
    def __init__(self):
        super(Tags,self).__init__()
        self.storage = []

    def AddTag(self,item):
        if self.has_key(item.name):
            raise Exception("Cant add the tag two times!!")

        self.storage.append(item)
        return item

    def CreateTag(self,name, errorIfAlreadyCreated=True):

        if self.has_key(name):
            if errorIfAlreadyCreated :
                raise Exception("Tag already exist")
            else:
                return self[name]
        else:
            return self.AddTag(Tag(name))


## function to act like a dict

    def keys(self):
        return [ x.name for x in self.storage ]

    def has_key(self, k):
        for item in self.storage:
            if item.name == k:
                return True
        return False

    def __getitem__(self,key):
        for item in self.storage:
            if item.name == key:
                return item
        raise KeyError("Tag '"+ str(key) + "' not found")
        #return self.storage[k]

    def __iter__(self):
        return iter(self.storage)

class MeshBase(BaseOutputObject):

    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = Tags()
        #self.nodesTags = {}

    def GetNodalTag(self, tagName):
        if not self.nodesTags.has_key(tagName):
           res = Tag(tagName)
           self.nodesTags.AddTag(res)
           return res
        else:
           return self.nodesTags[tagName]

    def GetNumberOfNodes(self):
        raise Exception()# pragma: no cover

    def PrepareForOutput(self):
        pass    # pragma: no cover

    def IsConstantRectilinear(self): return False
    def IsRectilinear(self): return False
    def IsStructured(self): return False
    def IsUnstructured(self): return False

def CheckIntegrity():
    obj = MeshBase()
    tag = obj.GetNodalTag("toto")
    tag.AddToTag(0)
    len(tag)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover