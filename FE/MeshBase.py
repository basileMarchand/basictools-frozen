# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class Tag(object):
    def __init__(self,tagname):
        self.name = tagname
        self._id = np.empty(0,dtype=np.int)
        self.cpt =0


    def AddToTag(self,tid):
        if type(tid).__module__ == np.__name__:
            self._id.resize(self.cpt+tid.size)
            self._id[self.cpt:] = tid

            self.cpt += tid.size
        else:
            if len(self._id) <= self.cpt:
                self._id.resize(self.cpt*2+1)

            self._id[self.cpt] = tid
            self.cpt += 1

    def __len__(self):
        return self.cpt

    def Tighten(self):
        if self._id.shape[0] != self.cpt:
            self._id = np.resize(self._id, (self.cpt,))

    def SetIds(self, ids):
        self._id = np.array(ids,dtype=np.int)
        self.cpt = len(ids)

    def SetId(self, pos, ids):
        self._id[pos] = ids

    def Allocate(self,l):
        self.cpt = l
        self.Tighten()

    def GetIds(self):
        self.Tighten()
        return self._id

    def __str__(self):
        res = ''
        res  = str(self.name) + " " + str(self.cpt) + " "
        return res

class Tags(BaseOutputObject):
    def __init__(self):
        super(Tags,self).__init__()
        self.storage = []

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

    def GetNodalTag(self, tagName):
        if tagName not in self.nodesTags:
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
    tag.Allocate(1)
    tag.SetId(0,0)
    tag.AddToTag(1)
    len(tag)
    print(obj.nodesTags)
    print('toto' in obj.nodesTags)
    obj.nodesTags.RenameTag('toto','newtoto')
    obj.nodesTags.RenameTag('toto','newtoto', noError= True)

    obj.nodesTags.DeleteTags(["newtoto"])

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover