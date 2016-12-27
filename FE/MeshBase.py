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
            
        
class MeshBase(BaseOutputObject):

    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = {}

    def GetNodalTag(self, tagName):
        if not self.nodesTags.has_key(tagName):
           res = Tag(tagName)
           self.nodesTags[tagName] = res 
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