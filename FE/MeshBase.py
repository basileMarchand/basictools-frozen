# -*- coding: utf-8 -*-

import numpy as np
from OTTools.Helpers.BaseOutputObject import BaseOutputObject

class Tag():
        def __init__(self,tagname):
            self.name = tagname
            self.id = np.empty((0,1),dtype=np.int)
        def AddToTag(self,tid):
            self.id = np.append(self.id,tid)
        def __len__(self):
            return len(self.id)

class MeshBase(BaseOutputObject):

    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = {}

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
    tag = Tag("toto")
    tag.AddToTag(0)
    len(tag)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover