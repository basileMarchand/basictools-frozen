# -*- coding: utf-8 -*-

import numpy as np

class Tag():
        def __init__(self,tagname):
            self.name = tagname
            self.id = np.empty((0,1),dtype=np.int)
        def AddToTag(self,tid):
            self.id = np.append(self.id,tid)
            
class MeshBase(object):
    
    def __init__(self):
        super(MeshBase,self).__init__()
        self.nodesTags = {}

    def GetNumberOfNodes(self):
        raise
        
    def PrepareForOutput(self):
        pass    
      
    def IsConstantRectilinear(self): return False
    def IsRectilinear(self): return False
    def IsStructured(self): return False
    def IsUnstructured(self): return False
        