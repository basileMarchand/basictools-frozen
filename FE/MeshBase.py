# -*- coding: utf-8 -*-

class MeshBase(object):
    
    def __init__(self):
          pass
      
    def isConstantRectilinear(self): return False
    def isRectilinear(self): return False
    def isStructured(self): return False
    def isUnstructured(self): return False