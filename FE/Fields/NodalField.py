# -*- coding: utf-8 -*-
import numpy as np


class NodalField(object):
    def __init__(self,name=None,mesh=None,space=None,numbering=None,data=None):
        super(NodalField,self).__init__()
        self.name = name
        self.data = data
        self.mesh = mesh
        self.space = space
        self.numbering = numbering

    def Allocate(self,val=0):
        if val == 0:
            self.data = np.zeros(self.numbering["size"],dtype=np.float)
        else:
            self.data = np.ones(self.numbering["size"],dtype=np.float)*val

#    def GetValueAt(self,id):
#        return self.data[id]

    def GetValueAtIP(self,elemtype,el,ip):
        sp = self.space[elemtype]
        num = self.numbering[elemtype][el,:]
        return sp.Eval_FieldI(ip,self.data[num],None,None,der=-1)

def CheckIntegrity():
    obj = NodalField("temp")
    return "ok"
