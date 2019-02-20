# -*- coding: utf-8 -*-
import numpy as np
from BasicTools.Helpers.TextFormatHelper import TFormat


class FEField(object):
    def __init__(self,name=None,mesh=None,space=None,numbering=None,data=None):
        super(FEField,self).__init__()
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

    def __str__(self):
        TFormat.II()
        res =  TFormat.InBlue("NodalField")+"\n"
        if self.name is not None:
          res += TFormat.GetIndent()
          res += TFormat.InGreen("Name : ") + self.name
        return res

    def GetPointRepresentation(self,fillvalue=0.):

        if fillvalue==0.:
            res = np.zeros(self.mesh.GetNumberOfNodes(),dtype=float)
        else:
            res = np.ones(self.mesh.GetNumberOfNodes(),dtype=float)*fillvalue

        res[self.numbering["doftopointLeft"]] = self.data[self.numbering["doftopointRight"]]

        return res

    def SetDataFromPointRepresentation(self,userdata, fillvalue=0.):
        if fillvalue==0.:
           self.data = np.zeros(self.numbering["size"])
        else:
           self.data = np.ones(self.numbering["size"])*fillvalue

        self.data[self.numbering["doftopointRight"]] = userdata[self.numbering["doftopointLeft"]]

    def __repr__(self):
        res = "FEField " + self.name
        return res

def CheckIntegrity():
    obj = FEField("temp")
    return "ok"
