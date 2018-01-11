# -*- coding: utf-8 -*-
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class ConstantField(BaseOutputObject):
    def __init__(self,name,val):
        super(ConstantField,self).__init__()
        self.val = val
        self.name = name


    # signature for elementtype, element number, ipoint
    def GetValueAtIP(self,name=None,el=None,ip=None):
        return self.val

def CheckIntegrity(GUI=False):
    obj = ConstantField("toto",5)

    if obj.GetValueAtIP("",0,0) == 5:
        return "OK"
    else:
        return "Not OK"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
