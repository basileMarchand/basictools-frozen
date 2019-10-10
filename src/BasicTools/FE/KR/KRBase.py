# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

class KRBase(BOO):
    def __init__(self):
        super(KRBase,self).__init__()
        self.args = []
        self.on = []
        self.blockDirections = [False,False,False]


    def AddArg(self,name):
        self.args.append(name)
        return self

    def On(self,zone):
        if type(zone) is list:
            self.on.extend(zone)
        else:
            self.on.append(zone)
        self.on = list(set(self.on))
        return self

    def Fix0(self,val=True):
        self.blockDirections[0] = val
        return self
    def Fix1(self,val=True):
        self.blockDirections[1] = val
        return self
    def Fix2(self,val=True):
        self.blockDirections[2] = val
        return self



def CheckIntegrity(GUI=False):
    obj = KRBase()
    obj.AddArg("u").On("Z0").Fix0().Fix1(False).Fix2(True)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
