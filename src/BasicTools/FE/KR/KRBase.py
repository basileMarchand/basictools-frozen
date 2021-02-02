# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
from BasicTools.Linalg.ConstraintsHolder import ConstraintsHolder
from BasicTools.FE.ProblemData import Transform


class KRBase(BOO):
    def __init__(self):
        super(KRBase,self).__init__()
        self.args = []
        self.on = []

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

    def _GetConstraintHolder(self,CH):
        if CH is None:
            return ConstraintsHolder()
        return CH

    def _ComputeOffsets(self, fields):

        cpt = 0
        offsets = []
        fieldOffsets = { }

        for field in fields:
            offsets.append(cpt)
            fieldOffsets[field.name] = cpt
            cpt += field.numbering["size"]

        return offsets, fieldOffsets

class KRBaseScalar(KRBase):
    def __init__(self):
        super(KRBaseScalar,self).__init__()
        self.value = lambda x: 0.0

    def SetFunction(self,func):
        self.value = func
        return self

    def SetValue(self,val):
        self.SetFunction(lambda x: val)
        return self

class KRBaseVector(KRBase):
    def __init__(self):
        super(KRBaseVector,self).__init__()
        self.blockDirections = [False,False,False]
        self.constraintDiretions = "Target" # ["Global","Local","Origin","Target","Normal"]

        self.originSystem = Transform()
        self.targetSystem = Transform()

    def Fix(self,direc,val=True):
        self.blockDirections[direc] = val
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

    def GetConstrainedDirections(self,pos=None,direction=None):
        res = []

        for x,y in zip([0,1,2],self.blockDirections ):
            if y:
                res.append(self.GetDirections(x,pos,direction))
        return res


    def GetDirections(self,i,pos=None,direction=None):
        if self.constraintDiretions == "Global":
            res = np.zeros(3)
            res[i] = 1
            return res
        elif self.constraintDiretions == "Origin":
            return self.originSystem.GetOrthoNormalBase().GetDirection(i,pos,direction)
        elif self.constraintDiretions == "Target":
            return self.targetSystem.GetOrthoNormalBase().GetDirection(i,pos,direction)
        elif self.constraintDiretions == "Local":
            if i == 0:
                return direction / np.linalg.norm(direction)
            #Create a local base based on the direction vector
            # need a second vector to generate a connsistent base(???)

            raise
        else:
            raise

def CheckIntegrity(GUI=False):

    obj = KRBaseScalar()
    obj.AddArg("u").On("Z0")
    obj.SetValue(1.0)

    obj = KRBaseVector()
    obj.AddArg("u").On("Z0").Fix0().Fix1(False).Fix2(True)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
