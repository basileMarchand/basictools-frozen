# -*- coding: utf-8 -*-
import numpy as np

from OTTools.Helpers.BaseOutputObject import BaseOutputObject

class FeaBase(BaseOutputObject):
     def __init__(self,dim=3, size= 1):
        super(FeaBase,self).__init__()

def deleterowcol(A, delrow, delcol, fixedValues ):
    # Assumes that matrix is in symmetric csc form !

    rhs = A*fixedValues
    #keep = np.delete (np.arange(0, m), delrow)
    A = A[np.logical_not(delrow) , :]
    #keep = np.delete (np.arange(0, m), delcol)
    A = A[:, np.logical_not(delcol)]

    return [A, rhs]



def CheckIntegrity():
    return "ok"