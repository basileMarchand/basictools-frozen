# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:52:41 2016

@author: d584808
"""

import OTTools.FE.ElementNames as ElementsNames
from OTTools.FE.Tri3 import Tri3 as Tri3

def GetElementFromName(name):
    if name == ElementsNames.Triangle_3:
        return Tri3()
#    elif name == ElementsNames.:
#        return ()
    raise        
    
def CheckIntegrity():
    GetElementFromName(ElementsNames.Triangle_3)
    return "OK"
    