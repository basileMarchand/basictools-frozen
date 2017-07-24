# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import BasicTools.FE.ElementNames as ElementsNames
from BasicTools.FE.Tri3 import Tri3 as Tri3

def GetElementFromName(name):
    if name == ElementsNames.Triangle_3:
        return Tri3()
#    elif name == ElementsNames.:
#        return ()
    raise #pragma: no cover

def CheckIntegrity():
    GetElementFromName(ElementsNames.Triangle_3)
    return "OK"
