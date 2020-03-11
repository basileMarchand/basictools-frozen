# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from BasicTools.Helpers.Factory import Factory

def RegisterClass(name, classtype, constructor=None, withError = True):
    return ImplicitGeometryFactory.RegisterClass(name,classtype, constructor=constructor, withError = withError )

def Create(name,ops=None):
   return ImplicitGeometryFactory.Create(name,ops)


class ImplicitGeometryFactory(Factory):
    _Catalog = {}

    def __init__(self):
        super(ImplicitGeometryFactory,self).__init__()

def InitAllImplicitGeometry():
    import BasicTools.ImplicitGeometry.ImplicitGeometryBase
    import BasicTools.ImplicitGeometry.ImplicitGeometryObjects
    import BasicTools.ImplicitGeometry.ImplicitGeometryOperators


def CheckIntegrity():
    obj = ImplicitGeometryFactory()
    InitAllImplicitGeometry()
    return "OK"

if __name__ == '__main__': # pragma: no cover
        print(CheckIntegrity())
