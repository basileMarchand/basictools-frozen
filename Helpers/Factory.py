# -*- coding: utf-8 -*-

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.Helpers.ParserHelper as PH
"""
to use this factory you must create a class like this

def RegisterClass(name, classtype, constructor=None, withError = True):
    return ImplicitGeometryFactory.RegisterClass(name,classtype, constructor=constructor, withError = withError )

def Create(name,ops=None):
   return ImplicitGeometryFactory.Create(name,ops)

class ImplicitGeometryFactory(Factory):
    _Catalog = {}

    def __init(self):
        super(ImplicitGeometryFactory,self).__init__()

and use the functions to register and create objects

RegisterClass("Offset",ImplicitGeometryOffset,CreateImplicitGeometryOffset)
or (with a special constructor)
RegisterClass("Offset",ImplicitGeometryOffset)

and to create a class
res = Create(name,ops)


"""

class Factory(BaseOutputObject):

    def __init__(self):
        super(Factory,self).__init__()

    @classmethod
    def keys(cls):
        return cls._Catalog.keys()

    @classmethod
    def RegisterClass(cls, name, classtype, constructor=None, withError = True):
        #cls().PrintDebug(str(name) + " -> " +  str(classtype) )
        if name in cls._Catalog and withError:
           raise (Exception ("Class "+ str(name) +" already in the catalog") )
        cls._Catalog[name] = (classtype,constructor)

    @classmethod
    def Create(cls,name,ops=None,propertiesAssign=True):

        res = None
        if name in cls._Catalog:
           classType, classConstructor = cls._Catalog[name]
           #cls().PrintDebug(str(classType)+ " : " + str(ops) )
           if classConstructor is None:
               try:
                   res = classType()
               except Exception as e:
                   print(classType)
                   print("Error creating class (" +name+ "):" + str(classType) + ". ")
                   raise(e)
               if propertiesAssign:
                   PH.ReadProperties(ops, ops, res)
           else:
               if ops is None:
                   return classType()
               res = classConstructor(ops)
           return res

        raise(Exception("Unable to create object of type " + str(name) +"\n Possible object are :"+ str(list(cls._Catalog.keys()))  ))# pragma: no cover

def CheckIntegrity(GUI=False):
    fact = Factory()
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
