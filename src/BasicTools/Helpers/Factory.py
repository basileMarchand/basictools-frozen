#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
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
    def GetClass(cls,name):
        classType, classConstructor = cls._Catalog[name]
        return classType

    @classmethod
    def GetConstructor(cls,name):
        classType, classConstructor = cls._Catalog[name]
        return classConstructor

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

    @classmethod
    def PrintAvailable(cls,fullDetails=False):
        def PrintDoctring(name, obj,doc_or_type="doc",indent=0):
            if doc_or_type == "doc":
                doc = obj.__doc__
                if doc is None:
                    print(" "*indent,name, ": ")
                else:
                    print(" "*indent,name," : ", doc)
            else:
                print(" "*indent + name+ " : (", type(obj) + ")")

        for name,data in cls._Catalog.items():
            if data[1] is not None:
                if data[1].__doc__ is not None:
                    obj = data[1]
                else:
                    obj = data[0]()
            else:
                obj = data[0]()

            PrintDoctring(name,obj,indent=0)
            if fullDetails:
                for propName in obj:
                    if propName[0] == "_" : continue
                    PrintDoctring(propName,obj.__dict__[propName],"type",6)



def CheckIntegrity(GUI=False):
    fact = Factory()
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
