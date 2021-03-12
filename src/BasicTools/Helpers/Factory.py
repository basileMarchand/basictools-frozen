# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


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
    def AllEntries(cls):
        return cls._SetCatalog

    @classmethod
    def RegisterClass(cls, name, classtype, constructor=None, withError = True):
        #cls().PrintDebug(str(name) + " -> " +  str(classtype) )
        if name in cls._Catalog and withError:
           raise (Exception ("Class "+ str(name) +" already in the catalog") )
        cls._Catalog[name] = (classtype,constructor)
        if hasattr(cls,"_SetCatalog"):
            cls._SetCatalog.add( (name,classtype,constructor))


    @classmethod
    def GetClass(cls,name):
        classType, classConstructor = cls._Catalog[name]
        return classType

    @classmethod
    def GetConstructor(cls,name):
        classType, classConstructor = cls._Catalog[name]
        return classConstructor

    @classmethod
    def GetAvailablesFor(cls, name):
        return [( obj,const) for key,obj,const in cls._SetCatalog if key == name]

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
               res = classConstructor(ops)
           return res

        raise(Exception("Unable to create object of type " + str(name) +"\n Possible object are :"+ str(list(cls._Catalog.keys()))  ))

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
                print(" "*indent + name+ " : (", str(type(obj)) + ")")

        for name,cl,cons in cls._SetCatalog:

            if cons is not None:
                if cons.__doc__ is not None:
                    obj = cons
                else:
                    obj = cl()
            else:
                obj = cl()
            print("---------------------------------------------------------")
            if (cl,cons) in cls._Catalog.values():
                print(" vvvvvv  default  for '"+name+"' vvvvv ")
            PrintDoctring(name,obj,indent=0)
            if fullDetails:
                if hasattr(obj, '__dict__') and  hasattr(obj.__dict__, '__iter__'):
                    for propName in obj.__dict__:
                        if propName[0] == "_" : continue
                        PrintDoctring(propName,obj.__dict__[propName],"type",6)



def CheckIntegrity(GUI=False):
    class DummyFactory(Factory):
        _Catalog = {}
        _SetCatalog = set()
        def __init__(self):
            super(DummyFactory,self).__init__()

    fact = DummyFactory()
    print(fact.keys())
    print(fact.AllEntries())
    fact.RegisterClass("str",str)
    ok = True
    fact.RegisterClass("str",str,withError=False)

    try:
        fact.RegisterClass("str",str)
        ok = False#pragma: no cover
    except:
        pass
    assert ok

    print(fact.GetClass("str"))
    print(fact.GetConstructor("str"))
    print(fact.GetAvailablesFor("str"))

    print(fact.Create("str"))

    ok = True
    try:
        print(fact.Create("Million Dollars"))
        ok = False#pragma: no cover
    except:
        pass
    assert ok

    def DummyConstructor(ops):
        """documentation for this dummy constructor
        """
        return str(round(ops))

    class DummyClass():
        def __init__(self):
            self.str = ""
    def DummyWithoutDocumentation(ops):
        return str(round(ops))

    fact.RegisterClass("rounded_str",str,DummyConstructor )
    print(fact.Create("rounded_str",ops=3.5))
    fact.RegisterClass("rounded_str",DummyClass,DummyWithoutDocumentation,withError=False )
    print(fact.Create("rounded_str",ops=3.5))

    fact.PrintAvailable(fullDetails=True)




    ok = True
    try:
        a = 5
        fact.RegisterClass("numpy",a)
        print(fact.Create("numpy"))
        ok = False#pragma: no cover
    except:
        pass
    assert ok

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
