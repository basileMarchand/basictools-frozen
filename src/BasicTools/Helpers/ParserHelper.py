# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" Functions to parse text into different types
"""
import numpy as np
import math

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.NumpyDefs import PBasicFloatType


class LocalVariables(BaseOutputObject):
    def __init__(self, prePostChars=("{","}")):
        super(LocalVariables,self).__init__()
        self.prePostChars = prePostChars
        self.variables = {}

    def SetVariable(self,key,value):
        self.variables[str(key)] = str(value)

    def UnsetVariable(self,key):
        self.variables.pop(str(key))

    def Apply(self,string):
        if type(string) != str: return string
        for key,value in self.variables.items():
            string= string.replace(self.prePostChars[0]+key+self.prePostChars[1],value)
        return string

globalDict = LocalVariables()

def AddCommonConstants():
    from scipy import constants
    siPrefixes = {"yotta":constants.yotta,
                  "zetta":constants.zetta,
                  "exa"  :constants.exa,
                  "peta" :constants.peta,
                  "tera" :constants.tera,
                  "giga" :constants.giga,
                  "mega" :constants.mega,
                  "kilo" :constants.kilo,
                  "hecto":constants.hecto,
                  "deka" :constants.deka,
                  "deci" :constants.deci,
                  "centi":constants.centi,
                  "milli":constants.milli,
                  "micro":constants.micro,
                  "nano" :constants.nano,
                  "pico" :constants.pico,
                  "femto":constants.femto,
                  "atto" :constants.atto,
                  "zepto":constants.zepto,
                  "yocto":constants.yocto,

                  "pi"   : constants.pi,
                  "GradToDegree": (180/constants.pi),
                  "DegreeToGrad": (constants.pi/180)}

    for name,value in siPrefixes.items():
        AddToGlobalDictionary(name,value)


def AddToGlobalDictionary(key,value):
    globalDict.SetVariable(key,value)

def RemoveFromGlobalDictionary(key):
    globalDict.UnsetVariable(key)

def ApplyGlobalDictionary(string):
    return globalDict.Apply(string)


def Read(inputData, inout):
    if type(inout).__module__ == np.__name__:
        return ReadVector(inputData, inout.dtype.type)
    else:
        if inout is None:
            inputData = ApplyGlobalDictionary(inputData)
            return inputData
        else:
            return ReadScalar(inputData, type(inout) )

def ReadScalar(inputData, inputType):
    inputData = ApplyGlobalDictionary(inputData)
    if inputType is bool:
        return ReadBool(inputData)

    try:
        if type(inputData) is str:
            return inputType(inputData.lstrip().rstrip())
        else:
            return inputType(inputData)
    except ValueError:
        from BasicTools.Containers.SymExpr import SymExprBase
        try:
            return inputType(SymExprBase(inputData).GetValue())
        except TypeError as e:
            print("Cannot convert expression to  {}, expr : {}".format(inputType,inputData))
            raise(e)

def ReadVector(string, dtype):

    if isinstance(string,(list,tuple,np.ndarray)) :
        return np.array([ ReadScalar(x,dtype) for x in string] ,dtype=dtype )
    else:
        tmp = string.lstrip().rstrip()
        return np.array([ ReadScalar(x,dtype) for x in tmp.split()] ,dtype=dtype )

def ReadString(string, d=None):
    if d is not None:
         string = string.replace(d)
    string = ApplyGlobalDictionary(string)
    return str(string)

def ReadStrings(string):
    return ReadVector(string,str)

def ReadFloat(string):
    return ReadScalar(string,float)

def ReadFloats(string):
    return ReadVector(string,float)

def ReadInt(string):
    return ReadScalar(string,int)

def ReadInts(string):
    return ReadVector(string,int)

def ReadBool(string):
    if type(string) is bool:
        return bool(string)

    if type(string) is not str:
        return bool(string)

    string = ApplyGlobalDictionary(string)

    tmp = string.lstrip().rstrip().lower()
    if tmp == "false" or tmp == "no":
        return False

    if tmp == "true" or tmp == "yes" or tmp=="on":
        return True

    try:
        d = ReadInt(tmp)
        return bool(d)
    except:
        pass

    try:
        d = ReadFloat(tmp)
        return bool(d)
    except:
        pass

    raise ValueError("cant convert '" + string +"' into bool")

def ReadBools(string):
    return ReadVector(string,bool)

def ReadVectorXY(string, normalised=False):
    res = ReadFloats(string)
    if len(res) != 2:
        raise
    if normalised:
        res /= np.linalg.norm(res)
    return res

def ReadVectorXYZ(string,normalised=False):
    res = ReadFloats(string)
    if len(res) != 3:
        raise
    if normalised:
        res /= np.linalg.norm(res)
    return res

def ReadVectorPhiThetaMag(string,normalised=False):
    res = ReadFloats(string)
    if len(res) == 3:
        phi,theta,mag = res
    else:
        phi,theta = res
        mag = 1.

    if normalised:
        mag = 1.
        if len(res) == 3:
            print("Warning: mag ignored")

    phi = phi*math.pi/180.
    theta = theta*math.pi/180.
    res = np.array([math.sin(phi)*math.cos(theta), math.sin(phi)*math.sin(theta), math.cos(phi)])
    return res*mag


def ReadProperties(data, props ,obj_or_dic,typeConversion=True):
    if props is None:
        if type(obj_or_dic) == dict:
            props = obj_or_dic.get("_parseProps", None)
        else:
            props = getattr( obj_or_dic, "_parseProps", None)

    if props is None:
        return
    try:
      for prop in props:
        if prop in data:

           theSetter = getattr( obj_or_dic, "Set"+prop[0].upper()+ str(prop[1:]), None)
           inputData = ApplyGlobalDictionary(data[prop])
           if theSetter is None:
              #print(obj.__dict__)
              #conversion only if the target type is different of None
              #try:


                  if type(obj_or_dic) == dict:
                     if typeConversion and (obj_or_dic[prop] is not None) :
                         obj_or_dic[prop] = Read(inputData,obj_or_dic[prop])
                     else:
                         obj_or_dic[prop] = inputData
                  else:
                     if typeConversion and (obj_or_dic.__dict__[prop] is not None) :
                         obj_or_dic.__dict__[prop] = Read(inputData,obj_or_dic.__dict__[prop])
                     else:
                         obj_or_dic.__dict__[prop] = inputData

              #except:
              #    raise (ValueError("Error setting  '"+str(prop)+"'  to object of type " + str(type(obj_or_dic)) ) )
           else:
                theGetter = getattr( obj_or_dic, "Get"+prop[0].upper()+ str(prop[1:]), None)
                if theGetter is not None and typeConversion:
                    theSetter(Read(inputData,theGetter()))
                else:
                    theSetter(inputData)
    except KeyError as e:
        print(" object of type " +str(type(obj_or_dic)) + " does not have attribute {0}: ".format( str(e) ))
        raise



def TestFunction(func,string,correctVal):
    val = func(string)


    print(str(func.__name__) + "("+  str(string) + ") = " +  str(val) + " =? " +str(correctVal))

    if type(val) != type(correctVal):
        raise Exception("returned values does not have the correct type")

    if np.any(val != correctVal):
        raise Exception("returned value does not match")

def CheckIntegrity():

    TestFunction(ReadBool,"true",True)
    TestFunction(ReadBool,True,True)
    TestFunction(ReadBool,"false",False)
    TestFunction(ReadBool,"0",False)
    TestFunction(ReadBool,"1",True)
    TestFunction(ReadBool,"1.1",True)
    TestFunction(ReadBool," no",False)
    TestFunction(ReadBool,"YES ",True)

    TestFunction(ReadBools,"YES no 2 1 0 True FALSe ", np.array([True,False,True,True,False,True,False]))

    TestFunction(ReadInt,"24",int(24) )
    TestFunction(ReadInt,24.0 ,int(24))
    TestFunction(ReadInts,"1 2 3 ", np.array([1,2,3]))

    TestFunction(ReadFloat,"3.14159", 3.14159 )
    TestFunction(ReadFloat,"3.14159*10/5/2", 3.14159 )
    TestFunction(ReadFloat,"exp(pi)", math.exp(math.pi) )
    TestFunction(ReadFloat,3.14159, 3.14159 )
    TestFunction(ReadFloats,"1 2 3 ", np.array([1,2,3],dtype=float))
    TestFunction(ReadFloats,"1 4/2 9/3 ", np.array([1,2,3],dtype=float))
    TestFunction(ReadInts,"1 4/2 9/3 ", np.array([1,2,3],dtype=int))
    AddCommonConstants()
    TestFunction(ReadFloat,"{pico}",1e-12)
    TestFunction(ReadFloat,"{pi}*{GradToDegree}*{DegreeToGrad}", 3.1415926535897927)
    TestFunction(ReadFloat,"arctan(1)*{GradToDegree}",45.)

    TestFunction(ReadStrings,"Hola Chao", np.array(["Hola","Chao"],dtype=str))

    # this call must fail
    try:
        ReadBool("toto")
        raise # pragma: no cover
    except:
        pass

    #### Reading data into a class of dictionary with type conversion
    data = {"monInt":"2.2","monFloat":"3.14159"}

    class Options():
        def __init__(self):
            self.monInt = 1
            self.monFloat = 0.1
        def SetMonFloat(self,data):
            self.monFloat = ReadFloat(data)


    ops = Options()
    ReadProperties(data,data.keys(),ops)
    if type(ops.monInt) != int or  ops.monInt != 2:
        raise
    if type(ops.monFloat) != float or  ops.monFloat != 3.14159:
        raise

    outputData = {"monInt":00,"monFloat":00.00}
    ReadProperties(data,data.keys(),outputData)

    if type(outputData['monInt']) != int or  outputData['monInt'] != 2:
        raise
    if type(outputData['monFloat']) != float or  outputData['monFloat'] != 3.14159:
        raise


    AddToGlobalDictionary("FILENAME","toto.xml")

    print(ReadString("my FILENAME is '{FILENAME}' "))
    RemoveFromGlobalDictionary("FILENAME")
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
