# -*- coding: utf-8 -*-
""" Functions to parse text into different types
"""
__author__ = "Felipe Bordeu"

import numpy as np
import math

def Read(inputData, inout):
    if type(inout).__module__ == np.__name__:
        return ReadVector(inputData, inout.dtype.type)
    else:
        if inout is None:
            return inputData
        else:
            return ReadScalar(inputData, type(inout) )

def ReadScalar(inputData,inputtype):

    if inputtype is bool:
        return ReadBool(inputData)

    if type(inputData) is str:
        return inputtype(inputData.lstrip().rstrip())
    else:
        return inputtype(inputData)

def ReadVector(string,dtype):

    if isinstance(string,list) or  isinstance(string,np.ndarray):
        return np.array([ ReadScalar(x,dtype) for x in string] ,dtype=dtype )
    else:
        tmp = string.lstrip().rstrip()
        return np.array([ ReadScalar(x,dtype) for x in tmp.split()] ,dtype=dtype )

def ReadString(string):
    return str(string)

def ReadStrings(string):
    return ReadVector(string,str)

def ReadFloat(string):
    return ReadScalar(string,float)

def ReadFloats(string):
    return ReadVector(string,np.float)

def ReadInt(string):
    return ReadScalar(string,np.int)

def ReadInts(string):
    return ReadVector(string,np.int)

def ReadBool(string):
    if type(string) is bool:
        return bool(string)

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

    raise Exception("cant conver '" + string +"' into bool")

def ReadBools(string):
    return ReadVector(string,bool)

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

    phi = phi*math.pi/180.
    theta = theta*math.pi/180.
    res = np.array([math.sin(phi)*math.cos(theta), math.sin(phi)*math.sin(theta), math.cos(phi)])
    return res*mag


def ReadProperties(data, props ,obj,typeConversion=True):
    if props is None:
        props = getattr( obj, "_parseProps", None)
    if props is None:
        return
    try:
      for prop in props:
        if prop in data:
           theSetter = getattr( obj, "Set"+prop[0].upper()+ str(prop[1:]), None)
           if theSetter is None:
              #print(obj.__dict__)
              #conversion only if the target type is different of None
              try:
                  if typeConversion and (obj.__dict__[prop] is not None) :
                      obj.__dict__[prop] = Read(data[prop],obj.__dict__[prop])
                  else:
                      obj.__dict__[prop] = data[prop]
              except:
                  print(obj )
                  print(prop)
                  print(data)
                  raise (ValueError("Error setting  '"+str(prop)+"'  to object of type " + str(type(obj)) ) )
           else:
              theSetter(data[prop])
    except KeyError as e:
        print(" object of type " +str(type(obj)) + " does not have atribute {0}: ".format( str(e) ))
        raise



def TestFunction(func,string,correctVal):
    val = func(string)
    print( str(string) + " : " +str(val))

    if type(val) != type(correctVal) or np.any(val != correctVal):
        raise

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
    TestFunction(ReadFloat,3.14159, 3.14159 )
    TestFunction(ReadFloats,"1 2 3 ", np.array([1,2,3],dtype=float))

    TestFunction(ReadStrings,"Hola Chao", np.array(["Hola","Chao"],dtype=str))

    # this call must fail
    try:
        ReadBool("toto")
        raise # pragma: no cover
    except:
        pass

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
