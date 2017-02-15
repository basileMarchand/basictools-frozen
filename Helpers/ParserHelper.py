  # -*- coding: utf-8 -*-
import numpy as np

def Read(string, inout):
    if type(inout).__module__ == np.__name__:
        return ReadVector(string, inout.dtype)
    else:
        return ReadScalar(string, type(inout) )

def ReadScalar(string,dtype):
    if dtype is bool:
        return ReadBool(string)

    if type(string) is str:
        return dtype(string.lstrip().rstrip())
    else:
        return dtype(string)

def ReadVector(string,dtype):

    tmp = string.lstrip().rstrip()

    return np.array([ ReadScalar(x,dtype) for x in tmp.split()] ,dtype=dtype )

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

    if tmp == "true" or tmp == "yes":
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

def TestFunction(func,string):
    val = func(string)
    print( str(string) + " : " +str(val))

def CheckIntegrity():

    TestFunction(ReadBool,"true")
    TestFunction(ReadBool,True)
    TestFunction(ReadBool,"false")
    TestFunction(ReadBool,"0")
    TestFunction(ReadBool,"1")
    TestFunction(ReadBool,"1.1")
    TestFunction(ReadBool," no")
    TestFunction(ReadBool,"YES ")

    TestFunction(ReadBools,"YES no 2 1 0 True FALSe ")

    TestFunction(ReadInt,"24" )
    TestFunction(ReadInt,24.0 )
    TestFunction(ReadInts,"1 2 3 ")

    TestFunction(ReadFloat,"3.14159" )
    TestFunction(ReadFloat,3.14159 )
    TestFunction(ReadFloats,"1 2 3 ")

    # this call must fail
    try:
        ReadBool("toto")
        raise # pragma: no cover
    except:
        pass

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
