# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import os
from os.path import exists
from collections import defaultdict
import hashlib

import numpy as np

from BasicTools.Helpers.Tests import TestTempDir
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

def CachedResultDecorator(name=None, path=None):
    def f(function):
        return GetFunctionWithCache(function,  name=name, path=path)
    return f

def HashFunction(v):
    if type(v) is str:
        return  hashlib.sha256(v.encode())
    elif v is None:
        return HashFunction("None")
    else:
        try :
            return hashlib.sha256(v)
        except:
            import pickle 
            return HashFunction(pickle.dumps(v, fix_imports=False))


def GetFunctionWithCache(function,  name=None, path=None):
    """ Helper function to conserve the output of a funciton for later invocation
        with the same arguments. The function must be a pure function (ie. the 
        function must not rely on any persisten or internal state). The user can 
        use a dummy argument to force the function to be calculated. 
        All the argument used must support the == operator and a sha256 or pickle 
    """    
    # warning in the case lambda function are used without name
    if function.__name__ == "<lambda>" and name is None:
        print("Warning using lambda function with name=None can lead to bugs ")

    # checking if the user try to pass a already cached function
    if function.__name__.find("_withcache") > -1: 
        raise(Exception("Cant create a cached funtion of a FuncWithCache function")) 

    # default name 
    if name is None:
        name = function.__name__
    
    # default path 
    if path is None:
        path = TestTempDir.GetTempPath()

    # make sure the path exist and is created
    os.makedirs(path, exist_ok=True)


    
    def FuncWithCache(*args,**kwargs):
        # generation of the hash to create a unique filename
        hash = ""
        # for the name of the function
        hash += HashFunction(function.__name__).hexdigest()

        # for each argument
        for a in args:
            hash += HashFunction(a).hexdigest()

        # for each keyword argument
        for k,v in kwargs.items():
            hash += HashFunction(k).hexdigest()
            hash += HashFunction(v).hexdigest()
        finalhash =HashFunction(hash).hexdigest()
        
        filenameinputs = path+"Cached_inputs_" + name + "_" + finalhash + ".pickle"
        filenameoutputs = path+"Cached_outputs_" + name + "_" + finalhash + ".pickle"

        # detection of old cache and handeling hash clash       
        if os.path.exists(filenameinputs) and os.path.exists(filenameoutputs):
            try:
                from BasicTools.IO.PickleTools import LoadData
                data = LoadData(filenameinputs)
               
                if len(data.unamed) != len(args):
                    raise(Exception("Invalid number of arguments"))

                for v0,v1 in zip(data.unamed,args):
                    if np.all(v0 == v1):
                        continue
                    raise(Exception("Error on arguments"))

                if len(data.named) != len(kwargs):
                    raise(Exception("Invalid number of keyword arguments"))

                for k in data.named:
                    v0 = data.named[k]
                    v1 = kwargs[k]
                    if np.all(v0 == v1):
                        continue
                    raise(Exception("Error on keyword arguments"))

                # all ok send back the cached result
                data = LoadData(filenameoutputs)
                return data.unamed[0]
                
            except Exception as e:
                BOO().PrintDebug("error reading data from or cached data not valid (rebuilding cache) " ) 
                BOO().PrintDebug(filenameinputs)
                BOO().PrintDebug(filenameoutputs)
                pass
        else:
            BOO().PrintDebug("cache not found")    
        # do the heavy computation 
        res = function(*args,**kwargs)

        # try to save the input and the ouput data
        try: 
            from BasicTools.IO.PickleTools import SaveData
            SaveData(filenameinputs,*args,**kwargs)
            SaveData(filenameoutputs,res)
        except:
            BOO().PrintDebug("Error saving cache data.")

        return res

    # change the name of the cached function to a more explicite name 
    FuncWithCache.__name__ = function.__name__ + "_withcache"
    return FuncWithCache

def CheckIntegrity():
    cpt = 0
    
    def CountTheNumberOfExecutions(arg):
        import time
        return time.localtime()

    print(CountTheNumberOfExecutions("hola"))

    f = GetFunctionWithCache(CountTheNumberOfExecutions)
    print(f("hola") )
    print(f("hola") )


    import numpy as np
    a = np.arange(5)
    def plus1(arg):
        return arg+1

    f = GetFunctionWithCache(plus1)

    print("Original Function ", plus1(a))    
    print("First call with cache",f(a))    
    print("Second call with cache", f(a))    
    print("First call with cache new argument",f(1))    
    print("Secont call with cache new argument ",f(1))    

    def return2():
        return 2

    @CachedResultDecorator(name="superFunction")
    def return3():
        import time
        #time.sleep(3)
        return 3

    f = GetFunctionWithCache(return2)
    print("With no args return 2  ",f())    

    print("With no args return 3 ",return3())
    print("With no args return 3 ",return3.__name__)


    @CachedResultDecorator(name="superFunction")
    def returnAddwithArgs(toto, tata=4):
        import time
        #time.sleep(toto)
        print(f"runnig function toto {toto}, tata {tata}**********************************", )
        return toto+tata

    print("----------------------------------------------------")
    print("With args return (3) ->",returnAddwithArgs(3))
    print("----------------------------------------------------")


    print("With args return (3,tata=4) ->",returnAddwithArgs(3,tata=4))
    print("----------------------------------------------------")
    print("With args return (3,tata=4) ->",returnAddwithArgs(3,tata=4))
    #print("----------------------------------------------------")
    #print("With args return (3,tata=6) ->",returnAddwithArgs(3,tata=6))
    print("----------------------------------------------------")
    #print(return2.__name__)
    #print(GetFunctionWithCache(lambda x: 5).__name__)
    return "ok"
    

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
