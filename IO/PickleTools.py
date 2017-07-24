# -*- coding: utf-8 -*-
import pickle as __pickle

__author__ = "Felipe Bordeu"

class IOHelper:
    """helper Class that represent the data from a file  """

    def __init__(self,data):
        self.unamed = data[0]
        self.named = data[1]


    def __str__(self):
        res  = " named : "  + str(self.named)  + "\n"
        res += " unamed : " + str(self.unamed) + "\n"
        return res

def SaveData(filename,  *argv,**kwargs):
    """Save the variables into the disk and return 0 if all ok

       Save variables into the disk, you can use unamed or named variables (keyword)
    """
    with open(filename,'wb') as pickle_file:
        pickler = __pickle.Pickler(pickle_file)
        pickler.dump([argv, kwargs])
        return 0
    return 1 # pragma: no cover

def LoadData(filename):
    """Load data from disk using pickle format

       Load data saved with the 'saveData' from file
       return an instance of IOHelper if ok
       return None if not ok
    """
    with open(filename,'rb') as pickle_file:
        unpickler = __pickle.Unpickler(pickle_file)
        data = unpickler.load()
        return  IOHelper(data)
    return None # pragma: no cover

def CheckIntegrity():
    """ AutoTest routine """

    from  BasicTools.Helpers.Tests import TestTempDir
    # create a temp file
    tempdir = TestTempDir.GetTempPath()
    try :
        # Save data
        SaveData(tempdir + "testFile.data","two", 3, (3,5),toto=10)
        # load data
        b = LoadData(tempdir + "testFile.data")
        # test correct data
        if(b.unamed[0] != "two"): raise Exception()
        if(b.unamed[1] != 3): raise Exception()
        if(b.unamed[2] != (3,5)): raise Exception()
        if(b.named['toto'] != 10): raise Exception()
        output = b.__str__()
        print(b)
        # delete temp directory
        return 'Ok'
    except:# pragma: no cover
        # delete temp directory
        raise

if __name__ == '__main__':
    #import time
    #stime = time.time()
    print(CheckIntegrity()) # pragma: no cover

    #print(time.time()-stime)
