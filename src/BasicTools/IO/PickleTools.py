# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import pickle as __pickle

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
        unpickler = __pickle.Unpickler(pickle_file,encoding='latin1')
        data = unpickler.load()
        return  IOHelper(data)
    return None # pragma: no cover

from BasicTools.IO.ReaderBase import ReaderBase

class PickleReader(ReaderBase):
    def __init__(self,fileName = None):
        super(PickleReader,self).__init__(fileName=fileName)
        self.canHandleTemporal = False

    def Read(self):
        internalReader = LoadData(self.fileName)
        self.output = internalReader.named["mesh"]
        return self.output

class PickleWriter(object):
    def __init__(self):
        super(PickleWriter,self).__init__()
        self.filename = ""
        self.canHandleBinaryChange = False

    def SetBinary(self,val=True):
        pass

    def SetFileName(self,filename):
        self.filename = filename;

    def Open(self,fileName=None):
        if not fileName is None:
            self.SetFileName(fileName)

    def Close(self):
        pass


    def Write(self,mesh, PointFields = None, CellFields = None, GridFields= None, PointFieldsNames = None, CellFieldsNames= None, GridFieldsNames=None):

        if PointFieldsNames is not None:
            nodeFields = {k:v for k,v in zip(PointFieldsNames,PointFields)}
        else:
            nodeFields = {}

        if CellFieldsNames is not None:
            elemFields = {k:v for k,v in zip(CellFieldsNames,CellFields)}
        else:
            elemFields = {}
        import copy
        cmesh = copy.copy(mesh)
        cmesh.nodeFields = nodeFields
        cmesh.elemFields = elemFields
        SaveData(self.filename,mesh= cmesh)


from BasicTools.IO.IOFactory import RegisterWriterClass, RegisterReaderClass
RegisterReaderClass(".pickle",PickleReader)
RegisterWriterClass(".pickle",PickleWriter)


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
    except:# pragma: no cover
        # delete temp directory
        raise

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateUniformMeshOfBars
    barmesh = CreateUniformMeshOfBars(0,8,10)
    import numpy as np
    PointFields = [np.arange(barmesh.GetNumberOfNodes())]
    PointFieldsNames = ["PointData"]

    CellFields = [np.arange(barmesh.GetNumberOfElements())]
    CellFieldsNames = ["PointData"]

    print(barmesh)
    pw = PickleWriter()
    pw.SetBinary()# this has no effect
    pw.Open(tempdir + "testFile.pickle")
    pw.Write(barmesh)
    pw.Write(barmesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFields=CellFields,CellFieldsNames=CellFieldsNames)
    pw.Close() # this has no effect
    pr = PickleReader()
    pr.SetFileName(tempdir + "testFile.pickle")
    barmeshII = pr.Read()

    from BasicTools.Containers.MeshTools import IsClose
    IsClose(barmesh,barmeshII)
    print(barmeshII)

    return 'Ok'

if __name__ == '__main__':
    #import time
    #stime = time.time()
    print(CheckIntegrity()) # pragma: no cover

    #print(time.time()-stime)
