#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-
import sys

from BasicTools.IO.IOFactory import RegisterWriterClass,RegisterReaderClass

class PipeReader(object):
    def __init__(self):
       super(PipeReader,self).__init__()

    def Open(self):
       self.inbuffer = sys.stdin.buffer

    def SetFileName(self,val):
        pass

    def SetBinary(self,val):
        pass

    def Read(self):
        import pickle
        return pickle.load(self.inbuffer,encoding = 'latin1')

##    def ReadStr(self):
##        import pickle
 #       return pickle.loads(self.inbuffer,encoding = 'latin1')


    def Open(self):
        pass

    def Close(self):
        pass


class PipeWriter(object):
    def __init__(self):
       super(PipeWriter,self).__init__()
       self.outbuffer = None

    def SetFileName(self,val):
        pass

    def SetBinary(self,val):
        pass

    def Write(self,outmesh,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None):
        import pickle
        pickle.dump(outmesh,self.outbuffer,0)
        self.outbuffer.flush()

#    def WriteStr(self,outmesh,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None):
#        import pickle
#        self.outbuffer.write(pickle.dumps(outmesh).decode('latin1'))

    def Open(self):
        import sys
        if int(sys.version_info.major) >= 3:
            self.outbuffer = sys.stdout.buffer
        else:
            self.outbuffer = sys.stdout

    def Close(self):
        pass

RegisterReaderClass(".PIPE",PipeReader)
RegisterWriterClass(".PIPE",PipeWriter)

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    mesh = CreateCube()
    print(mesh)
    PipeReader()
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
