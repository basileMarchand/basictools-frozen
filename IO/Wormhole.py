# -*- coding: utf-8 -*-
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Helpers.TextFormatHelper import TFormat as TFormat

import pickle as pickle
#import cPickle as pickle

import socket


"""Backport of importlib.import_module from 3.x."""
# While not critical (and in no way guaranteed!), it would be nice to keep this
# code compatible with Python 2.3.
#code to make work the import_module in abaqus python 2.6
# if python of abaqus gets updated please erase this function
import sys

def _resolve_name(name, package, level):
    """Return the absolute name of the module to be imported."""
    if not hasattr(package, 'rindex'):
        raise ValueError("'package' not set to a string")
    dot = len(package)
    for x in xrange(level, 1, -1):
        try:
            dot = package.rindex('.', 0, dot)
        except ValueError:
            raise ValueError("attempted relative import beyond top-level "
                              "package")
    return "%s.%s" % (package[:dot], name)


def import_module(name, package=None):
    """Import a module.

    The 'package' argument is required when performing a relative import. It
    specifies the package to use as the anchor point from which to resolve the
    relative import to an absolute import.

    """
    if name.startswith('.'):
        if not package:
            raise TypeError("relative imports require the 'package' argument")
        level = 0
        for character in name:
            if character != '.':
                break
            level += 1
        name = _resolve_name(name[level:], package, level)
    __import__(name)
    return sys.modules[name]
####################################




class WormholeBase(BaseOutputObject):
    def __init__(self,port = None):
        super(WormholeBase,self).__init__()
        self.socket = None
        self.otherSide =  None

    def _internalReciveData(self):
          sizestream = ""
          while (len(sizestream) < 8):
              sizestream += self.otherSide.recv(1)
          #print(sizestream)
          size = int(sizestream)
          #print(size)
          datastream = self.otherSide.recv(size)
          data = pickle.loads(datastream)
          return data

    def _internalSendData(self,data):
        streamdata = pickle.dumps(data)
        self.otherSide.send( str(len(streamdata)).zfill(8)  )
        self.otherSide.send(streamdata)


class WormholeServer(WormholeBase):
    def __init__(self,port = None,dry=False):
       super(WormholeServer,self).__init__()
       self.globals = {}

       # no code is executes only print it
       self.drymode = dry
       if port is not None:
           self.Listen(port)


    def Listen(self,port=None):
        if port is None:
              port = 12345

        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.bind(('', port))

        self.MainLoop()
        self.socket.close()



    def MainLoop(self):


        self.socket.listen(0)
        self.otherSide, address = self.socket.accept()
        print "{0} connected".format( address )
        while(True):
            action = self.otherSide.recv(1)
            if action == "s":
                #print("receiving data")
                key = self._internalReciveData()
                value = self._internalReciveData()
                #print(key + " = value")
                if self.drymode:
                    print(str(key) +" = " + str(value))
                else:
                    self.globals[key] = value
                #print(self.globals)
                #eval(key + " = value",self.globals)

#            elif action == "e":
#                print("evaluation expression")
#                key = self._internalReciveData()
#                expretion = self._internalReciveData()
#                self.globals[key] =  eval(expretion,self.globals)
#                print(self.globals)
            elif action == "c":
                expression = self._internalReciveData()
                if self.drymode:
                    print("Exec: " + expression )
                else:
                    exec expression in self.globals

                #print(self.globals)
#            elif action == "i":
#                print("import module")
#                module = self._internalReciveData()
#                asname = self._internalReciveData()
#                print(module)
#                print(asname)
#
#                #try:
#                if True:
#                    #import importlib
#                    #self.globals[asname] =  importlib.import_module(module)
#                    #import importlib
#                    self.globals[asname] =  import_module(module)
#                    print("-------------------------------------")
#                    print(self.globals[asname])
#
#                    #print(self.globals)
#                    self._internalSendData("OK")
#                #except :
#                #    self._internalSendData("KO")
#

            elif action == "r":
                print("Sending data back")

                if self.drymode:
                    print("Sending back: = " + str(key))
                    self._internalSendData(None)
                else:
                    key = self._internalReciveData()
                    self._internalSendData(self.globals[key])

            elif action  == "x":
                print("exit")
                self.otherSide.close()
                return
            else:
                print("Dont know how to treat " + str(action))
                return


class WormholeClient(WormholeBase):
    def __init__(self):
       super(WormholeClient,self).__init__()

    def Connect(self,port=None, host=None):
        if port is None:
            port = 12345

        if host is None:
            host = "localhost"

        self.otherSide = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        print("connecting to "),
        print((host, port))
        self.otherSide.connect((host, port))
        print "Connection on {0}".format(port)

    def SendData(self,key,data):

        self.otherSide.send("s")
        self._internalSendData(key)
        self._internalSendData(data)

    #def RemoteEval(self,key, expression):
    #    self.otherSide.send("e")
    #    self._internalSendData(key)
    #    self._internalSendData(expression)

    def RemoteExec(self,expression):
        self.otherSide.send("c")
        self._internalSendData(expression)


    def RetrieveData(self,variable):
        self.otherSide.send("r")
        self._internalSendData(variable)
        return self._internalReciveData()

    #def ImportModule(self,module,asname):
    #    self.otherSide.send("i")
    #    self._internalSendData(module)
    #    self._internalSendData(asname)
    #    print(self._internalReciveData())

    def Exit(self):
        self.otherSide.send("x")
        self.otherSide.close()


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1 and  sys.argv[1]  == "-s":
        print("Server Side")
        WormholeServer(12345)
    else:
        print("Client Side")
        client = WormholeClient()
        client.Connect(12345)
        client.SendData("Hola",5)
        client.RemoteExec("Hola += 3")
        newhola = client.RetrieveData("Hola")
        #client.ImportModule("BasicTools.IO.XdmfWriter","XdmfWriter")
        client.RemoteExec("from BasicTools.IO import XdmfWriter")
        client.RemoteExec("from BasicTools.IO.XdmfWriter import WriteMeshToXdmf")
        client.SendData("p",[5, 10, 4])
        #client.RemoteEval("v","XdmfWriter.ArrayToString(p)")
        client.RemoteExec("v =  XdmfWriter.ArrayToString(p)")
        client.RemoteExec("a =2 ")
        v = client.RetrieveData("v")
        print("v " + str(v))
        client.Exit()

