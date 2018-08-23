# -*- coding: utf-8 -*-
import pickle as pickle
import socket

__author__ = "Felipe Bordeu"
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

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
          while (len(sizestream) < 64):
              sizestream += self.otherSide.recv(1).decode()
          #print(sizestream)
          size = int(sizestream)
          #if self.socket is not None:
          #    print("s : " + str(size))
          #else:
          #    print("c : " + str(size))

          datastream = self.otherSide.recv(size)
          ldata = len(datastream )
          while ldata < size:
              datastream += self.otherSide.recv(size-ldata)
              ldata = len(datastream )

          data = pickle.loads(datastream)
          return data

    def _internalSendData(self,data):
        streamdata = pickle.dumps(data,self.proto)
        data = str(len(streamdata)).zfill(64)
        self.otherSide.send( data.encode())
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

    def ProtocolNegotiation(self):
        ClientHighestProtocol = int(self.otherSide.recv(1).decode())
        self.proto = min(pickle.HIGHEST_PROTOCOL, ClientHighestProtocol)
        print("(s) Using Protocol " + str(self.proto))
        self.otherSide.send(str(self.proto).encode() )


    def MainLoop(self):


        self.socket.listen(0)
        self.otherSide, address = self.socket.accept()
        print("(s) {0} connected".format( address ))





        while(True):
            rdata = self.otherSide.recv(1)
            if len(rdata) == 0 :
                continue
            action = rdata.decode()

            self.PrintDebug(action)
            if action == "p":
                self.ProtocolNegotiation()
            elif action == "s":
                #print("receiving data")
                key = self._internalReciveData()
                value = self._internalReciveData()
                #print(key + " = value")
                if self.drymode:
                    print("(s)" + str(key) +" = " + str(value))
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
                print(expression)
                if self.drymode:
                    print("(s) Exec: " + expression )
                else:
                    exec(str(expression), self.globals)

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

                print("(s) Sending data back : "+ str(key))

                if self.drymode:
                    self._internalSendData(None)
                else:
                    key = self._internalReciveData()
                    self._internalSendData(self.globals[key])

            elif action  == "x":
                print("(s) exit")
                self.otherSide.close()
                return
            else:
                print("Dont know how to treat " + str(action))
                return


class WormholeClient(WormholeBase):
    def __init__(self):
       super(WormholeClient,self).__init__()
       self.proto = 2

    def Connect(self,port=None, host=None):
        if port is None:
            port = 12345

        if host is None:
            host = "localhost"

        self.otherSide = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        print("(c) connecting to "),
        print((host, port))
        self.otherSide.connect((host, port))
        print("(c) Connection on {0}".format(port))
        self.ProtocolNegotiation()
        print("(c) Using protocol " + str(self.proto) )

    def ProtocolNegotiation(self):
        #request protocol negotiation
        self.otherSide.send(b"p")
        self.otherSide.send(str(pickle.HIGHEST_PROTOCOL)[0].encode() )
        self.proto = int(self.otherSide.recv(1).decode())


    def SendData(self,key,data):

        print("(c) sending " + str(key) + " :: " + str(data))
        self.otherSide.send(b"s")
        self._internalSendData(key)
        self._internalSendData(data)

    #def RemoteEval(self,key, expression):
    #    self.otherSide.send("e")
    #    self._internalSendData(key)
    #    self._internalSendData(expression)

    def RemoteExec(self,expression):
        self.otherSide.send(b"c")
        self._internalSendData(expression)


    def RetrieveData(self,variable):
        self.otherSide.send(b"r")
        self._internalSendData(variable)
        return self._internalReciveData()

    #def ImportModule(self,module,asname):
    #    self.otherSide.send("i")
    #    self._internalSendData(module)
    #    self._internalSendData(asname)
    #    print(self._internalReciveData())

    def Exit(self):
        self.otherSide.send(b"x")
        self.otherSide.close()


def CheckIntegrity():


   testport = 12349

   import threading
   try:

     def runServer():
         print("(s) Starting Server Side")
         WormholeServer(testport,dry=False)
     TT = threading.Thread(target=runServer )
     TT.start()


     client = WormholeClient()
     client.Connect(testport)
     client.SendData("Hola",5)
     client.RemoteExec("Hola += 3")
     newhola = client.RetrieveData("Hola")
     client.Exit()

     if newhola == 8:
         return 'ok'
   except:
     TT.join(0)
     return "Not OK"


def GetAnFreePortNumber():
    import socket;
    s=socket.socket();
    s.bind(("", 0));
    portNumber = s.getsockname()[1];
    s.close()
    return portNumber


if __name__ == '__main__':

  def RunClient(testport ):
    print("Client Side")
    client = WormholeClient()
    client.Connect(testport)
    client.SendData("Hola",5)
    client.RemoteExec("Hola += 3")
    newhola = client.RetrieveData("Hola")
    #client.ImportModule("BasicTools.IO.XdmfWriter","XdmfWriter")
    client.RemoteExec("from BasicTools.IO import XdmfWriter")
    client.RemoteExec("from BasicTools.IO.XdmfWriter import WriteMeshToXdmf")
    import numpy as np
    client.SendData("p",np.array([5, 10, 4]))
    #client.RemoteEval("v","XdmfWriter.ArrayToString(p)")
    client.RemoteExec("v =  XdmfWriter.ArrayToString(p)")
    client.RemoteExec("a =2 ")
    v = client.RetrieveData("v")
    print(type(v))
    print("v " + str(v))
    p = client.RetrieveData("p")
    print(type(p))
    print("p " + str(p))
    client.Exit()

  testport = 12348
  if len(sys.argv) > 1 and  sys.argv[1]  == "-s":
     print("Server Side")
     WormholeServer(testport,dry=False)
  elif len(sys.argv) > 1 and  sys.argv[1]  == "-c":
     RunClient(testport)
  else:
    testport = GetAnFreePortNumber()

    import threading
    t = threading.Thread(target=WormholeServer,name="WormHoleServer",args=(testport,))
    t.daemon = True
    t.start()
    import time

    time.sleep(0.1)

    RunClient(testport)


    t.join(5)
