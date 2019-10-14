# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       
import pickle as pickle
import socket
import signal


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
    for x in range(level, 1, -1):
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
    def __init__(self,timeout=3600):
        super(WormholeBase,self).__init__()
        self.socket = None
        self.otherSideR =  None
        self.otherSideS =  None
        self.proto = 0
        self.timeout = timeout

        # solution for time out form https://stackoverflow.com/questions/492519/timeout-on-a-function-call
        def TimeOutHandler(signum, frame):
           print("Connection TimeOut")
           exit(1)
        if not (timeout is None):
            signal.signal(signal.SIGALRM, TimeOutHandler)

    def Receive(self):
      if not (self.timeout is None):
          signal.alarm(self.timeout)

      try:
          sizestream = ""
          while (len(sizestream) < 64):
              if self.socket is None:
                  sizestream += self.otherSideR.read(1).decode('utf8')
              else:
                  sizestream += self.otherSideR.recv(1).decode('utf8')
          size = int(sizestream)

          if self.socket is None:
              datastream = self.otherSideR.read(size)
          else:
              datastream = self.otherSideR.recv(size)
          ldata = len(datastream )
          while ldata < size:
              if self.socket is None:
                  datastream += self.otherSideR.read(size-ldata)
              else:
                  datastream += self.otherSideR.recv(size-ldata)
              ldata = len(datastream )
          if int(sys.version_info.major) >= 3:
              data = pickle.loads(datastream,encoding = 'latin1')
          else:
              data = pickle.loads(datastream,)
          return data
      except Exception:
          exit(1)


    def Send(self,data):
        print("Sending data")
        if int(sys.version_info.major) >= 3:
            streamdata = pickle.dumps(data,self.proto,fix_imports=True)
        else:
            streamdata = pickle.dumps(data,self.proto)

        data = str(len(streamdata)).zfill(64)
        if self.socket is None:
            self.otherSideS.write( data.encode('utf8'))
            self.otherSideS.write(streamdata)
            self.otherSideS.flush()
        else:
            self.otherSideS.send( data.encode('utf8'))
            self.otherSideS.send(streamdata)

    def Close(self):

        if not (self.otherSideR is self.otherSideS):
            self.otherSideS.close()
        self.otherSideR.close()

class WormholeServer(BaseOutputObject):
    def __init__(self,port = None, cmd=None ,dry=False,timeout=3600):
       super(WormholeServer,self).__init__()
       self.globals = {}

       # no code is executes only print it
       self.drymode = dry
       if port is not None:
           self.communicator = WormholeBase(timeout=timeout)
           self.ListenUsingPort(port)
           self.MainLoop()
           self.communicator.socket.close()
       elif cmd is not None:
           #from BasicTools.IO.Proxy import ServerProxy
           self.communicator = WormholeBase(timeout=timeout)
           from BasicTools.Helpers.PrintBypass import PrintBypass
           self.printBypass = PrintBypass()
           from BasicTools.Helpers.Tests import GetUniqueTempFile
           out = GetUniqueTempFile(".log","output_")[1]
           err = GetUniqueTempFile(".log2","output_")[1]
           self.printBypass.ToDisk(out,err)
           #self.printBypass.ToSink()
           self.StartUsingPipe()
           self.MainLoop()
           self.printBypass.Restore()

    def StartUsingPipe(self):
        if int(sys.version_info.major) >= 3:
            self.communicator.otherSideR = self.printBypass.stdin_.buffer
            self.communicator.otherSideS = self.printBypass.stdout_.buffer
        else:
            self.communicator.otherSideR = self.printBypass.stdin_
            self.communicator.otherSideS = self.printBypass.stdout_

    def ListenUsingPort(self,port=None):
        if port is None:
              port = 12345

        self.communicator.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.communicator.socket.bind(('', port))
        self.communicator.socket.listen(0)
        self.communicator.otherSideR, address = self.communicator.socket.accept()
        self.communicator.otherSideS = self.communicator.otherSideR
        print("(s) {0} connected".format( address ))


    def ProtocolNegotiation(self):
        ClientHighestProtocol = self.communicator.Receive()
        proto = min(pickle.HIGHEST_PROTOCOL, ClientHighestProtocol)
        print("(s) Using Protocol " + str(proto))
        self.communicator.Send(proto)
        self.communicator.proto = proto


    def MainLoop(self):

        while(True):
            action = self.communicator.Receive()
            #if len(rdata) == 0 :
            #    continue
            #action = rdata.decode()

            self.PrintDebug(action)
            if action == "p":
                self.ProtocolNegotiation()

            elif action == "s":
                #print("receiving data")
                key = self.communicator.Receive()
                #print("receiving data for ",key)
                value = self.communicator.Receive()
                print(key ," = ", value)
                if self.drymode:
                    print("(s)" + str(key) +" = " + str(value))
                else:
                    self.globals[key] = value
                #print(self.globals)
                #eval(key + " = value",self.globals)

            elif action == "c":
                expression = self.communicator.Receive()
                print(expression)
                if self.drymode:
                    print("(s) Exec: " + expression )
                else:
                    exec(str(expression), self.globals)

            elif action == "r":
                print("(s) Sending data back : "+ str(key))
                if self.drymode:
                    self.Send(None)
                else:
                    key = self.communicator.Receive()
                    self.communicator.Send(self.globals[key])

            elif action  == "x":
                print("(s) exit")
                self.communicator.Close()
                return

            else:
                print("Dont know how to treat " + str(action))
                return

class WormholeClient(BaseOutputObject):
    def __init__(self,port = None,host=None,proc=None):
       super(WormholeClient,self).__init__()
       self.communicator = None

       if port is not None:
           self.communicator = WormholeBase()
           self.Connect(port=port, host=host)
       elif proc is not None:
           #from BasicTools.IO.Proxy import ServerProxy
           self.communicator = WormholeBase()
           self.StartUsingPipe(proc)

    def StartUsingPipe(self,proc):
        self.communicator.otherSideS = proc.stdin
        self.communicator.otherSideR = proc.stdout
        self.ProtocolNegotiation()
        self.communicator.proto = self.proto

    def Connect(self,port=None, host=None):
        if port is None:
            port = 12345

        if host is None:
            host = "localhost"

        self.communicator.otherSideS = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.communicator.socket = True

        self.communicator.otherSideR  = self.communicator.otherSideS
        print("(c) connecting to "),
        print((host, port))
        self.communicator.otherSideS.connect((host, port))
        print("(c) Connection on {0}".format(port))
        self.ProtocolNegotiation()
        print("(c) Using protocol " + str(self.proto) )
        self.communicator.proto = self.proto

    def ProtocolNegotiation(self):
        #request protocol negotiation
        print("pickle.HIGHEST_PROTOCOL ",pickle.HIGHEST_PROTOCOL )
        self.communicator.Send("p")
        self.communicator.Send(pickle.HIGHEST_PROTOCOL )
        #self.Communicator.send(str(pickle.HIGHEST_PROTOCOL)[0].encode() )
        self.proto = self.communicator.Receive()
        print("self.proto ",self.proto )

    def SendData(self,key,data):

        print("(c) sending " + str(key) + " :: " + str(data))
        self.communicator.Send("s")
        self.communicator.Send(key)
        self.communicator.Send(data)

    def RemoteExec(self,expression):
        print("Remote Exec : '" + str(expression) +"'")
        self.communicator.Send("c")
        self.communicator.Send(expression)


    def RetrieveData(self,variable):
        self.communicator.Send("r")
        self.communicator.Send(variable)
        return self.communicator.Receive()

    def Exit(self):
        print("Sending Exit")
        self.communicator.Send("x")
        self.communicator.Close()


def CheckIntegrityNetWork():

   import time
   testport = GetAnFreePortNumber()

   import threading
   #try:
   if True:
     def runServer():
         print("(s) Starting Server Side ",testport)
         WormholeServer(testport,dry=False,timeout=None)
     TT = threading.Thread(target=runServer )

     TT.start()
     time.sleep(3.)
     print("(c) Starting Client Side ",testport)
     client = WormholeClient(testport)
     client.SendData("Hola",5)
     client.RemoteExec("Hola += 3")
     newhola = client.RetrieveData("Hola")
     client.Exit()
     print("Done")
     if newhola == 8:
         return 'ok'
     return "Not ok"
   try:
      pass
   except:
     TT.join(0)
     return "Not  OK"

def GetPipeWrormholeScript():
    from BasicTools.Helpers.Tests import TestTempDir

    return """
from BasicTools.IO.Wormhole import WormholeServer
from BasicTools.Helpers.Tests import TestTempDir
TestTempDir.SetTempPath("{0}")
a = WormholeServer(cmd="")
exit()
""".format(TestTempDir.GetTempPath())

def CheckIntegrityPipe():

   import time

   #try:
   if True:
     def runServerPipe(cmd):

         script = GetPipeWrormholeScript()
         print("(s) Starting Server Side")
         import subprocess, os
         proc = subprocess.Popen([cmd, '-c', script], cwd=os.getcwd(),
                 stdout=subprocess.PIPE, stdin=subprocess.PIPE)
         return proc

     proc = runServerPipe(sys.executable)
     print("(c) Starting Client Side")
     time.sleep(1.)

     client = WormholeClient(proc=proc)
     client.SendData("Hola",5)
     client.RemoteExec("Hola += 3")
     newhola = client.RetrieveData("Hola")
     client.Exit()
     print("Done")
     if newhola == 8:
         return 'ok'
     return "Not ok"
   try:
       pass
   except:
     TT.join(0)
     return "Not  OK"


def GetAnFreePortNumber():
    import socket
    s=socket.socket()
    s.bind(("", 0))
    portNumber = s.getsockname()[1]
    s.close()
    return portNumber

def CheckIntegrity(GU=False):
    res = CheckIntegrityNetWork()
    if str(res).lower() != "ok":
        return res
    return CheckIntegrityPipe()

if __name__ == '__main__':

  def RunClient(testport ):
    print("Client Side")
    client = WormholeClient(testport)

    client.SendData("Hola",5)
    client.RemoteExec("Hola += 3")
    newhola = client.RetrieveData("Hola")
    print("new Hola :",newhola )
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
    t = threading.Thread(target=WormholeServer,name="WormHoleServer",kwargs={"port":testport,"timeout":None})
    t.daemon = True
    t.start()
    import time
    time.sleep(1)

    RunClient(testport)
    t.join(5)
