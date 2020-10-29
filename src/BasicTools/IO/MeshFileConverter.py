# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


from BasicTools.Helpers.TextFormatHelper import TFormat as TF
from BasicTools.IO.UniversalReader import ReadMesh
from BasicTools.IO.UniversalWriter import WriteMesh
import  BasicTools.IO.IOFactory as IOF

IOF.InitAllReaders()
IOF.InitAllWriters()


def PrintHelp():
  print( 'python  MeshFileConverter -i <inputfile> -o <outputfile>')
  print( 'options :')
  print( '       -i    Input file name')
  print( '       -o    output file name')
  print( '       -h    this help')
  print( '       -t    time to read (if the input file can handle time) (default last time step is writen)')
  print( '       -p    print availables time to read ')
  print( '       -v    more verbose output ')
  print( '       -m    Activate MeshIO Readers and Writers ')
  print( '       -s    Plot mesh before continue (press key "q" exit)')

  print("Available Readers : ", IOF.GetAvailableReaders())
  print("Available Writers : ", IOF.GetAvailableWriter())

#MeshFileConverter -i meshfile.meshb -o .PIPE > toto

from BasicTools.Helpers.PrintBypass import PrintBypass

def Convert(inputfilename,outputfilename,ops):
      IOF.InitAllReaders()
      IOF.InitAllWriters()

      with PrintBypass() as f:

          if ".PIPE" in outputfilename :
              f.ToDisk("MeshFileConverter.log")

          print("Start Reading...", inputfilename)

          if ops["printTimes"]:
              from BasicTools.IO.IOFactory import InitAllReaders
              InitAllReaders()
              import os
              basename,extention = os.path.splitext(os.path.basename(inputfilename))
              from BasicTools.IO.IOFactory import CreateReader
              reader = CreateReader("."+inputfilename.split(".")[-1].lower())
              reader.SetFileName(inputfilename)
              if reader.canHandleTemporal :
                  reader.ReadMetaData()
                  print("Available Times in files:")
                  print(reader.GetAvilableTimes())
                  import sys
                  sys.exit(0)

          mesh = ReadMesh(inputfilename,timeToRead = ops["timeToRead"])
          if ops["PlotOnScreen"]:
              from BasicTools.Containers.vtkBridge import PlotMesh
              PlotMesh(mesh)
          if len(outputfilename) == 0:
              PrintMeshInformation(mesh)
              print("No output file name")
              print("Done")
              return
          else:
              print(mesh)

          print("Start Writing to "+  str(outputfilename))
          writer = None


          from BasicTools.IO.IOFactory import CreateWriter
          if ".PIPE" in outputfilename :
              writer = CreateWriter("."+outputfilename.split(".")[-1])
              writer.outbuffer = f.stdout_.buffer

          WriteMesh(outputfilename,mesh,writer=writer)
          print("DONE")

def PrintMeshInformation(mesh):

      def L25(text):
          return TF.Left(text ,fill=" ",width=25)
      def L15(text):
          return TF.Left(text ,fill=" ",width=15)

      print(TF.Center("check the integrity of the mesh"))
      print(L25("nodes.shape:") , str(mesh.nodes.shape) )
      mesh.ComputeBoundingBox()
      print(L25("boundingMin:") + str(mesh.boundingMin) )
      print(L25("boundingMax:") + str(mesh.boundingMax) )
      print(L25("originalIDNodes.shape:") , str(mesh.originalIDNodes.shape) )
      print(L25("min(originalIDNodes):") , str(min(mesh.originalIDNodes) )  )
      print(L25("max(originalIDNodes):") , str(max(mesh.originalIDNodes) )  )
      print("---- node tags ---- " )
      print(mesh.nodesTags)
      for tag in mesh.nodesTags:
          res = L15(tag.name) + L15(" size:"+ str(len(tag)))
          if len(tag):
              res +=L15(" min:"+str(min(tag.GetIds()) ) )+ "  max:"+str(max(tag.GetIds()) )
          print(res)


      print("---- Elements ---- " )
      for name,data in mesh.elements.items():
          res = L25(str(type(data)).split("'")[1].split(".")[-1]+ " " + name+" ")
          res += L15("size:" + str(data.GetNumberOfElements()))
          res += L25("min(connectivity):" + str(min(data.connectivity.ravel())))
          res += L25("max(connectivity):" +  str(max(data.connectivity.ravel())))
          print(res)
          #if len(data.tags):
          #    print("  ---- Element tags for "+name+" ---- " )
          for tag in data.tags:
              res = L15(tag.name) + L15(" size:"+ str(len(tag)))
              if len(tag):
                  res +=L15(" min:"+str(min(tag.GetIds()) ) )+ "  max:"+str(max(tag.GetIds()) )
              print("  Tag: " +res)

      print(TF.Center("check the integrity of the mesh DONE"))

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.TestData import GetTestDataPath

    inputfiles = ["coneAscii.stl",
                  "coneBinary.stl",
                  "cube.geof",
                  "GCodeTest.gcode",
                  "mesh1.msh",
                  #"Structured.xmf",
                  #"Unstructured.xmf"
            ]

    outputext = [ "geof",
                  "mesh",
                  "msh",
                  "stl",
                  "xdmf"
            ]

    for iff in inputfiles:
        for off in outputext:
            inputfilename = GetTestDataPath() + iff
            outputfilename = TestTempDir().GetTempPath()+iff+"." + off
            ops= {}
            ops["timeToRead"] = -1
            ops["printTimes"] = False
            ops["PlotOnScreen"] = GUI
            Convert(inputfilename,outputfilename,ops)

    return "ok"

def Main():
    import sys, getopt
    if len(sys.argv) == 1:
        PrintHelp()
        sys.exit()
    else:
      #try:
      if True:
          opts, args = getopt.getopt(sys.argv[1:],"svphmt:i:o:")
      #except getopt.GetoptError:
      #    PrintHelp()
      #    sys.exit(2)

      outputfilename = ""
      ops= {}
      ops["timeToRead"] = -1
      ops["printTimes"] = False
      ops["PlotOnScreen"] = False
      for opt, arg in opts:
         if opt == '-h':
             PrintHelp()
             sys.exit()
         elif opt in ("-i"):
            inputfilename = arg
         elif opt in ("-o"):
            outputfilename = arg
         elif opt in ("-t"):
            ops["timeToRead"] = float(arg)
            print(ops)
            #raise
         elif opt in ("-p"):
            ops["printTimes"] = bool(True)
         elif opt in ("-v"):
            from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
            BOO.SetGlobalDebugMode(True)
         elif opt in ("-c"):
            print(CheckIntegrity(GUI=True))
            sys.exit(0)
         elif opt in ("-m"):
             from BasicTools.Containers.MeshIOBridge import InitAllReaders,InitAllWriters,AddReadersToBasicToolsFactory,AddWritersToBasicToolsFactory
             InitAllReaders()
             InitAllWriters()
             AddReadersToBasicToolsFactory()
             AddWritersToBasicToolsFactory()
         elif opt in ("-s"):
             ops["PlotOnScreen"] = True


    Convert(inputfilename,outputfilename,ops )

if __name__ == '__main__' :
    Main()
