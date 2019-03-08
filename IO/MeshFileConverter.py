# -*- coding: utf-8 -*-

from BasicTools.IO.UniversalReader import ReadMesh
from BasicTools.IO.UniversalWriter import WriteMesh
import  BasicTools.IO.IOFactory as IOF

def PrintHelp():
  print( 'python  MeshFileConverter -i <inputfile> -o <outputfile>')
  print( 'options :')
  print( '       -i    Input file name')
  print( '       -o    output file name')
  print( '       -h    this help')
  IOF.InitAllReaders()
  print("Available Readers : ", IOF.GetAvailableReaders())

  IOF.InitAllWriters()
  print("Available Writers : ", IOF.GetAvailableWriter())
  sys.exit(2)

def Convert(inputfilename,outputfilename):
      IOF.InitAllReaders()
      IOF.InitAllWriters()

      print("Start Reading...")
      mesh = ReadMesh(inputfilename)
      print(mesh)
      print("Start Writing to ", outputfilename)
      WriteMesh(outputfilename,mesh)
      print("DONE")


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
            Convert(inputfilename,outputfilename)

    return "ok"

if __name__ == '__main__' :
    import sys, getopt
    if len(sys.argv) == 1:
        PrintHelp()
        sys.exit()
    else:
      try:
          opts, args = getopt.getopt(sys.argv[1:],"thi:o:")
      except getopt.GetoptError:
          PrintHelp()
          sys.exit(2)

      for opt, arg in opts:
         if opt == '-h':
             PrintHelp()
             sys.exit()
         elif opt in ("-i"):
            inputfilename = arg
         elif opt in ("-o"):
            outputfilename = arg
         elif opt in ("-t"):
            print(CheckIntegrity(GUI=True))
            exit(0)


    Convert(inputfilename,outputfilename)