# -*- coding: utf-8 -*-

from BasicTools.IO.UniversalReader import ReadMesh
from BasicTools.IO.UniversalWriter import WriteMesh

def PrintHelp():
  print( 'python  MeshFileConverter -i <inputfile> -o <outputfile>')
  print( 'options :')
  print( '       -i    Input file name')
  print( '       -o    output file name')
  print( '       -h    this help')
  sys.exit(2)

if __name__ == '__main__' :
    import sys, getopt
    if len(sys.argv) == 1:
        PrintHelp()
        sys.exit()
    else:
      try:
          opts, args = getopt.getopt(sys.argv[1:],"hi:o:")
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

      print("Start Reading...")
      mesh = ReadMesh(inputfilename)
      print("Start Writing...")
      WriteMesh(outputfilename,mesh)
      print("DONE")


def CheckIntegrity(GUI=False):
    return "ok"
