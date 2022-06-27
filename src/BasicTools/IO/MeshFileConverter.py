# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


from BasicTools.IO.UniversalReader import ReadMesh
from BasicTools.IO.UniversalWriter import WriteMesh
import  BasicTools.IO.IOFactory as IOF

def LoadReadersAndWriters(ops = None):
    if ops is not None and ops.get("OnlyAbaqusReader",False):
        from BasicTools.IO.OdbReader import OdbReader
        import BasicTools.IO.PickleTools
        import BasicTools.IO.XdmfWriter
    else:
        IOF.InitAllReaders()
        IOF.InitAllWriters()

    if ops is not None and ops.get("MeshIO",False):
        from BasicTools.Bridges.MeshIOBridge import InitAllReaders,InitAllWriters,AddReadersToBasicToolsFactory,AddWritersToBasicToolsFactory
        InitAllReaders()
        InitAllWriters()
        AddReadersToBasicToolsFactory()
        AddWritersToBasicToolsFactory()

def PrintHelp(ops):
    LoadReadersAndWriters(ops)

    print( 'python  MeshFileConverter -i <inputfile> -o <outputfile>')
    print( 'options :')
    print( '       -i    Input file name')
    print( '       -o    output file name')
    print( '       -h    this help')
    print( '       -t    time to read (if the input file can handle time) (default last time step is writen)')
    print( '       -p    print availables time to read ')
    print( '       -T    Convert all the time steps')
    print( '       -v    more verbose output ')
    print( '       -m    Activate MeshIO Readers and Writers ')
    print( '       -s    Plot mesh before continue (press key "q" exit)')
    print( '       -a    Abaqus Mode (python2). Load only the abaqus reader and pickle writer')
    print( '       -b    Force binary output on the writer')
    print( '       -c    (reserved)')

    print("Available Readers : ", IOF.GetAvailableReaders())
    print("Available Writers : ", IOF.GetAvailableWriter())

#MeshFileConverter -i meshfile.meshb -o .PIPE > toto

from BasicTools.Helpers.PrintBypass import PrintBypass

def Convert(inputfilename,outputfilename,ops):

    LoadReadersAndWriters(ops)

    with PrintBypass() as f:

          if ".PIPE" in outputfilename :
              f.ToDisk("MeshFileConverter.log")

          print("Start Reading...", inputfilename)

          from BasicTools.IO.IOFactory import CreateReader
          if ops["printTimes"]:
              from BasicTools.IO.IOFactory import InitAllReaders
              InitAllReaders()
              import os
              basename,extention = os.path.splitext(os.path.basename(inputfilename))

              reader = CreateReader("."+inputfilename.split(".")[-1].lower())
              reader.SetFileName(inputfilename)
              if reader.canHandleTemporal :
                  reader.ReadMetaData()
                  print("Available Times in files:")
                  print(reader.GetAvailableTimes())
                  import sys
                  sys.exit(0)

          if ops["timeToRead"] == -10:
            reader = CreateReader("."+inputfilename.split(".")[-1].lower())
            reader.SetFileName(inputfilename)
            if reader.canHandleTemporal :
                reader.ReadMetaData()
            times = reader.GetAvailableTimes()
          else:
            mesh = ReadMesh(inputfilename,timeToRead = ops["timeToRead"])
            if ops["PlotOnScreen"]:
                from BasicTools.Bridges.vtkBridge import PlotMesh
                PlotMesh(mesh)
            if len(outputfilename) == 0:
                from BasicTools.Containers.UnstructuredMeshInspectionTools import PrintMeshInformation
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

          if ops["timeToRead"] == -10:
              writer = CreateWriter("."+outputfilename.split(".")[-1])
              writer.SetFileName(outputfilename)
              writer.SetTemporal()
              writer.SetBinary(ops["binary"])
              writer.Open(outputfilename)
              reader = CreateReader("."+inputfilename.split(".")[-1].lower())
              reader.SetFileName(inputfilename)
              reader.ReadMetaData()

              for t in times:
                print("Treating time : {}".format(t))
                reader.SetTimeToRead(t)

                mesh = reader.Read()

                PointFields = None
                PointFieldsNames = None
                if hasattr(mesh,"nodeFields"):
                    PointFieldsNames = list(mesh.nodeFields.keys())
                    PointFields = list(mesh.nodeFields.values())

                CellFields = None
                CellFieldsNames = None
                if hasattr(mesh,"elemFields"):
                    CellFieldsNames = list(mesh.elemFields.keys())
                    CellFields = list(mesh.elemFields.values())

                writer.Write(mesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields, Time=t )
              writer.Close()
          else:
            WriteMesh(outputfilename,mesh,writer=writer,binary=ops["binary"])
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
            for binary in [True,False]:
                inputfilename = GetTestDataPath() + iff
                outputfilename = TestTempDir().GetTempPath()+iff+"." + off
                ops= {}
                ops["timeToRead"] = -1
                ops["printTimes"] = False
                ops["PlotOnScreen"] = GUI
                ops["OnlyAbaqusReader"] = False
                ops["OnHelp"] = False
                ops["MeshIO"] = False
                ops["binary"] = binary
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
          opts, args = getopt.getopt(sys.argv[1:],"bsvphmaTt:i:o:")
      #except getopt.GetoptError:
      #    PrintHelp()
      #    sys.exit(2)

      outputfilename = ""
      ops= {}
      ops["timeToRead"] = -1
      ops["printTimes"] = False
      ops["PlotOnScreen"] = False
      ops["OnlyAbaqusReader"] = False
      ops["OnHelp"] = False
      ops["MeshIO"] = False
      ops["binary"] = None
      for opt, arg in opts:
         if opt == '-h':
             ops["OnHelp"] = True
         elif opt in ("-i",):
            inputfilename = arg
         elif opt in ("-o",):
            outputfilename = arg
         elif opt in ("-t",):
            ops["timeToRead"] = float(arg)
         elif opt in ("-T",):
            ops["timeToRead"] = -10
         elif opt in ("-p",):
            ops["printTimes"] = bool(True)
         elif opt in ("-v",):
            from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
            BOO.SetGlobalDebugMode(True)
         elif opt in ("-c",):
            print(CheckIntegrity(GUI=True))
            sys.exit(0)
         elif opt in ("-m",):
             ops["MeshIO"] = True

         elif opt in ("-s",):
             ops["PlotOnScreen"] = True
         elif opt in ("-a",):
             ops["OnlyAbaqusReader"] = True
         elif opt in ("-b",):
             ops["binary"] = True

    if ops["OnHelp"]:
        PrintHelp(ops)
        sys.exit()
    else:
        Convert(inputfilename,outputfilename,ops )

if __name__ == '__main__' :
    Main()
