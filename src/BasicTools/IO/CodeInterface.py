# -*- coding: utf-8 -*-
""" Class to help the execution of an external program
"""

import subprocess
import os

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class Interface(BaseOutputObject):

    def __init__(self, workingDirectory = os.getcwd()):
        super(Interface,self).__init__()

        # Working Folder
        """Working folder must contain:
            a folder 'self.tpl_directory' containing :
                'self.tpl_filename'
        """
        self.SetWorkingDirectory(workingDirectory)

        # Model parameters
        self.parameters = {}
        self.ptab       = {}
        self.bc         = {}


        # Template
        self.tplFilename = 'template.tpl'
        try:
          self.tpl = self.ReadFile(self.workingDirectory + os.sep + self.tplFilename)
        except IOError:# pragma: no cover
          True

        # Temporary files folder creation
        self.processDirectory = self.workingDirectory + os.sep

        # Output file name
        self.inputFilename      = 'calcul'
        self.inputFileExtension = '.inp'

        # Code command
        self.codeCommand = 'Zrun'
        self.options = []
        self.lastCommandExecuted = None
        self.openExternalWindows = False
        self.keepExternalWindows = False

        self.withFilename = True

    def WriteFile(self, idProc):

        # Write code input file
        try:
            def expand_vars(string, vars):
                while True:
                    expanded = string.format(**vars)
                    if expanded == string:
                        break
                    string = expanded
                return string
            inpString = expand_vars( self.tpl,self.parameters)
            #inpString = self.tpl.format(**self.parameters)
        except KeyError as e: # pragma: no cover
            print("The user must supply the key: %s" % str(e))
            raise

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension

        with open(self.processDirectory + os.sep + inpFilename, 'w') as inpFile:

            inpFile.write(inpString)

    def SetOptions(self, opts):
        self.options = opts

    def GenerateCommandToRun(self,idProc=0):
        # Command to execute
        cmd = []
        cmd.append(self.codeCommand)

        for opt in self.options:
            cmd.append(opt)

        if self.withFilename:
            inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension
            cmd.append(inpFilename)


        for i in range(len(cmd)):
           cmd[i] = cmd[i].format(**self.parameters)
        return cmd

    def SingleRunComputation(self, idProc,stdout = None):


        cmd = self.GenerateCommandToRun(idProc)

        if stdout is None:
            out = open(os.devnull, 'w')
        else:
            out = stdout

        # Commande execution

        if self.openExternalWindows:
            import platform
            if platform.system() == "Windows":
                if self.keepExternalWindows:
                   cmd.insert(0,"/K")
                else:
                   cmd.insert(0,"/C")
                cmd.insert(0,"cmd")
                cmd.insert(0,"/wait")
                cmd.insert(0,"start")
            else:
                if self.keepExternalWindows:
                    cmd = [ " ".join(cmd)+" ; bash " ]
                cmd.insert(0,"-e")
                cmd.insert(0,"/usr/bin/xterm")

        self.lastCommandExecuted = cmd
        print(cmd)
        proc = subprocess.Popen(cmd , cwd=self.processDirectory , stdout=out, shell=False)

        return proc

    def SingleRunComputationAndReturnOutput(self, idProc=0):

        cmd = self.GenerateCommandToRun(idProc)

        # Commande execution
        self.lastCommandExecuted = cmd;
        self.PrintVerbose(self.processDirectory)
        self.PrintVerbose(cmd)
        out = subprocess.check_output(cmd, cwd=self.processDirectory, shell=False ).decode("utf-8","ignore")

        return out

    def ReadFile(self, filenameDir):

        # Template file read
        with open(filenameDir, 'r') as File:
            string = File.read()

        return string


    def SetWorkingDirectory(self,Dir):
        self.workingDirectory = os.path.dirname(Dir);
        if len(self.workingDirectory) and self.workingDirectory[-1] != os.sep:
            self.workingDirectory += os.sep

    def SetProcessDirectory(self,Dir):
        self.processDirectory = os.path.dirname(Dir);
        if len(self.processDirectory) and self.processDirectory[-1] != os.sep:
            self.processDirectory += os.sep

    def SetCodeCommand(self,ccommand):          # Code command
        self.codeCommand = ccommand

    def SetTemplateFile(self,filename):
        self.tplFilename = filename;

    def ReadTemplateFile(self, filename = None):
        if filename is not None:
            self.SetTemplateFile(filename)
        self.tpl = open(self.workingDirectory + os.sep + self.tplFilename).read()

    def CopyFile(self,filetocopy):
        import shutil

        shutil.copy(self.workingDirectory + os.sep + filetocopy,
                   self.processDirectory + os.sep +filetocopy.split('/') [-1])


def CheckIntegrity(GUI=False):



    import BasicTools.Helpers.Tests as T
    import BasicTools.TestData as BasicToolsTestData

    interface = Interface(BasicToolsTestData.GetTestDataPath())
    interface.SetGlobalDebugMode(True)
    interface.keepExternalWindows = GUI
    interface.parameters['calcul']        = 'thermal_transient'
    interface.parameters['Ti']            = 1000.0
    interface.parameters['mesh']          = 'Cube_3D.geof'
    interface.parameters['solver']        = 'mumps'
    interface.parameters['python_script'] = 'reduction'
    interface.parameters['sequence']      = 1
    interface.parameters['time']          = 2000.
    interface.parameters['increment']     = 20
    interface.parameters['iteration']     = 1000
    interface.parameters['ratio']         = 0.001
    interface.parameters['algorithm']     = 'p1p2p3'

    interface.parameters['bc']            = 'myBC'
    interface.parameters['table']         = 'myTable'

    interface.parameters['conductivity']  = '280.+100.*atan(0.01*(1100-temperature));'
    interface.parameters['coefficient']   = '8.*(430.+40.*atan(0.01*(temperature-500.)))*(1.5-0.5*exp(-200./((temperature-1200.)*(temperature-1200.))));'
    import sys

    interface.SetProcessDirectory(T.TestTempDir.GetTempPath())

    interface.SetCodeCommand("dir")
    interface.ReadTemplateFile('template.tpl')
    interface.WriteFile(1)

    interface.openExternalWindows = GUI
    interface.keepExternalWindows = GUI
    interface.SingleRunComputation(1,sys.stdout).wait()
    print("lastCommandExecuted: " + str(interface.lastCommandExecuted))


    interface.SetCodeCommand("dir")
    interface.SetOptions(["{filter}"])
    interface.parameters['filter']        = 'calcul1.inp'
    interface.withFilename = False


    print("output is :" + str(interface.SingleRunComputationAndReturnOutput(1).encode("ascii","ignore") ))
    print("lastCommandExecuted: " + str(interface.lastCommandExecuted))
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover

