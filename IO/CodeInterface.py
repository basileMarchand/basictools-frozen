# -*- coding: utf-8 -*-
""" Class to help the execution of an external program
"""

import subprocess
import os
__author__ = "Felipe Bordeu, Fabien Casenave"

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
        self.lastCommandExecuted = None
        self.openExternalWindows = True
        self.keepExternalWindows = True

        self.withFilename = True

    def WriteFile(self, idProc):

        # Write code input file
        try:
            inpString = self.tpl.format(**self.parameters)
        except KeyError as e: # pragma: no cover
            print("The user must supply the key: %s" % str(e))
            raise

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension

        with open(self.processDirectory + os.sep + inpFilename, 'w') as inpFile:

            inpFile.write(inpString)
    def GenerateCommandToRun(self,idProc=0):
        # Command to execute
        cmd = self.codeCommand

        if self.withFilename:
            inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension
            cmd += ' ' + inpFilename

        cmd = cmd.format(**self.parameters)
        return cmd

    def SingleRunComputation(self, idProc,stdout = None):


        cmd = self.GenerateCommandToRun(idProc)

        if stdout is None:
            out = open(os.devnull, 'w')
        else:
            out = stdout

        # Commande execution

        if self.openExternalWindows:
            if self.keepExternalWindows:
                self.lastCommandExecuted = "/usr/bin/xterm -e '"+ cmd + ";bash'"
            else:
                self.lastCommandExecuted = "/usr/bin/xterm -e '"+ cmd + "'"
        else:
            self.lastCommandExecuted = cmd;

        proc = subprocess.Popen(cmd , cwd=self.processDirectory , stdout=out, shell=True)

        return proc

    def SingleRunComputationAndReturnOutput(self, idProc=0):

        return proc


        cmd = self.GenerateCommandToRun(idProc)

        # Commande execution
        self.lastCommandExecuted = cmd;
        print(cmd)
        out = subprocess.check_output(cmd, cwd=self.processDirectory, shell=True ).decode("utf-8","ignore")

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


def CheckIntegrity():



    import BasicTools.Helpers.Tests as T
    import BasicTools.TestData as BasicToolsTestData

    interface = Interface(BasicToolsTestData.GetTestDataPath())

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

    interface.SetCodeCommand("dir ")
    interface.ReadTemplateFile('template.tpl')
    interface.WriteFile(1)

    import sys
    interface.SingleRunComputation(1,sys.stdout).wait()
    print("lastCommandExecuted: " + str(interface.lastCommandExecuted))

    interface.SetCodeCommand("dir {filter}")
    interface.parameters['filter']        = '*.inp'
    interface.withFilename = False


    print("output is :" + str(interface.SingleRunComputationAndReturnOutput(1).encode("ascii","ignore") ))
    print("lastCommandExecuted: " + str(interface.lastCommandExecuted))
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover


