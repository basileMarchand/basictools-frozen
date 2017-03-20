# -*-coding:Latin-1 -*
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
        self.lastCommandExecuted = None

    def WriteFile(self, idProc):

        # Write code input file
        try:
            inpString = self.tpl.format(**self.parameters)
        except KeyError as e: # pragma: no cover
            print("The user must supply the key: %s" % str(e))
            raise

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension
        with open(self.processDirectory + inpFilename, 'w') as inpFile:
            inpFile.write(inpString)

    def SingleRunComputation(self, idProc,stdout = None):# pragma: no cover

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension

        # Command to execute
        cmd = 'cd ' + self.processDirectory + ' && ' + self.codeCommand + ' ' + inpFilename

        if stdout is None:
            out = open(os.devnull, 'w')
        else:
            out = stdout

        # Commande execution
        self.lastCommandExecuted = cmd;
        proc = subprocess.Popen(cmd , stdout=out, shell=True)

        return proc

    def SingleRunComputationAndReturnOutput(self, idProc):# pragma: no cover

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension

        # Command to execute
        cmd = self.codeCommand + ' ' + inpFilename

        # Commande execution
        self.lastCommandExecuted = cmd;
        out = subprocess.check_output(cmd, cwd=self.processDirectory, shell=True ).decode("utf-8")

        return out

    def ReadFile(self, filenameDir):

        # Template file read
        with open(filenameDir, 'r') as File:
            string = File.read()

        return string


    def SetWorkingDirectory(self,Dir):
        self.workingDirectory = Dir;

    def SetProcessDirectory(self,Dir):
        self.processDirectory = Dir;

    def SetCodeCommand(self,ccommand):          # Code command
        self.codeCommand = ccommand

    def SetTemplateFile(self,filename):
        self.tplFilename = filename;

    def ReadTemplateFile(self, filename = None):
        if filename is not None:
            self.SetTemplateFile(filename)
        self.tpl = open(self.workingDirectory + os.sep + self.tplFilename).read()



def CheckIntegrity():
    import BasicTools.IO.CodeInterface as CI
    import BasicTools.Helpers.Tests as T
    import BasicTools.TestData as T2
    dataPath = T2.GetTestDataPath()

    interface = CI.Interface(dataPath)

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

    interface.SetProcessDirectory(T.TestTempDir.GetTempPath())
    interface.SetCodeCommand("ls -l ")
    #interface.SetTemplateFile('template.tpl')
    interface.ReadTemplateFile('template.tpl')
    interface.WriteFile(1)
    import sys
    interface.SingleRunComputation(1,sys.stdout).wait()

    print("output is :" + str(interface.SingleRunComputationAndReturnOutput(1).encode("ascii","ignore") ))
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover


