# -*-coding:Latin-1 -*
import subprocess
import numpy as np
import os
import time

class Interface(object):

    def __init__(self, workingDirectory = os.getcwd()):

        # Working Folder
        """Working folder must contain:
            a folder 'self.tpl_directory' containing :
                'self.tpl_filename'
        """
        self.workingDirectory = workingDirectory

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

    def WriteFile(self, idProc):

        # Write code input file
        inpString = self.tpl.format(**self.parameters)

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension
        with open(self.processDirectory + inpFilename, 'w') as inpFile:
            inpFile.write(inpString)

    def SingleRunComputation(self, idProc):# pragma: no cover

        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension

        # Command to execute
        cmd = 'cd ' + self.processDirectory + ' && ' + self.codeCommand + ' ' + inpFilename

        # Commande execution
        FNULL = open(os.devnull, 'w')
        proc = subprocess.Popen(cmd , stdout=FNULL, shell=True)

        return proc

    def ReadFile(self, filenameDir):

        # Template file read
        with open(filenameDir, 'r') as File:
            string = File.read()

        return string


def CheckIntegrity():
    import OTTools.IO.CodeInterface as CI
    import OTTools.Helpers.Tests as T
    import OTTools.TestData as T2
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

    interface.processDirectory  = T.TestTempDir.GetTempPath()
    interface.WriteFile(1)
    #interface.SingleRunComputation(1)
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover


