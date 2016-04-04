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
        self.parameterNames  = []
        self.parameterValues = []
        
        # Template
        self.tplDirectory = self.workingDirectory + os.sep
        if not os.path.exists(self.tplDirectory): os.mkdir(self.tplDirectory)
        self.tplFilename = 'template.tpl'
        self.tpl = self.ReadFile(self.tplDirectory + self.tplFilename)
        
        # Temporary files folder creation
        self.processDirectory = self.tplDirectory
        if not os.path.exists(self.processDirectory): os.mkdir(self.processDirectory)
        
        # Output file name
        self.inputFilename      = 'calcul'
        self.inputFileExtension = '.inp'
        
        # Code command
        self.codeCommand = 'Zrun'
        
    def WriteFile(self, idProc):
        
        parametersDict = dict((name, self.parameterValues[i]) for i, name in enumerate(self.parameterNames))
        
        # Write code input file
        inpString = self.tpl.format(**parametersDict)
        
        inpFilename = self.inputFilename + str(idProc) + self.inputFileExtension
        
        with open(self.processDirectory + inpFilename, 'w') as inpFile:
            inpFile.write(inpString)
        
    def SingleRunCode(self, idProc):# pragma: no cover
      
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
    interface = CI.Interface(T.GetTestDataPath())
    interface.parameterNames  = ['calcul', 'Ti', 'mesh', 'solver', 'python_script', 'sequence', 'time', 'increment', 'iteration', 'ratio', 'algorithm', 'table', 'bc', 'conductivity', 'coefficient']
    interface.parameterValues = ['thermal_transient', 1000.0, 'Cube_3D.geof', 'mumps', 'reduction', 1, 2000., 20, 1000, 0.001, 'p1p2p3', '**name ptab\n*time  0.0  10000000.0  \n*value  1.0  1.0', '**surface_heat_flux\nface_x_0  1000.0  ptab\nface_x_1  1000.0  ptab', '280.+100.*atan(0.01*(1100-temperature));', '8.*(430.+40.*atan(0.01*(temperature-500.)))*(1.5-0.5*exp(-200./((temperature-1200.)*(temperature-1200.))));']
    interface.processDirectory  = T.TestTempDir.GetTempPath()
    interface.WriteFile(1)
    #interface.SingleRunCode(1)
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
