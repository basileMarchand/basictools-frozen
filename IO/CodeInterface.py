# -*-coding:Latin-1 -*
import subprocess
import numpy as np
import os
import time

class Interface:
    
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
        
    def SingleRunCode(self, idProc):
      
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
    interface.parameterNames  = ['a', 'b']
    interface.parameterValues = [1., 2.]
    interface.WriteFile(1)
    #interface.SingleRunCode(1)
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
