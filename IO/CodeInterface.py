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
        self.tpl = self.ReadFile(self.workingDirectory + os.sep + self.tplFilename)
        
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
      
    def WriteCubeZebulonMesher(self, mu1, mu2, mu3, nCuts):
        with open(self.processDirectory + 'cube.mast', 'w') as mast_file:
            mast_file.write("****master\n***geometry\n**point point0\n*position 0. 0. 0.\n**point point1\n*position "+str(mu1)+" 0. 0.")
            mast_file.write("\n**point point2\n*position  "+str(mu1)+" "+str(mu2)+" 0.\n**point point3\n*position 0. "+str(mu2)+" 0.\n")
            mast_file.write("**line line0\n*p1 point0\n*p2 point1\n*cuts "+str(nCuts)+"\n**line line1\n*p1 point1\n*p2 point2\n*cuts "+str(nCuts)+"\n")
            mast_file.write("**line line2\n*p1 point2\n*p2 point3\n*cuts "+str(nCuts)+"\n**line line3\n*p1 point3\n*p2 point0\n*cuts "+str(nCuts)+"\n")
            mast_file.write("**domain domain1 Ruled0\n*name R\n*element c2d4\n*method 1\n*b1 line0\n*b2 line1\n*b3 line2\n*b4 line3\n")
            mast_file.write("***mesher\n**extension\n*elset ALL_ELEMENT\n*prog 1.00000\n*fusion 0.00100000\n*distance "+str(mu3)+"\n*num "+str(nCuts)+"\n")
            mast_file.write("*dir ( 0.00000 0.00000 1.00000 )\n***plotting\n**default_options 1 1 1\n****return")
        mast_file.close()  
      
    def WriteZebulonTablesBC(self):
        self.parameters['table'] = ""
        for i in xrange(len(self.ptab['names'])):
          self.parameters['table']  = self.parameters['table'] + '**name '+self.ptab['names'][i]+'\n'
          self.parameters['table']  = self.parameters['table'] + '*time  '+str(self.ptab['times'][i][0])+'  '+str(self.ptab['times'][i][1])+'\n'
          self.parameters['table']  = self.parameters['table'] + '*value  '+str(self.ptab['values'][i][0])+'  '+str(self.ptab['values'][i][1])+'\n'
     
        self.parameters['bc'] = "**surface_heat_flux\n"
        for i in xrange(len(self.bc['surface_heat_flux']['bsets'])):
          self.parameters['bc']  = self.parameters['bc']+self.bc['surface_heat_flux']['bsets'][i]+'  '+str(self.bc['surface_heat_flux']['values'][i])+'  '+self.bc['surface_heat_flux']['ptabs'][i]+'\n'         
                   

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

    interface.ptab['names']               = ['ptab']
    interface.ptab['times']               = [[0., 10000000.]]
    interface.ptab['values']              = [[1., 1.]]

    interface.bc['surface_heat_flux']     = {}
    interface.bc['surface_heat_flux']['bsets']  = ['face_x_0', 'face_x_1']
    interface.bc['surface_heat_flux']['values'] = [1000., 1000.]
    interface.bc['surface_heat_flux']['ptabs']  = ['ptab', 'ptab']

    interface.WriteZebulonTablesBC()
        
    interface.parameters['conductivity']  = '280.+100.*atan(0.01*(1100-temperature));'
    interface.parameters['coefficient']   = '8.*(430.+40.*atan(0.01*(temperature-500.)))*(1.5-0.5*exp(-200./((temperature-1200.)*(temperature-1200.))));'
    
    interface.processDirectory  = T.TestTempDir.GetTempPath()
    interface.WriteFile(1)
    #interface.SingleRunComputation(1)
    interface.WriteCubeZebulonMesher(1., 1., 1., 3)
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   

    
