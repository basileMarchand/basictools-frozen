# -*- coding: utf-8 -*-

from OTTools.Helpers.TextFormatHelper import TFormat
import numpy as np

def ArrayToString(data):
    return " ".join(str(x) for x in data)

def WriteMeshToXdmf(filename, baseMeshObject, PointFields = None, CellFields = None,GridFields= None, PointFieldsNames = None, CellFieldsNames=None, GridFieldsNames=None ):
    
    if PointFields is None:
        PointFields = [];
     
    if CellFields  is None:
        CellFields   = [];
    
    if GridFields is None:
        GridFields  = [];
    
    if PointFieldsNames is None:
        PointFieldsNames  = [];
    
    if CellFieldsNames is None:
        CellFieldsNames  = [];
    
    if GridFieldsNames is None:
        GridFieldsNames  = [];

    writer = XdmfWriter(filename)
    writer.SetBinary()
    writer.Open()
    writer.Write(baseMeshObject, 
                 PointFields= PointFields,
                 CellFields = CellFields,
                 GridFields = GridFields,
                 PointFieldsNames = PointFieldsNames,
                 CellFieldsNames = CellFieldsNames,
                 GridFieldsNames = GridFieldsNames 
                 )    
    writer.Close()
#

class XdmfWriter:
    def __init__(self, fileName = None):
        self.fileName = None;
        self.timeSteps = [];
        self.currentTime = 0;
        self.__XmlSizeLimit = 10
        
        self.__isBinary = False;
        self.__isTemporalOutput = False
        self.__binFileName = None;
        self.__filePointer = None;        
        self.__isOpen = False;
        self.__binarycpt = 0;
        
        self.SetFileName(fileName)
        
    def __str__(self):
        res  = 'XdmfWriter : \n'
        res += '   FileName : '+self.fileName+'\n'
        if self.__isBinary:
           res += '   Binary output Active \n'
           res += '   Binary FileName : '+ self.__binFileName +'\n'
        if self.__isTemporalOutput:
           res += '   Temporal output Active \n'
           res += '   TimeSteps : '+ str(self.timeSteps) + '\n'
        if self.__isOpen:
           res += '   The File is Open!! \n' 
        return res

    def SetFileName(self, fileName ):
        
        if fileName is None :
            self.fileName = None;    
            self.__path = None;    
            return 

        #print('Setting filename '+ fileName)
        import os
        
        self.fileName = fileName;
        self.__path  = os.path.abspath(os.path.dirname(fileName));
        self.__binFileName = os.path.splitext(os.path.abspath(fileName))[0] + ".bin"
        self.__binFileNameOnly = os.path.basename(self.__binFileName)

    def Step(self, dt = 1):
        self.currentTime += dt;
        self.timeSteps.append(self.currentTime);
    
    def SetTemporal(self, val = True):
        if self.__isOpen :
            print(TFormat.InRed("SetTemporal before opening"))
            raise Exception
        self.__isTemporalOutput = val
        
    def SetBinary(self, val = True):
        if self.__isOpen :
            print(TFormat.InRed("SetBinary before opening"))
            raise Exception
        self.__isBinary = val
        
    def SetXmlSizeLimit(self,val):
        self.__XmlSizeLimit= val
        
    def Open(self, filename = None):
        
        if self.__isOpen :
            print(TFormat.InRed("The file is already open !!!!!"))
            raise Exception
            
        
        if filename is not None:
            self.SetFileName(filename)
        
        ## we use unbuffered so we can repaire broken files easily 
        try :
            self.filePointer = open(self.fileName, 'w',0)
        except:
            print(TFormat.InRed("Error File Not Open"))
            raise  
        
        #print("File : '" + self.fileName + "' is open.")         
        self.__isOpen = True
        self.filePointer.write('<?xml version="1.0" encoding="utf-8"?>\n')
        self.filePointer.write('<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.92">\n')
        self.filePointer.write('<Domain>\n')
        
        if self.__isTemporalOutput:
            self.filePointer.write('<Grid Name="Grid_T" GridType="Collection" CollectionType="Temporal"  >\n')
        #clean the binary file  =
        if self.__isBinary:
            self.__binaryFilePointer = open (self.__binFileName, "wb")
                
    def Close(self):
        if self.__isOpen:
            if self.__isTemporalOutput:
                self.__WriteTime()
                self.filePointer.write('    </Grid>\n')  
            self.filePointer.write('  </Domain>\n')
            self.filePointer.write('</Xdmf>\n')
            self.filePointer.close()
            self.__isOpen = False
            if self.__isBinary:
                self.__binaryFilePointer.close()
            #print("File : '" + self.fileName + "' is close.")         
    
    def __WriteGeoAndTopo(self,baseMeshObject):
        if baseMeshObject.isConstantRectilinear() :
            origin = baseMeshObject.GetOrigin()
            spacing = baseMeshObject.GetSpacing()
            dims = baseMeshObject.GetDimensions() ## number of nodes per 
            self.filePointer.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
            self.filePointer.write('      <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'+ArrayToString(reversed(origin)) +'</DataItem>\n')
            self.filePointer.write('      <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'+ArrayToString(reversed(spacing)) +'</DataItem>\n')
            self.filePointer.write('    </Geometry>\n')
            self.filePointer.write('    <Topology Dimensions="'+ArrayToString(reversed(dims))  +'" Type="3DCoRectMesh"/>\n')
        elif baseMeshObject.isRectilinear() :# pragma: no cover 
            print(TFormat.InRed("Mesh Type Rectilinear Not Supported"))        # pragma: no cover 
            raise Exception                                                    # pragma: no cover 
        elif baseMeshObject.isStructured() :                                   # pragma: no cover 
            print(TFormat.InRed("Mesh Type Structured Not Supported"))         # pragma: no cover 
            raise Exception                                                    # pragma: no cover 
        elif baseMeshObject.isUnstructured() :                                 # pragma: no cover 
            print(TFormat.InRed("Mesh Type Unstructured Not Supported"))       # pragma: no cover 
            raise Exception                                                    # pragma: no cover 
        else:                                                                  # pragma: no cover 
            print(TFormat.InRed("Mesh Type Not Supported"))                    # pragma: no cover 
            raise Exception                                                    # pragma: no cover 

    def __WriteAttribute(self,data,name,center,shape,baseMeshObject):
        
       ndata = np.prod(shape)
       if data.size == ndata:
           attype = "Scalar"
           #print(shape)
           #print(data.shape)
           if baseMeshObject.isConstantRectilinear() and len(shape)>1 :
           
               if len(data.shape) <= 2:
                   data.shape = tuple(shape)
                   data = data.T
               else:
                   data = data.T
                   
       elif data.size == ndata*3:
           
           attype = "Vector"
           #print(attype)
           #print(data.shape)
           #print(shape)
           if baseMeshObject.isConstantRectilinear() and len(shape)>1:
               shape = (shape[0], shape[1],shape[2],3)
               if len(data.shape) <= 2:
                   shape = (shape[0], shape[1],shape[2],3)
                   data.shape = shape
                   data = data.transpose(2,1,0,3)
               else:
                   #print(data.shape)
                   data = data.transpose(2,1,0,3)
       else:
           print(TFormat.InRed("I dont kow how to treak fields with " +str(data.size/ndata) +" components"))  # pragma: no cover       
           raise Exception                                                                                    # pragma: no cover
           
       self.filePointer.write('    <Attribute Center="'+center+'" Name="'+name+'" Type="'+attype+'">\n')#

       self.__WriteDataItem(data.flatten(),shape  )

       self.filePointer.write('    </Attribute>\n')

        
    def Write(self,baseMeshObject, PointFields = None, CellFields = None, GridFields= None, PointFieldsNames = None, CellFieldsNames= None, GridFieldsNames=None , Time= None, TimeStep = None ):
        
         if PointFields is None:
             PointFields = [];
             
         if CellFields  is None:
            CellFields   = [];
            
         if GridFields is None:
            GridFields  = [];
            
         if PointFieldsNames is None:
            PointFieldsNames  = [];
            
         if CellFieldsNames is None:
            CellFieldsNames  = [];
            
         if GridFieldsNames is None:
            GridFieldsNames  = [];
            
         if not self.__isOpen :
            print(TFormat.InRed("Please Open The writer First!!!"))
            raise Exception
        
         dt = 1
         if Time is not None:
             dt = Time - self.currentTime
         elif TimeStep is not None:
             dt = TimeStep
             
         if self.__isTemporalOutput:
             self.Step(dt)
             
         self.filePointer.write('    <Grid Name="Grid_'+str(len(self.timeSteps))+'">\n')
         
         self.__WriteGeoAndTopo(baseMeshObject)

         
         for i in range(len(PointFields)): 
           name = 'PField'+str(i)
           if len(PointFields)  == len(PointFieldsNames):
               name = PointFieldsNames[i]
               
           self.__WriteAttribute(np.array(PointFields[i]), name, "Node",baseMeshObject.GetDimensions(),baseMeshObject)
           
         for i in range(len(CellFields)): 
           name = 'CField'+str(i)
           if len(CellFields) == len(CellFieldsNames):
               name = CellFieldsNames[i]

           self.__WriteAttribute(np.array(CellFields[i]), name, "Cell",baseMeshObject.GetDimensions()-1,baseMeshObject)
           
          
           
#         for i in range(len(CellFields)): 
#           name = 'CField'+str(i)
#           if len(CellFields) == len(CellFieldsNames):
#               name = CellFieldsNames[i]
#           
#           self.filePointer.write('    <Attribute Center="Cell" Name="'+name+'" Type="Scalar">\n')
#           self.filePointer.write('      <DataItem DataType="Float" Dimensions="'+" ".join(str(x-1) for x in np.flipud(baseMeshObject.GetDimensions())) +'" Format="XML" Precision="4">\n')
#           field = CellFields[i].view()
#           field.shape = tuple((x-1) for x in baseMeshObject.GetDimensions())
#           self.filePointer.write(" ".join(str(x) for x in field.T.flatten()))
#           self.filePointer.write('    </DataItem>\n')
#           self.filePointer.write('    </Attribute>\n')
           
           
           
         for i in range(len(GridFields)): 
             
           name = 'GField'+str(i)
           if len(GridFields) == len(GridFieldsNames):
               name = GridFieldsNames[i]
               
           self.__WriteAttribute(np.array(GridFields[i]), name, "Grid",[1],baseMeshObject)
           
          
         self.filePointer.write('    </Grid>\n')  
            
    def __WriteTime(self):
        
        if self.__isOpen:
            #self.filePointer.write('<Time TimeType="List">\n')
            #self.__WriteDataItem(self.timeSteps)        
            #self.filePointer.write('</Time>\n')
            self.filePointer.write('<Time Value="'+ (" ".join(str(x) for x in self.timeSteps)) +'"/>\n')
            
    def __WriteDataItem(self,_data, _shape):
        import numpy as np
        data = np.array(_data)
        shape = np.array(_shape)
        if self.__isOpen:
            if type(data[0]) == np.float64:
                typename = 'Float'
                s = data.dtype.itemsize
            elif type(data[0]) == np.float32:
                typename = 'Float'
                s = data.dtype.itemsize
            elif type(data[0]) == np.int32:
                typename = 'Int'
                s = data.dtype.itemsize
            elif type(data[0]) == np.int64:
                typename = 'Int'
                s = data.dtype.itemsize
            else:
                print('Output Not implemented for data of type ')              # pragma: no cover 
                print(type(data[0]))                                           # pragma: no cover 
                raise                                                          # pragma: no cover 
                
            dimension = ArrayToString(shape)
            
            if self.__isBinary and len(data) > self.__XmlSizeLimit:
                #bindata = bytearray(data.flatten())
                
                #self.__binaryFilePointer.write(bindata)  
                data.flatten().tofile(self.__binaryFilePointer)
                
                self.filePointer.write(' <DataItem Format="Binary"'+
                ' NumberType="'+typename+'"'+
                ' Dimensions="'+dimension+'" '+
                ' Seek="'+str(self.__binarycpt)+'" '+
                ' Endian="Native" '+
                ' Precision="'+str(s)+'" '+

                ' Compression="Raw" >')               
                self.filePointer.write(self.__binFileNameOnly)
                self.__binarycpt += s*len(data)
            else:
                self.filePointer.write(' <DataItem Format="XML" NumberType="'+typename+'" Dimensions="'+dimension+'">')
                self.filePointer.write(" ".join(str(x) for x in data.flatten()))

                
            self.filePointer.write('</DataItem>\n')


         


def WriteTest(tempdir,Temporal, Binary):
     
    from OTTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh
    
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,3,4]);
    myMesh.SetSpacing([0.1, 0.1, 0.1]);
    myMesh.SetOrigin([-2.5,-1.2,-1.5]);
    
    dataT = np.array(range(24),dtype=np.float32)
    dataT.shape = (2,3,4)
    dataDep = np.array(range(24*3))+0.3
    
    dataDep.shape = (2,3,4,3)
    
    writer = XdmfWriter(tempdir + '\TestOutput_Bin_'+str(Binary)+'_Temp_'+str(Temporal)+'.xmf')
    writer.SetTemporal(Temporal)
    writer.SetBinary(Binary)
    writer.Open()
    print(writer)
    writer.Write(myMesh,PointFields=[dataT, dataDep], PointFieldsNames=["Temp","Dep"],CellFields=[np.array(range(6))],CellFieldsNames=['S'], Time=0);
    writer.Write(myMesh,GridFields=[0, 1], GridFieldsNames=['K','P'], TimeStep = 1);
    writer.Close()
   
def CheckIntegrity():
    from OTTools.Helpers.Tests import TestTempDir
    from OTTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh


    tempdir = TestTempDir.GetTempPath()
    
    WriteTest(tempdir,False, False)
    WriteTest(tempdir,False, True)
    WriteTest(tempdir,True, False)
    WriteTest(tempdir,True, True)
    
    WriteMeshToXdmf(tempdir+'testdirect.xdmf',ConstantRectilinearMesh() )
    
    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetFileName(tempdir+'testerros.xdmf')
    writer.SetXmlSizeLimit(0)
    writer.Open();

    ## test of the erros
    try:
        writer.SetTemporal()
        return 'Not ok'# pragma: no cover 
    except:
        pass

    try:
        writer.SetBinary()
        return 'Not ok'# pragma: no cover 
    except:
        pass
    
    
    
    try:
        writer.Open();
    except:
        pass
    
    writer.Close()
    
    try:
        writer.Write(ConstantRectilinearMesh())
        return "Not ok"# pragma: no cover 
    except:
        pass
    
    # for the moment we comment this test (verification of "unable to open file")
    # need to find a way to raise an exception in linux and in windows with the 
    # same filename
    #
    #try:
    #    writer.Open("c:\\")  # in windows this will raise an exception
    #    writer.Open("\")     # in linux this will raise an exception
    #    return 'Not ok'# pragma: no cover 
    #except:
    #    pass
    
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover 
    