# -*- coding: utf-8 -*-

from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.IO.WriterBase import WriterBase as WriterBase

import numpy as np
import os
import BasicTools.FE.ElementNames as EN

def ArrayToString(data):
    return " ".join(str(x) for x in data)

XdmfName = {}
XdmfName[EN.Point_1] = 'Polyvertex'
XdmfName[EN.Bar_2] = 'Polyline'
XdmfName[EN.Triangle_3] = 'Triangle'
XdmfName[EN.Quadrangle_4] = 'Quadrilateral'
XdmfName[EN.Tetrahedron_4] ="Tetrahedron"
XdmfName[EN.Pyramid_5] = 'Pyramid'
XdmfName[EN.Wedge_6] = 'Wedge'
XdmfName[EN.Hexaedron_8] = 'Hexahedron'


XdmfName[EN.Bar_3] = "Edge_3"
XdmfName[EN.Triangle_6] = 'Triangle_6'
XdmfName[EN.Quadrangle_9] = 'Quadrilateral_9'
XdmfName[EN.Quadrangle_8] = 'Quadrilateral_8'
XdmfName[EN.Tetrahedron_10] = 'Tetrahedron_10'
XdmfName[EN.Pyramid_13] = 'Pyramid_13'
XdmfName[EN.Wedge_15] = 'Wedge_13'
XdmfName[EN.Wedge_18] = 'Wedge_18'
XdmfName[EN.Hexaedron_20] = 'Hexahedron_20'

XdmfNumber = {}
XdmfNumber[EN.Point_1] = 0x1
XdmfNumber[EN.Bar_2] = 0x2
XdmfNumber[EN.Triangle_3] = 0x4
XdmfNumber[EN.Quadrangle_4] = 0x5
XdmfNumber[EN.Tetrahedron_4] = 0x6
XdmfNumber[EN.Pyramid_5] = 0x7
XdmfNumber[EN.Wedge_6] = 0x8
XdmfNumber[EN.Hexaedron_8] = 0x9

XdmfNumber[EN.Bar_3] = 0x22
XdmfNumber[EN.Triangle_6] = 0x24
XdmfNumber[EN.Quadrangle_9] = 0x23
XdmfNumber[EN.Quadrangle_8] = 0x25
XdmfNumber[EN.Tetrahedron_10] = 0x26
XdmfNumber[EN.Pyramid_13] = 0x27
XdmfNumber[EN.Wedge_15] = 0x28
XdmfNumber[EN.Wedge_18] = 0x29
XdmfNumber[EN.Hexaedron_20] = 0x30

#* Xdmf supports the following topology types:
# *   NoTopologyType
# *   Polyvertex - Unconnected Points
# *   Polyline - Line Segments
# *   Polygon - N Edge Polygon
# *   Triangle - 3 Edge Polygon
# *   Quadrilateral - 4 Edge Polygon
# *   Tetrahedron - 4 Triangular Faces
# *   Wedge - 4 Triangular Faces, Quadrilateral Base
# *   Hexahedron - 6 Quadrilateral Faces
# *   Edge_3 - 3 Node Quadratic Line
# *   Triangle_6 - 6 Node Quadratic Triangle
# *   Quadrilateral_8 - 8 Node Quadratic Quadrilateral
# *   Quadrilateral_9 - 9 Node Bi-Quadratic Quadrilateral
# *   Tetrahedron_10 - 10 Node Quadratic Tetrahedron
# *   Pyramid_13 - 13 Node Quadratic Pyramid
# *   Wedge_15 - 15 Node Quadratic Wedge
# *   Wedge_18 - 18 Node Bi-Quadratic Wedge
# *   Hexahedron_20 - 20 Node Quadratic Hexahedron
# *   Hexahedron_24 - 24 Node Bi-Quadratic Hexahedron
# *   Hexahedron_27 - 27 Node Tri-Quadratic Hexahedron
# *   Hexahedron_64 - 64 Node Tri-Cubic Hexahedron
# *   Hexahedron_125 - 125 Node Tri-Quartic Hexahedron
# *   Hexahedron_216 - 216 Node Tri-Quintic Hexahedron
# *   Hexahedron_343 - 343 Node Tri-Hexic Hexahedron
# *   Hexahedron_512 - 512 Node Tri-Septic Hexahedron
# *   Hexahedron_729 - 729 Node Tri-Octic Hexahedron
# *   Hexahedron_1000 - 1000 Node Tri-Nonic Hexahedron
# *   Hexahedron_1331 - 1331 Node Tri-Decic Hexahedron
# *   Hexahedron_Spectral_64 - 64 Node Spectral Tri-Cubic Hexahedron
# *   Hexahedron_Spectral_125 - 125 Node Spectral Tri-Quartic Hexahedron
# *   Hexahedron_Spectral_216 - 216 Node Spectral Tri-Quintic Hexahedron
# *   Hexahedron_Spectral_343 - 343 Node Spectral Tri-Hexic Hexahedron
# *   Hexahedron_Spectral_512 - 512 Node Spectral Tri-Septic Hexahedron
# *   Hexahedron_Spectral_729 - 729 Node Spectral Tri-Octic Hexahedron
# *   Hexahedron_Spectral_1000 - 1000 Node Spectral Tri-Nonic Hexahedron
# *   Hexahedron_Spectral_1331 - 1331 Node Spectral Tri-Decic Hexahedron
#  *   Mixed - Mixture of Unstructured Topologies


def WriteMeshToXdmf(filename, baseMeshObject, PointFields = None, CellFields = None,GridFields= None, PointFieldsNames = None, CellFieldsNames=None, GridFieldsNames=None , Binary= True):

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
    writer.SetBinary(Binary)
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

class XdmfWriter(WriterBase):
    def __init__(self, fileName = None):
        super(XdmfWriter,self).__init__()

        self.fileName = None;
        self.timeSteps = [];
        self.currentTime = 0;
        self.__XmlSizeLimit = 10
        self.automaticOpen = False;

        self.SetBinary(False)
        self.__isTemporalOutput = False
        self.__binFileName = None;
        self.__filePointer = None;
        #self.__isOpen = False;
        self.__binarycpt = 0;
        self.__binfilecounter = 0
        self.__keepXmlFileInSaneState = True;

        #set to off is you what to put the time at the end of the temporal grid
        #keep this option True to be compatible with the XMDF3 reader of Paraview
        self.__printTimeInsideEachGrid = True
        self.SetFileName(fileName)
        self.appendMode = False

    def SetAppendMode(self,mode = True):
        self.appendMode = mode

    def __str__(self):
        res  = 'XdmfWriter : \n'
        res += '   FileName : '+self.fileName+'\n'
        if self.isBinary():
           res += '   Binary output Active \n'
           res += '   Binary FileName : '+ self.__binFileName +'\n'
        if self.__isTemporalOutput:
           res += '   Temporal output Active \n'
           res += '   TimeSteps : '+ str(self.timeSteps) + '\n'
        if self.isOpen():
           res += '   The File is Open!! \n'
        return res

    def SetFileName(self, fileName ):

        if fileName is None :
            self.fileName = None;
            self.__path = None;
            return

        #print('Setting filename '+ fileName)


        self.fileName = fileName;
        self.__path  = os.path.abspath(os.path.dirname(fileName));
        self.binfilecounter = 0
        self.NewBinaryFilename()

    def NewBinaryFilename(self):

        self.__binFileName = os.path.splitext(os.path.abspath(self.fileName))[0] +"_" +str(self.__binfilecounter) +".bin"
        self.__binFileNameOnly = os.path.basename(self.__binFileName)
        self.__binfilecounter +=1

    def Step(self, dt = 1):
        self.currentTime += dt;
        self.timeSteps.append(self.currentTime);

    def SetTemporal(self, val = True):
        if self.isOpen() :
            print(TFormat.InRed("SetTemporal before opening"))
            raise Exception
        self.__isTemporalOutput = val


    def SetXmlSizeLimit(self,val):
        self.__XmlSizeLimit= val

    def Open(self, filename = None):

        # we dont use the open from WriterBase because the set binary is used
        # for the .bin file and not for the .xdmf file

        if self.isOpen() :
            print(TFormat.InRed("The file is already open !!!!!"))
            return
            #raise Exception


        if filename is not None:
            self.SetFileName(filename)

        ## we use unbuffered so we can repaire broken files easily
        try :
            # in python 3 we cant use unbuffered  text I/O (bug???)
            #self.filePointer = open(self.fileName, 'w',0)
            if self.appendMode:
                self.filePointer = open(self.fileName, 'r+')
                self.filePointer.seek(-100,2)

                lines = self.filePointer.readlines()
                #print("lines")
                for line in lines:
                    #print(line)
                    if line.find("<!--Temporal") != -1:
                         #pos (for seek)
                         l = line.find('pos="')
                         r = line.find('" ',l+1)
                         pos = int(line[l+5:r])


                         #__binfilecounter
                         l = line.find('ter="')
                         r = line.find('" ',l+1)
                         binfilecounter = int(line[l+5:r])

                         self.NewBinaryFilename()
                         self.__binaryFilePointer = open (self.__binFileName, "wb")
                         self.__binarycpt = 0

                         #time
                         l = line.find('ime="')
                         r = line.find('" ',l+1)
                         currentTime = float(line[l+5:r])

                         self.filePointer.seek(pos)
                         self._isOpen = True
                         self.__binfilecounter = binfilecounter
                         self.currentTime = currentTime
                         return
                raise(Exception("Unable Open file in append mode "))
            else:
                self.filePointer = open(self.fileName, 'w')

        except: # pragma: no cover
            print(TFormat.InRed("Error File Not Open"))
            raise

        self._isOpen = True

        self.filePointer.write('<?xml version="1.0" encoding="utf-8"?>\n')
        self.filePointer.write('<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.92">\n')
        self.filePointer.write('<Domain>\n')

        if self.__isTemporalOutput:
            self.filePointer.write('<Grid Name="Grid_T" GridType="Collection" CollectionType="Temporal"  >\n')
        #clean the binary file  =
        if self.isBinary():
            self.__binaryFilePointer = open (self.__binFileName, "wb")
        if self.__keepXmlFileInSaneState:
            self.WriteTail()

    def Close(self):
        if self.isOpen():
            self.WriteTail()
            self.filePointer.close()
            self._isOpen = False
            if self.isBinary():
                self.__binaryFilePointer.close()
            #print("File : '" + self.fileName + "' is close.")

    def WriteTail(self):
        if self.isOpen():
            filepos = self.filePointer.tell()

            if self.__isTemporalOutput:
                self.__WriteTime()
                self.filePointer.write('    </Grid><!--Temporal pos="'+str(filepos)+'" __binfilecounter="'+str(self.__binfilecounter)+'" time="'+str(self.currentTime)+'" -->\n')
            self.filePointer.write('  </Domain>\n')
            self.filePointer.write('</Xdmf>\n')
            # we put the pointer just before the tail so we can continue writting
            # to the file
            self.filePointer.seek(filepos)

    def __WriteGeoAndTopo(self,baseMeshObject):
        if baseMeshObject.IsConstantRectilinear() :
            origin = baseMeshObject.GetOrigin()
            spacing = baseMeshObject.GetSpacing()
            dims = baseMeshObject.GetDimensions() ## number of nodes per
            dimensionality = baseMeshObject.GetDimensionality()
            if dimensionality != 3:
                origin = np.append(origin,0)
                spacing = np.append(spacing,1)
                dims = np.append(dims,1)

            dimensionality = 3
            #if dimensionality == 3:
            self.filePointer.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
            #else:
            #    self.filePointer.write('    <Geometry Type="ORIGIN_DXDY">\n')


            self.filePointer.write('      <DataItem DataType="Float" Dimensions="'+str(dimensionality)+'" Format="XML" Precision="8">'+ArrayToString(reversed(origin)) +'</DataItem>\n')
            self.filePointer.write('      <DataItem DataType="Float" Dimensions="'+str(dimensionality)+'" Format="XML" Precision="8">'+ArrayToString(reversed(spacing)) +'</DataItem>\n')
            self.filePointer.write('    </Geometry>\n')
            self.filePointer.write('    <Topology Dimensions="'+ArrayToString(reversed(dims))  +'" Type="'+str(dimensionality)+'DCoRectMesh"/>\n')
        elif baseMeshObject.IsRectilinear() :# pragma: no cover
            print(TFormat.InRed("Mesh Type Rectilinear Not Supported"))        # pragma: no cover
            raise Exception                                                    # pragma: no cover
        elif baseMeshObject.IsStructured() :
            dimensionality = baseMeshObject.GetDimensionality()
            dims = baseMeshObject.GetDimensions() ## number of nodes per

            self.filePointer.write('    <Geometry Type="XYZ">\n')
            self.__WriteDataItem(baseMeshObject.GetPosOfNodes().ravel(), (baseMeshObject.GetNumberOfNodes(),3)  )
            self.filePointer.write('    </Geometry>\n')
            self.filePointer.write('    <Topology Dimensions="'+ArrayToString(reversed(dims))  +'" Type="'+str(dimensionality)+'DSMesh"/>\n')
        elif baseMeshObject.IsUnstructured() :
            self.filePointer.write('    <Geometry Type="XYZ">\n')
            if ( baseMeshObject.GetDimensionality()  == 2 ):
                nodes = baseMeshObject.GetPosOfNodes()
                nodes = np.concatenate((nodes,np.zeros((baseMeshObject.GetNumberOfNodes(),1))), axis=1 );
                self.__WriteDataItem(nodes.ravel(), (baseMeshObject.GetNumberOfNodes(),3)  )
            else:
                self.__WriteDataItem(baseMeshObject.GetPosOfNodes().ravel(), (baseMeshObject.GetNumberOfNodes(),3)  )

            self.filePointer.write('    </Geometry>\n')
            if len(baseMeshObject.elements) > 1:
                self.filePointer.write('    <Topology TopologyType="Mixed" NumberOfElements="{}">\n'.format(baseMeshObject.GetNumberOfElements()))
                ntotalentries = 0
                for ntype, data in baseMeshObject.elements.items():
                    ntotalentries += data.GetNumberOfElements()*(data.GetNumberOfNodesPerElement()+1)
                    if data.elementType == 'bar2' or data.elementType == 'point1':
                        ntotalentries += data.GetNumberOfElements()



                dataarray = np.empty((ntotalentries,),dtype=np.int)
                cpt =0;
                for ntype, data in baseMeshObject.elements.items():
                   self.PrintVerbose("printing elements of type : {}".format(data.elementType) )
                   self.PrintDebug("printing {}  elements".format(data.GetNumberOfElements()) )
                   self.PrintDebug("with  {}  nodes per element ".format(data.GetNumberOfNodesPerElement()) )
                   elemtype = XdmfNumber[ntype]
                   for i in range(data.GetNumberOfElements() ):
                       dataarray[cpt] = elemtype
                       cpt += 1;
                       if elemtype == 0x2 :
                           dataarray[cpt] = 2
                           cpt += 1;
                       elif elemtype == 0x1 :
                           dataarray[cpt] = 1
                           cpt += 1;

                       for j in range(data.GetNumberOfNodesPerElement()):
                           dataarray[cpt] = data.connectivity[i,j]
                           cpt += 1;
                self.PrintDebug("Number Of Entries {}".format(ntotalentries))
                self.PrintDebug("counter {}".format(cpt))

                self.__WriteDataItem(dataarray)
            elif len(baseMeshObject.elements):
                elements = list(baseMeshObject.elements.keys())[0]
                elementType = XdmfName[elements]
                self.filePointer.write('    <Topology TopologyType="{}" NumberOfElements="{}" '.format(elementType,baseMeshObject.GetNumberOfElements()))
                if XdmfNumber[elements] == 0x2:
                    self.filePointer.write('NodesPerElement="2"  ')
                if XdmfNumber[elements] == 0x1:
                    self.filePointer.write('NodesPerElement="1"  ')
                self.filePointer.write(' >\n')
                self.__WriteDataItem(baseMeshObject.elements[elements].connectivity.ravel())
            else:
                self.filePointer.write('    <Topology TopologyType="mixed" NumberOfElements="0"  >\n')

            self.filePointer.write('    </Topology> \n')

        else:                                                                  # pragma: no cover
            print(TFormat.InRed("Mesh Type Not Supported"))                    # pragma: no cover
            raise Exception                                                    # pragma: no cover

    def __WriteAttribute(self,data,name,center,baseMeshObject):
       shape = None
       if center == "Node":
           ndata = baseMeshObject.GetNumberOfNodes()
       elif center == "Cell":
           ndata = baseMeshObject.GetNumberOfElements()
       elif center == "Grid":
           ndata = 1;
       else:
          raise Exception('Cant treat this type of field support' + center )    # pragma: no cover

       if data.size == ndata:
           attype = "Scalar"
           #print(shape)
           #print(data.shape)
           if baseMeshObject.IsConstantRectilinear():
             if center == "Node":
                  shape = baseMeshObject.GetDimensions()
             elif center == "Cell":
                  shape = baseMeshObject.GetDimensions() -1
             else:
                  shape = [1]

             if len(shape)>1 :
               if len(data.shape) <= 2:
                   data.shape = tuple(shape)
                   data = data.T
               else:
                   data = data.T


           if baseMeshObject.IsConstantRectilinear():
             shape = np.array(shape)
             if center == "Node":
                  shape = shape.T
             elif center == "Cell":
                  shape = shape.T

       elif data.size == ndata*3:

           attype = "Vector"
           #print(attype)
           #print(data.shape)
           #print(shape)
           if baseMeshObject.IsConstantRectilinear() :
             if center == "Node":
                  shape = baseMeshObject.GetDimensions()
             elif center == "Cell":
                  shape = baseMeshObject.GetDimensions() -1
             else:
                  shape = [1]


             if len(shape)>1:
               if baseMeshObject.GetDimensionality() == 2:
                   shape = tuple(shape) + (1,3,)
               else:
                  shape = tuple(shape) + (3,)
               #print(shape)
               #if baseMeshObject.GetDimensionality() == 3:
               if len(data.shape) <= 2:
                  #shape = (shape[0], shape[1],shape[2],3)
                  data.shape = shape
               data = data.transpose(2,1,0,3)
               #else:
               #    if len(data.shape) <= 2:
               #        data.shape = shape
               #    data = data.transpose(1,0,2)


           if baseMeshObject.IsConstantRectilinear():
             shape = np.array(shape)
             if center == "Node":
                  shape = shape[[2,1,0,3]]
             elif center == "Cell":
                  shape = shape[[2,1,0,3]]

       elif data.size == ndata*2:
           if baseMeshObject.IsConstantRectilinear() :
             if center == "Node":
                  shape = baseMeshObject.GetDimensions()
             elif center == "Cell":
                  shape = baseMeshObject.GetDimensions() -1
             else:
                  shape = [1]

           # add an extra zeros to the data and print it out as Vector 3
           if baseMeshObject.GetDimensionality() == 2:
                 shape = tuple(shape) + (1,2,)
           else:
                 shape = tuple(shape) + (2,)



           if len(data.shape) <= 2:
               data.shape = shape

           #data = data.transpose(1,0,2)

           shape = list(data.shape)
           shape [-1] = 3
           shape = tuple(shape)
           data1 = np.zeros(shape ,dtype=data.dtype)

           data1[...,0:2] = data

           self.__WriteAttribute(data1.ravel(),name,center,baseMeshObject)

           return
       else:
           print(TFormat.InRed("I dont kow how to treak fields '"+ str(name)+"' with " +str(data.size/ndata) +" components"))  # pragma: no cover
           print(TFormat.InRed("Data has size : " + str(data.size) ))  # pragma: no cover
           print(TFormat.InRed("But support has size : " + str(ndata) ))  # pragma: no cover

           raise Exception                                                                                    # pragma: no cover
       #print(data.shape)
       #print(shape)
       self.filePointer.write('    <Attribute Center="'+center+'" Name="'+name+'" Type="'+attype+'">\n')#
       self.PrintDebug("Writing field '"+name +"' at '"+center+ "' of type " + attype )
       self.__WriteDataItem(data.ravel(),shape)

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

         if not self.isOpen() :
            if self.automaticOpen:
                self.Open()
            else:
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
         if self.__isTemporalOutput and self.__printTimeInsideEachGrid :
             self.filePointer.write('    <Time Value="'+str(self.currentTime)+'" /> \n')

         self.__WriteGeoAndTopo(baseMeshObject)

         #Writing tags for points
         #try:
         for tag in baseMeshObject.nodesTags:
                 name = "Tag_" + tag.name
                 data = np.zeros((baseMeshObject.GetNumberOfNodes(),1),dtype=np.int8)
                 data[baseMeshObject.nodesTags[tag.name].GetIds()] = 1;
                 self.__WriteAttribute(np.array(data), name, "Node",baseMeshObject)
         #except:
         #   pass

         #Cell Tags
         baseMeshObject.PrepareForOutput();

         celtags = baseMeshObject.GetNamesOfElemTags()
         for tagname in celtags:
             name = "Tag_" + tagname
             data = baseMeshObject.GetElementsInTag(tagname)
             res = np.zeros((baseMeshObject.GetNumberOfElements(),1),dtype=np.int8)
             res[data] = 1;

             self.__WriteAttribute(np.array(res), name, "Cell", baseMeshObject)

         for i in range(len(PointFields)):
           name = 'PField'+str(i)
           if len(PointFields)  == len(PointFieldsNames):
               name = PointFieldsNames[i]

           self.__WriteAttribute(np.array(PointFields[i]), name, "Node",baseMeshObject)

         for i in range(len(CellFields)):
           name = 'CField'+str(i)
           if len(CellFields) == len(CellFieldsNames):
               name = CellFieldsNames[i]

           self.__WriteAttribute(np.array(CellFields[i]), name, "Cell",baseMeshObject)


         for i in range(len(GridFields)):

           name = 'GField'+str(i)
           if len(GridFields) == len(GridFieldsNames):
               name = GridFieldsNames[i]

           self.__WriteAttribute(np.array(GridFields[i]), name, "Grid",baseMeshObject)

         self.filePointer.write('    </Grid>\n')

         if(self.__keepXmlFileInSaneState):
             self.WriteTail()
         if self.isBinary() :# pragma: no cover
           self.__binaryFilePointer.flush()


    def __WriteTime(self):
        """ this function is called by the WriteTail, this function must NOT change
         the state of the instance, also no writting to binary neither hdf5 files """
        if self.isOpen() and self.__printTimeInsideEachGrid == False:
            #self.filePointer.write('<Time TimeType="List">\n')
            #self.__WriteDataItem(self.timeSteps)
            #self.filePointer.write('</Time>\n')
            self.filePointer.write('<Time Value="'+ (" ".join(str(x) for x in self.timeSteps)) +'"/>\n')

    def __WriteDataItem(self,_data, _shape= None):

        import numpy as np
        data = np.array(_data)
        if _shape is None:
            _shape = _data.shape
        shape = np.array(_shape)

        if self.isOpen():
            if data.dtype == np.float64:
                typename = 'Float'
                s = data.dtype.itemsize
            elif data.dtype == np.float32:
                typename = 'Float'
                s = data.dtype.itemsize
            elif data.dtype == np.int32:
                typename = 'Int'
                s = data.dtype.itemsize
            elif data.dtype == np.int64:
                typename = 'Int'
                s = data.dtype.itemsize
            elif data.dtype == np.int8:
                typename = 'Char'
                s = data.dtype.itemsize
            else:
                print('Output Not implemented for data of type ' + str(type(data[0])))              # pragma: no cover
                raise                                                          # pragma: no cover

            dimension = ArrayToString(shape)

            # to test this feature a big file must be created (so we dont test it)
            if self.isBinary() and len(data) > self.__XmlSizeLimit:# pragma: no cover
                if self.__binarycpt > (2**30) :
                    self.__binaryFilePointer.close()
                    self.NewBinaryFilename()
                    self.__binaryFilePointer = open (self.__binFileName, "wb")
                    self.__binarycpt = 0

                #bindata = bytearray(data.ravel())

                #self.__binaryFilePointer.write(bindata)
                data.ravel().tofile(self.__binaryFilePointer)

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
                self.filePointer.write(" ".join(str(x) for x in data.ravel()))


            self.filePointer.write('</DataItem>\n')

def WriteTest(tempdir,Temporal, Binary):

    from BasicTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,3,4]);
    myMesh.SetSpacing([0.1, 0.1, 0.1]);
    myMesh.SetOrigin([-2.5,-1.2,-1.5]);

    dataT = np.arange(24,dtype=np.float32)
    dataT.shape = (2,3,4)
    dataDep = np.arange(24*3)+0.3

    dataDep.shape = (2,3,4,3)

    writer = XdmfWriter(tempdir + 'TestOutput_Bin_'+str(Binary)+'_Temp_'+str(Temporal)+'.xmf')
    writer.SetTemporal(Temporal)
    writer.SetBinary(Binary)
    writer.Open()
    print(writer)
    writer.Write(myMesh,PointFields=[dataT, dataDep], PointFieldsNames=["Temp","Dep"],CellFields=[np.arange(6)],CellFieldsNames=['S'], Time=0);
    writer.Write(myMesh,GridFields=[0, 1], GridFieldsNames=['K','P'], TimeStep = 1);
    writer.Close()

def WriteTestAppend(tempdir,Temporal, Binary):

    from BasicTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,3,4]);
    myMesh.SetSpacing([0.1, 0.1, 0.1]);
    myMesh.SetOrigin([-2.5,-1.2,-1.5]);

    dataT = np.arange(24,dtype=np.float32)
    dataT.shape = (2,3,4)
    dataDep = np.arange(24*3)+0.3

    dataDep.shape = (2,3,4,3)

    writer = XdmfWriter(tempdir + 'TestOutput_Bin_'+str(Binary)+'_Temp_'+str(Temporal)+'.xmf')
    writer.SetAppendMode(True)
    writer.SetTemporal(Temporal)
    writer.SetBinary(Binary)
    writer.Open()
    print(writer)
    writer.Write(myMesh,PointFields=[dataT, dataDep], PointFieldsNames=["Temp","Dep"],CellFields=[np.arange(6)],CellFieldsNames=['S']);
    writer.Write(myMesh,GridFields=[0, 1], GridFieldsNames=['K','P'], TimeStep = 1);
    writer.Close()

def CheckIntegrity():
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh
    import BasicTools.FE.UnstructuredMesh as UM
    from BasicTools.FE.UnstructuredMeshTools import CreateMeshOfTriangles

    TestTempDir.SetTempPath("/tmp/BasicTools_Test_Directory_jMnw6z_safe_to_delete/")
    tempdir = TestTempDir.GetTempPath()

    res = CreateMeshOfTriangles([[0.,0.,0],[1.,2.,3],[1, 3, 2]], np.array([[0,1,2]]))
    print(res)
    res.SetGlobalDebugMode()

    WriteMeshToXdmf(tempdir+"TestUnstructured.xdmf", res, PointFields = [np.array([1.,2,3])], CellFields =[ np.array([1])] ,GridFields= [[0],  np.array([1,2,3]).astype(np.int64) ],
                                                                    PointFieldsNames = ["PS"],
                                                                    CellFieldsNames = ["CS"],
                                                                    GridFieldsNames = ["GS", "GV"] , Binary= True)

    elements = res.GetElementsOfType(EN.Bar_2)
    elements.AddNewElement([1,2],1)
    elements = res.GetElementsOfType(EN.Point_1)
    elements.AddNewElement([0],2)
    res.nodes = res.nodes[:,0:2]
    res.GetNodalTag("First_Point").AddToTag(0)

    res.AddElementToTagUsingOriginalId(1,"bars")

    WriteMeshToXdmf(tempdir+"TestUnstructured_multielements.xdmf", res, PointFields = [np.array([1.,2,3])], CellFields =[ np.array([1,0,4])] ,GridFields= [[0]],
                                                                    PointFieldsNames = ["PS"],
                                                                    CellFieldsNames = ["CS"],
                                                                    GridFieldsNames = ["GS"] , Binary= True)
    #----------------------
    res = UM.UnstructuredMesh()
    WriteMeshToXdmf(tempdir+"TestUnstructured_EmptyMesh.xdmf", res)

    res.nodes = np.array([[0,0,0],[1,0,0]],dtype=np.float32)
    elements = res.GetElementsOfType(EN.Point_1)
    elements.AddNewElement([0],1)
    WriteMeshToXdmf(tempdir+"TestUnstructured_OnlyPoints.xdmf", res)

    res = UM.UnstructuredMesh()
    res.nodes = np.array([[0,0,0],[1,0,0]],dtype=np.float32)
    elements = res.GetElementsOfType(EN.Bar_2)
    elements.AddNewElement([0,1],1)
    WriteMeshToXdmf(tempdir+"TestUnstructured_OnlyBars.xdmf", res)



    #----------------------

    WriteTest(tempdir,False, False)
    WriteTest(tempdir,False, True)
    WriteTest(tempdir,True, False)
    WriteTest(tempdir,True, True)
    WriteTestAppend(tempdir,True, True)



    WriteMeshToXdmf(tempdir+'testdirect.xdmf',ConstantRectilinearMesh() )

    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetFileName(tempdir+'testerros.xdmf')
    writer.SetXmlSizeLimit(0)
    writer.Open();

    ## test of the errors
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

    # no error anymore just a warning
    #try:
    writer.Open();
    #except:
    #    pass

    writer.Close()

    try:
        writer.Write(ConstantRectilinearMesh())
        return "Not   ok"# pragma: no cover
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

    print("ConstantRectilinearMesh in 2D")
    CRM2D = ConstantRectilinearMesh(2);

    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetXmlSizeLimit(0)
    writer.SetBinary(True)
    writer.Open(filename=tempdir+'testdirect.xdmf');
    writer.Write(CRM2D, PointFields = [ np.zeros((CRM2D.GetNumberOfNodes(),) ).astype(np.float32), np.zeros((CRM2D.GetNumberOfNodes(),) ).astype(np.int) ],
                                                            PointFieldsNames = ["Test", "testint"],
                                                            CellFields= [ np.arange(CRM2D.GetNumberOfElements()*2 ).astype(np.float64)],
                                                            CellFieldsNames = [ "TestV"] );
    writer.Close()



    from BasicTools.IO.XdmfReader import XdmfReader  as XR
    f = XR(filename = tempdir+'testdirect.xdmf' )
    f.lazy = False;
    f.Read();
    print(tempdir+'testdirect.xdmf')
    f.xdmf.GetDomain(0).GetGrid(0).GetFieldsOfType("Cell")
    print(f.xdmf.GetDomain(0).GetGrid(0).attributes )

    #print('testdirect.xdmf')
    #os.system('cmd /C C:\Users\D584808\Apps\ParaView-5.0.1-Qt4-OpenGL2-Windows-64bit\\bin\paraview.exe ' + 'testdirect.xdmf')

    print("Structured Mesh in 3D")
    import BasicTools.FE.StructuredMesh as SM

    SM3D = SM.StructuredMesh();


    WriteMeshToXdmf(tempdir+'StructuredMesh.xdmf',SM3D, PointFields = [ np.arange(SM3D.GetNumberOfNodes()) ],
                                                                       PointFieldsNames = ["Test"],
                                                                       CellFields= [ np.arange(SM3D.GetNumberOfElements()*3)],
                                                                       CellFieldsNames = [ "TestV"] );

    from BasicTools.IO.XdmfReader import XdmfReader  as XR
    f = XR(filename = tempdir+'StructuredMesh.xdmf' )
    f.lazy = False;
    f.Read();

    domain = f.xdmf.GetDomain(0)
    grid  = domain.GetGrid(0)
    grid.GetFieldsOfType("Cell")
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover

    #import BasicTools.IO.AscReader as AR
    #m =  AR.ReadAsc(fileName='C:\\Users\\D584808\\Documents\\Projects\\Python\\Topotools\\SUPPORT_VERIN_DATA1.ASC')
    #WriteMeshToXdmf('FromASCReader.xdmf',m ,Binary= False)

    #m.ComputeBoundingBox()
    #print(m.boundingMin)
    #print(m.boundingMax)
