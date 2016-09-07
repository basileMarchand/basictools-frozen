# -*- coding: utf-8 -*-

from OTTools.Helpers.TextFormatHelper import TFormat
import numpy as np
import os
import OTTools.FE.ElementNames as EN

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
XdmfName[EN.Triangle_6] = 'Triagle_6'
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
XdmfNumber[EN.Triangle_6] = 0x23
XdmfNumber[EN.Quadrangle_9] = 0x24
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
        self.__binfilecounter = 0
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
            return
            #raise Exception


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
            if dimensionality == 3:
                self.filePointer.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
            else:
                self.filePointer.write('    <Geometry Type="ORIGIN_DXDY">\n')


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
            self.__WriteDataItem(baseMeshObject.GetPosOfNodes().flatten(), (baseMeshObject.GetNumberOfNodes(),3)  )
            self.filePointer.write('    </Geometry>\n')
            self.filePointer.write('    <Topology Dimensions="'+ArrayToString(reversed(dims))  +'" Type="'+str(dimensionality)+'DSMesh"/>\n')
        elif baseMeshObject.IsUnstructured() :
            self.filePointer.write('    <Geometry Type="XYZ">\n')
            self.__WriteDataItem(baseMeshObject.GetPosOfNodes().flatten(), (baseMeshObject.GetNumberOfNodes(),3)  )
            self.filePointer.write('    </Geometry>\n')
            if len(baseMeshObject.elements) > 1:
                self.filePointer.write('    <Topology TopologyType="Mixed" NumberOfElements="{}">\n'.format(baseMeshObject.GetNumberOfElements()))
                ntotalentries = 0
                for ntype, data in baseMeshObject.elements.iteritems():
                    ntotalentries += data.GetNumberOfElements()*(data.GetNumberOfNodesPerElement()+1)
                    if data.elementType == 'bar2' or data.elementType == 'point1':
                        ntotalentries += data.GetNumberOfElements()



                dataarray = np.empty((ntotalentries,),dtype=np.int)
                cpt =0;
                for ntype, data in baseMeshObject.elements.iteritems():
                   print("printing {}  elements".format(data.elementType) )
                   print("printing {}  elements".format(data.GetNumberOfElements()) )
                   print("with  {}  nodes per element ".format(data.GetNumberOfNodesPerElement()) )
                   elemtype = XdmfNumber[ntype]
                   for i in xrange(data.GetNumberOfElements() ):
                       dataarray[cpt] = elemtype
                       cpt += 1;
                       if elemtype == 0x2 :
                           dataarray[cpt] = 2
                           cpt += 1;
                       elif elemtype == 0x1 :
                           dataarray[cpt] = 1
                           cpt += 1;

                       for j in xrange(data.GetNumberOfNodesPerElement()):
                           dataarray[cpt] = data.connectivity[i,j]
                           cpt += 1;
                print("Number Of Entries {}".format(ntotalentries))
                print("counter {}".format(cpt))

                self.__WriteDataItem(dataarray)
            else:
                #
                elementType = XdmfName[baseMeshObject.elements.keys()[0]]
                self.filePointer.write('    <Topology TopologyType="{}" NumberOfElements="{}"  >\n'.format(elementType,baseMeshObject.GetNumberOfElements()))
                self.__WriteDataItem(baseMeshObject.elements[baseMeshObject.elements.keys()[0]].connectivity.flatten())

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
          raise Exception('Cant treat this type of field support' + center )

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

               shape = tuple(baseMeshObject.GetDimensions()) + (3,)
               #print(shape)
               if baseMeshObject.GetDimensionality() == 3:
                   if len(data.shape) <= 2:
                       #shape = (shape[0], shape[1],shape[2],3)
                       data.shape = shape
                   data = data.transpose(2,1,0,3)
               else:
                   if len(data.shape) <= 2:
                       data.shape = shape
                   data = data.transpose(1,0,2)


       elif data.size == ndata*2:
           # add an extra zeros to the data and print it out as Vector 3
           shape = tuple(baseMeshObject.GetDimensions())+ (2,)

           if len(data.shape) <= 2:
               data.shape = shape

           #data = data.transpose(1,0,2)

           shape = list(data.shape)
           shape [-1] = 3
           shape = tuple(shape)
           data1 = np.zeros(shape ,dtype=data.dtype)

           data1[...,0:2] = data

           self.__WriteAttribute(data1.flatten(),name,center,baseMeshObject)

           return
       else:
           print(TFormat.InRed("I dont kow how to treak fields with " +str(data.size/ndata) +" components"))  # pragma: no cover
           raise Exception                                                                                    # pragma: no cover

       self.filePointer.write('    <Attribute Center="'+center+'" Name="'+name+'" Type="'+attype+'">\n')#


       self.__WriteDataItem(data.flatten(),shape)

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

         #Writing tags for points
         #try:
         for tagname in baseMeshObject.nodesTags:
                 name = "Tag_" + tagname
                 data = np.zeros((baseMeshObject.GetNumberOfNodes(),1),dtype=np.int)
                 data[baseMeshObject.nodesTags[tagname].id] = 1;
                 self.__WriteAttribute(np.array(data), name, "Node", (baseMeshObject.GetNumberOfNodes(),),baseMeshObject)
         #except:
         #   pass

         #Cell Tags
         baseMeshObject.PrepareForOutput();

         celtags = baseMeshObject.GetNamesOfCellTags()
         for tagname in celtags:
             name = "Tag_" + tagname
             #data = baseMeshObject.GetElementTag(tagname)
             data = baseMeshObject.GetElementsInTag(tagname)
             res = np.zeros((baseMeshObject.GetNumberOfElements(),1),dtype=np.int)
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

           self.__WriteAttribute(np.array(GridFields[i]), name, "Grid",baseMeshObject)

         self.filePointer.write('    </Grid>\n')

    def __WriteTime(self):

        if self.__isOpen:
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
                if self.__binarycpt > (2**30) :
                    self.__binaryFilePointer.close()
                    self.NewBinaryFilename()
                    self.__binaryFilePointer = open (self.__binFileName, "wb")
                    self.__binarycpt = 0

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

    writer = XdmfWriter(tempdir + 'TestOutput_Bin_'+str(Binary)+'_Temp_'+str(Temporal)+'.xmf')
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
    import OTTools.FE.UnstructuredMesh as UM

    tempdir = TestTempDir.GetTempPath()

    res = UM.CreateMeshOfTetras([[0.,0.,0],[1.,2.,3],[1, 3, 2]], np.array([[0,1,2]]))
    print(res)
    WriteMeshToXdmf(tempdir+"TestUnstructured.xdmf", res, PointFields = [np.array([1.,2,3])], CellFields =[ np.array([1])] ,GridFields= [[0]],
                                                                    PointFieldsNames = ["PS"],
                                                                    CellFieldsNames = ["CS"],
                                                                    GridFieldsNames = ["GS"] , Binary= True)

    #----------------------

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

    print("ConstantRectilinearMesh in 2D")
    CRM2D = ConstantRectilinearMesh(2);
    WriteMeshToXdmf('testdirect.xdmf',CRM2D, PointFields = [ np.zeros((CRM2D.GetNumberOfNodes(),) ) ], PointFieldsNames = ["Test"]);

    #print('testdirect.xdmf')
    #os.system('cmd /C C:\Users\D584808\Apps\ParaView-5.0.1-Qt4-OpenGL2-Windows-64bit\\bin\paraview.exe ' + 'testdirect.xdmf')

    print("Structured Mesh in 3D")
    import OTTools.FE.StructuredMesh as SM

    SM3D = SM.StructuredMesh();


    WriteMeshToXdmf(tempdir+'StructuredMesh.xdmf',SM3D, PointFields = [ np.arange(SM3D.GetNumberOfNodes()) ], PointFieldsNames = ["Test"]);
    print(tempdir)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover

    #import OTTools.IO.AscReader as AR
    #m =  AR.ReadAsc(fileName='C:\\Users\\D584808\\Documents\\Projects\\Python\\Topotools\\SUPPORT_VERIN_DATA1.ASC')
    #WriteMeshToXdmf('FromASCReader.xdmf',m ,Binary= False)

    #m.ComputeBoundingBox()
    #print(m.boundingMin)
    #print(m.boundingMax)
