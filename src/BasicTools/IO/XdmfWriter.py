# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
import os

from BasicTools.Helpers.TextFormatHelper import TFormat

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.XdmfTools import XdmfName,XdmfNumber
from BasicTools.IO.WriterBase import WriterBase as WriterBase
from BasicTools.Helpers.MPIInterface import MPIInterface as MPI

def ArrayToString(data):
    return " ".join(str(x) for x in data)


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


def WriteMeshToXdmf(filename,
                    baseMeshObject,
                    PointFields = None,
                    CellFields = None,
                    GridFields= None,
                    IntegrationPointData=None,
                    PointFieldsNames = None,
                    CellFieldsNames=None,
                    GridFieldsNames=None,
                    IntegrationPointDataNames=None,
                    IntegrationRule=None,
                    Binary= True, ):
    """
    Functional version of the Xdmf Writer
    """

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

    if IntegrationPointData is None:
        IntegrationPointData  = [];

    if IntegrationPointDataNames is None:
        IntegrationPointDataNames  = [];

    writer = XdmfWriter(filename)
    writer.SetBinary(Binary)
    writer.Open()
    writer.Write(baseMeshObject,
                 PointFields= PointFields,
                 CellFields = CellFields,
                 GridFields = GridFields,
                 IntegrationPointData=IntegrationPointData,
                 PointFieldsNames = PointFieldsNames,
                 CellFieldsNames = CellFieldsNames,
                 GridFieldsNames = GridFieldsNames,
                 IntegrationPointDataNames=IntegrationPointDataNames,
                 IntegrationRule=IntegrationRule,
                 )
    writer.Close()
#

class inmemoryfile():
   """
   Helper class to write the xmf part of the file into memory
   """
   def __init__(self,saveFilePointer):
       self.data = ""
       self.saveFilePointer = saveFilePointer
   def write(self,data):
       self.data += data
   def tell(self):
       return 0
   def seek(self,pos):
       return 0
   def close(self):
       pass

class BinaryStorage(object):

    def __init__(self,data=None, filePointer=None):
        self.filename = ""
        self.offset = 0
        self.type = None
        self.itemsize = 0
        self.vectorsize = 0
        self.usedByNInstances = 0
        if data is not None:
            self.usedByNInstances += 1
            self.type = data.dtype
            self.itemsize = data.dtype.itemsize
            self.vectorsize = data.size

        if filePointer is not None:
            self.filename = filePointer.name
            self.offset = filePointer.tell()

    def __disp__(self):
        return str(self.vectorsize) +":"+str(self.usedByNInstances)

    def ChangePathOfBinaryStorage(self,newpath):
        import os
        self.filename = newpath + os.sep + os.path.basename(self.filename)

    def GetData(self):
       f = open(self.filename,'rb')
       f.seek(self.offset,0)
       data = np.fromfile(f,self.type,self.vectorsize,sep="")
       f.close()
       return data

    def UpdateHeavyStorage(self,data):
        if self.usedByNInstances > 1:
            raise(Exception("This pointer is used for more than 1 field please (overright or setup the writer with the option maxStorageSize=0"))

        if data.size != self.vectorsize:
            raise(Exception("Size of data and storage not compatible"))
        f = open(self.filename,'r+b')
        f.seek(self.offset,0)
        data.astype(self.type).ravel().tofile(f)
        f.close()




class XdmfWriter(WriterBase):
    """
    Class to Write Xdmf files for:
        - classic finit element solutions
        - domain decomposition problem (multi mesh)
        - transient solution (the mesh changes in time)
        - solution written in parafac format (monolitic or in ddm mpi)
    """


    def __init__(self, fileName = None):
        super(XdmfWriter,self).__init__()
        self.canHandleTemporal = True
        self.canHandleAppend = True
        self.canHandleMultidomain = True

        self.fileName = None;
        self.timeSteps = [];
        self.parafacCpt = 0;
        self.ddmCpt = 0;
        self.currentTime = 0;
        self.__XmlSizeLimit = 0
        self.automaticOpen = False;

        self.SetBinary(False)
        self.__binFileName = None;
        self.__filePointer = None;
        #self.__isOpen = False;
        self.__binarycpt = 0;
        self.__binfilecounter = 0
        self.__keepXmlFileInSaneState = True;
        self.__isParafacFormat = False

        #set to off is you what to put the time at the end of the temporal grid
        #keep this option True to be compatible with the XMDF3 reader of Paraview
        self.__printTimeInsideEachGrid = True

        self.SetFileName(fileName)
        self.pointfieldsStorage = {}
        self.cellfieldsStorage = {}
        self.gridfieldsStorage = {}
        self.iptorage = {}
        self.globalStorage = {}
        self.maxStorageSize = 50

    def __str__(self):
        res  = 'XdmfWriter : \n'
        res += '   FileName : '+ str(self.fileName) +'\n'
        res += '   isParafacFormat : '+ ('True' if self.__isParafacFormat else "False") +'\n'
        if self.isBinary():
           res += '   Binary output Active \n'
           res += '   Binary FileName : '+ self.__binFileName +'\n'
        if self.IsTemporalOutput():
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

        self.fileName = fileName;
        self.__path  = os.path.abspath(os.path.dirname(fileName));
        self.binfilecounter = 0
        self.NewBinaryFilename()

    def SetParafac(self,val = True):
        self.__isParafacFormat = val

    def IsParafacOutput(self):
        return self.__isParafacFormat

    def NewBinaryFilename(self):
        name = os.path.splitext(os.path.abspath(self.fileName))[0]
        name += "" +str(self.__binfilecounter)
        if MPI.IsParallel():
            name += "D"+ str(MPI.Rank())
        name += ".bin"

        self.__binFileName = name
        self.__binFileNameOnly = os.path.basename(self.__binFileName)
        self.__binfilecounter +=1

    def Step(self, dt = 1):
        self.currentTime += dt;
        self.timeSteps.append(self.currentTime);


    def SetXmlSizeLimit(self,val):
        self.__XmlSizeLimit= val

    def Open(self, filename = None):

        # we dont use the open from WriterBase because the set binary is used
        # for the .bin file and not for the .xdmf file

        if self.isOpen() :
            print(TFormat.InRed("The file is already open !!!!!"))
            return

        if filename is not None:
            self.SetFileName(filename)

        ## we use unbuffered so we can repaire broken files easily
        try :
            # in python 3 we cant use unbuffered  text I/O (bug???)
            #self.filePointer = open(self.fileName, 'w',0)
            if self.InAppendMode():
                self.filePointer = open(self.fileName, 'r+')
                self.filePointer.seek(-100,2)

                lines = self.filePointer.readlines()
                for line in lines:
                    if line.find("<!--Temporal") != -1:
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
               mpi = MPI()
               if mpi.IsParallel():
                   self.__keepXmlFileInSaneState = False
                   if mpi.rank > 0:
                       self.filePointer = inmemoryfile(self.filePointer)
                   else:
                       self.filePointer = open(self.fileName, 'w')
               else:
                self.filePointer = open(self.fileName, 'w')

        except: # pragma: no cover
            print(TFormat.InRed("Error File Not Open"))
            raise

        self._isOpen = True

        self.filePointer.write('<?xml version="1.0" encoding="utf-8"?>\n')
        self.filePointer.write('<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.92">\n')
        self.filePointer.write('<Domain>\n')

        if self.IsTemporalOutput():
            self.filePointer.write('<Grid Name="Grid_T" GridType="Collection" CollectionType="Temporal"  >\n')
        if self.IsMultidomainOutput():
            self.filePointer.write('<Grid Name="Grid_S" GridType="Collection" CollectionType="Spatial" >\n')

        ## here we recover the output for the slave nodes (rank> 0) and then we
        ## write this information in the master node (rank = 0)
        self._OpenParallelWriter()

        if self.IsParafacOutput():
            self.filePointer.write('<Grid Name="Grid_P" GridType="Collection" CollectionType="None" >\n')

        if self.isBinary():
            self.__binaryFilePointer = open (self.__binFileName, "wb")

        if self.__keepXmlFileInSaneState:
            self.WriteTail()

    def _OpenParallelWriter(self):

         mpi = MPI()
         if mpi.mpiOK:
           self.__keepXmlFileInSaneState = False
           if mpi.rank > 0:
             self.filePointer = inmemoryfile(self.filePointer )

    def _CloseParallelWriter(self):

        mpi = MPI()
        if mpi.mpiOK and mpi.rank > 0:
             mpi.comm.send(self.filePointer.data, dest=0)
             self.filePointer = self.filePointer.saveFilePointer
        else:
             for i in range(1,mpi.size):
                 data = mpi.comm.recv( source=i)
                 self.filePointer.write(data)

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

            if self.IsParafacOutput():
                self.filePointer.write('</Grid> <!-- Parafac grid -->\n')

            self._CloseParallelWriter()

            if self.IsMultidomainOutput() :
                self.filePointer.write('</Grid> <!-- collection grid -->\n')



            if self.IsTemporalOutput():
                self.__WriteTime()
                self.filePointer.write('    </Grid><!--Temporal pos="'+str(filepos)+'" __binfilecounter="'+str(self.__binfilecounter)+'" time="'+str(self.currentTime)+'" -->\n')
            self.filePointer.write('  </Domain>\n')
            self.filePointer.write('</Xdmf>\n')
            # we put the pointer just before the tail so we can continue writting
            # to the file for a new time step
            self.filePointer.seek(filepos)

    def __WriteGeoAndTopo(self,baseMeshObject,name=None):

        if self.__isParafacFormat:
            if "ParafacDims" in baseMeshObject.props:
                self.filePointer.write('    <Information Name="Dims" Value="'+str(baseMeshObject.props["ParafacDims"])+'" /> \n')
                for i in range(baseMeshObject.props["ParafacDims"]):
                    self.filePointer.write('    <Information Name="Dim'+str(i)+'" Value="'+baseMeshObject.props["ParafacDim"+str(i)]+'" /> \n')
                if   "ParafacUnit0" in baseMeshObject.props:
                    for i in range(baseMeshObject.props["ParafacDims"]):
                        self.filePointer.write('    <Information Name="Unit'+str(i)+'" Value="'+baseMeshObject.props["ParafacUnit"+str(i)]+'" /> \n')

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
            self.__WriteDataItem(baseMeshObject.GetPosOfNodes().ravel(), (baseMeshObject.GetNumberOfNodes(),3) , name="GEO_S_"+str(name) )
            self.filePointer.write('    </Geometry>\n')
            self.filePointer.write('    <Topology Dimensions="'+ArrayToString(reversed(dims))  +'" Type="'+str(dimensionality)+'DSMesh"/>\n')
        elif baseMeshObject.IsUnstructured() :
            self.filePointer.write('    <Geometry Type="XYZ">\n')
            if ( baseMeshObject.GetDimensionality()  == 2 ):
                nodes = baseMeshObject.GetPosOfNodes()
                nodes = np.concatenate((nodes,np.zeros((baseMeshObject.GetNumberOfNodes(),1))), axis=1 );
                self.__WriteDataItem(nodes.ravel(), (baseMeshObject.GetNumberOfNodes(),3)  , name="GEO_U_"+str(name) )
            else:
                self.__WriteDataItem(baseMeshObject.GetPosOfNodes().ravel(), (baseMeshObject.GetNumberOfNodes(),3)  , name="GEO_U_"+str(name) )

            self.filePointer.write('    </Geometry>\n')
            if len(baseMeshObject.elements) > 1:
                self.filePointer.write('    <Topology TopologyType="Mixed" NumberOfElements="{0}">\n'.format(baseMeshObject.GetNumberOfElements()))
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

                self.__WriteDataItem(dataarray, name="Topo_U_"+str(name) )
            elif len(baseMeshObject.elements):
                elements = list(baseMeshObject.elements.keys())[0]
                elementType = XdmfName[elements]
                self.filePointer.write('    <Topology TopologyType="{}" NumberOfElements="{}" '.format(elementType,baseMeshObject.GetNumberOfElements()))
                if XdmfNumber[elements] == 0x2:
                    self.filePointer.write('NodesPerElement="2"  ')
                if XdmfNumber[elements] == 0x1:
                    self.filePointer.write('NodesPerElement="1"  ')
                self.filePointer.write(' >\n')
                self.__WriteDataItem(baseMeshObject.elements[elements].connectivity.ravel(), name="Topo_U_"+str(name) )
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



       elif data.size == ndata*9:

           attype = "Vector"
           if baseMeshObject.IsConstantRectilinear() :
             if center == "Node":
                  shape = baseMeshObject.GetDimensions()
             elif center == "Cell":
                  shape = baseMeshObject.GetDimensions()-1
             else:
                  shape = [1]
           #print(attype)
           #print(data.shape)
           #print(shape)

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
       self.filePointer.write('    <Attribute Center="'+center+'" Name="'+name+'" Type="'+attype+'">\n')#
       #self.PrintDebug("Writing field '"+name +"' at '"+center+ "' of type " + attype )
       try:
           self.__WriteDataItem(data.ravel(),shape,name=name)
       except:
           print("Error Writing heavy data of field: " + str(name))
           raise


       self.filePointer.write('    </Attribute>\n')

    def WriteIntegrationsPoints(self,datas):

        for elemName, data in datas.items():
            points, weight = data
            npdata = np.asarray(points,dtype=float).copy()
            nip = npdata.shape[0]
            if npdata.shape[1] < 3:
                b = np.zeros((nip,3))
                b[:,:-1] = npdata
                npdata = b

            self.filePointer.write('    <Information Name="QP" Value="')#
            self.filePointer.write(str(EN.numberOfNodes[elemName]) + " ")
            #from BasicTools.Containers.vtkBridge import vtkNumberByElementName
            #self.filePointer.write(str(vtkNumberByElementName[elemName]) + " ")
            self.filePointer.write(str(XdmfNumber[elemName]) + " ")


            self.filePointer.write(str(nip) + '" > \n')

            self.__WriteDataItem(npdata.ravel(),[npdata.size],name="ip")
            self.filePointer.write('</Information> \n')#

    def WriteIntegrationsPointDatas(self,names,datas):
        for i, name in enumerate(names):
            data = datas[i]
            self.filePointer.write('    <Information Name="IPF" Value="'+str(name)+'" > \n')#
            self.iptorage[name] = self.__WriteDataItem(data.ravel(),[data.size],name='IPD_'+name)
            self.filePointer.write('</Information> \n')#

    def NextDomain(self):
        self.parafacCpt = 0
        self.ddmCpt += 1
        if self.IsMultidomainOutput() :
            if self.IsParafacOutput():
                 self.filePointer.write('</Grid> <!-- Parafac grid -->\n')
                 self.filePointer.write('<Grid Name="Grid_P" GridType="Collection" CollectionType="None" >\n')
        else:
            raise(Exception("Cant make a new domain without a Multidomain output"))

    def MakeStep(self,Time=None,TimeStep=None):

        if self.IsMultidomainOutput():
            self.filePointer.write('</Grid> <!-- collection grid -->\n')
            self.filePointer.write('<Grid Name="Collection" GridType="Collection" CollectionType="Spatial" >\n')

        dt = 1
        if Time is not None:
           dt = Time - self.currentTime
        elif TimeStep is not None:
           dt = TimeStep

        if self.IsTemporalOutput():
           self.Step(dt)

        if self.IsTemporalOutput() and self.__printTimeInsideEachGrid and (self.IsMultidomainOutput() == False) :
             self.filePointer.write('    <Time Value="'+str(self.currentTime)+'" /> \n')

    def Write(self,baseMeshObject, PointFields = None, CellFields = None, GridFields= None, PointFieldsNames = None, CellFieldsNames= None, GridFieldsNames=None , Time= None, TimeStep = None,domainName=None, IntegrationRule=None, IntegrationPointData=None, IntegrationPointDataNames=None ):

         if (Time is not None or TimeStep is not None) and self.IsMultidomainOutput():
             raise(Exception("set time using MakeStep, not the Write option") )

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

         if IntegrationRule is None:
            IntegrationRule  = {};

         if IntegrationPointData is None:
            IntegrationPointData  = [];

         if IntegrationPointDataNames is None:
            IntegrationPointDataNames  = [];

         self.pointfieldsStorage = {}
         self.cellfieldsStorage = {}
         self.gridfieldsStorage = {}
         self.iptorage = {}

         if not self.isOpen() :
            if self.automaticOpen:
                self.Open()
            else:
                print(TFormat.InRed("Please Open The writer First!!!"))
                raise Exception

         if not self.IsMultidomainOutput():
             dt = 1
             if Time is not None:
                 dt = Time - self.currentTime
             elif TimeStep is not None:
                 dt = TimeStep
             if self.IsTemporalOutput():
                 self.Step(dt)

         if self.IsParafacOutput ():
             sufix = "P" + str(self.parafacCpt)
             self.parafacCpt += 1
         else:
             sufix = str(len(self.timeSteps))

         if self.IsMultidomainOutput():
             if MPI.IsParallel():
               sufix += "_D" + str(MPI.Rank())
             else:
               sufix += "_D" + str(self.ddmCpt)

         #if we have a constantRectilinear mesh with more than only "bulk"
         # elements, we add a collection to add all the rest of the elements

         if baseMeshObject.IsConstantRectilinear() and len(baseMeshObject.elements) > 1:
             self.filePointer.write('<Grid Name="Grid_S'+sufix+'" GridType="Collection" CollectionType="Spatial" >\n')
             if self.IsTemporalOutput() and self.__printTimeInsideEachGrid and (self.IsMultidomainOutput() == False) :
                self.filePointer.write('    <Time Value="'+str(self.currentTime)+'" /> \n')
             self.filePointer.write('    <Grid Name="Bulk">\n')
         else:
             self.filePointer.write('    <Grid Name="Grid_'+sufix+'">\n')
             if self.IsTemporalOutput() and self.__printTimeInsideEachGrid and (self.IsMultidomainOutput() == False) :
                 self.filePointer.write('    <Time Value="'+str(self.currentTime)+'" /> \n')

         self.__WriteGeoAndTopo(baseMeshObject,name=sufix)
         self.__WriteNodesTagsElementsTags(baseMeshObject,PointFieldsNames,CellFieldsNames)
         self.__WriteNodesFieldsElementsFieldsGridFields(baseMeshObject,
                                                   PointFieldsNames,PointFields,
                                                   CellFieldsNames,CellFields,
                                                   GridFieldsNames,GridFields)

         self.WriteIntegrationsPoints(IntegrationRule)
         self.WriteIntegrationsPointDatas(IntegrationPointDataNames,IntegrationPointData)

         if baseMeshObject.IsConstantRectilinear() and len(baseMeshObject.elements) > 1:
             self.filePointer.write('    </Grid>\n')
             from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
             tempmesh = UnstructuredMesh()
             tempmesh.nodes = baseMeshObject.nodes
             for name,data in baseMeshObject.elements.items():
                 if data.mutable :
                     tempmesh.elements[name] =  data
             self.filePointer.write('    <Grid Name="Sets">\n')
             self.__WriteGeoAndTopo(tempmesh)
             self.__WriteNodesTagsElementsTags(tempmesh,[],[])
             self.filePointer.write('    </Grid>\n')

         self.filePointer.write('    </Grid>\n')

         if(self.__keepXmlFileInSaneState):
             self.WriteTail()

         self.filePointer.flush()

         if self.isBinary() :# pragma: no cover
           self.__binaryFilePointer.flush()
    def __WriteNodesFieldsElementsFieldsGridFields(self,baseMeshObject,
                                                   PointFieldsNames,PointFields,
                                                   CellFieldsNames,CellFields,
                                                   GridFieldsNames,GridFields):

         for i in range(len(PointFields)):
           name = 'PField'+str(i)
           if len(PointFields)  == len(PointFieldsNames):
               name = PointFieldsNames[i]
           self.pointfieldsStorage[name] = self.__WriteAttribute(np.array(PointFields[i]), name, "Node",baseMeshObject)

         for i in range(len(CellFields)):
           name = 'CField'+str(i)
           if len(CellFields) == len(CellFieldsNames):
               name = CellFieldsNames[i]

           self.cellfieldsStorage[name] = self.__WriteAttribute(np.array(CellFields[i]), name, "Cell",baseMeshObject)

         for i in range(len(GridFields)):

           name = 'GField'+str(i)
           if len(GridFields) == len(GridFieldsNames):
               name = GridFieldsNames[i]

           self.gridfieldsStorage[name] = self.__WriteAttribute(np.array(GridFields[i]), name, "Grid",baseMeshObject)

    def __WriteNodesTagsElementsTags(self,baseMeshObject,PointFieldsNames,CellFieldsNames):
         for tag in baseMeshObject.nodesTags:
             if tag.name in PointFieldsNames:
                 name = "Tag_" + tag.name
             else:
                 name = tag.name

             data = np.zeros((baseMeshObject.GetNumberOfNodes(),1),dtype=np.int8)
             data[baseMeshObject.nodesTags[tag.name].GetIds()] = 1;
             self.__WriteAttribute(np.array(data), name, "Node",baseMeshObject)

         #Cell Tags
         baseMeshObject.PrepareForOutput();

         if baseMeshObject.IsConstantRectilinear():
             celtags = baseMeshObject.GetNamesOfElemTagsBulk()
             GetElementsInTag = baseMeshObject.GetElementsInTagBulk
         else:
             celtags = baseMeshObject.GetNamesOfElemTags()
             GetElementsInTag = baseMeshObject.GetElementsInTag

         for tagname in celtags:
             if tagname in CellFieldsNames:
                 name = "Tag_" + tagname
             else:
                 name = tagname
             data = GetElementsInTag(tagname)
             res = np.zeros((baseMeshObject.GetNumberOfElements(),1),dtype=np.int8)
             res[data] = 1;

             self.__WriteAttribute(np.array(res), name, "Cell", baseMeshObject)

    def __WriteTime(self):
        """ this function is called by the WriteTail, this function must NOT change
         the state of the instance, also no writting to binary neither hdf5 files """
        if self.isOpen() and self.__printTimeInsideEachGrid == False:
            #self.filePointer.write('<Time TimeType="List">\n')
            #self.__WriteDataItem(self.timeSteps)
            #self.filePointer.write('</Time>\n')
            self.filePointer.write('<Time Value="'+ (" ".join(str(x) for x in self.timeSteps)) +'"/>\n')

    def __WriteDataItem(self,_data, _shape= None,name=None):



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
                raise(Exception('Output Not implemented for data of type ' + str(type(data[0]))))                                                         # pragma: no cover

            dimension = ArrayToString(shape)


            if self.isBinary() and len(data) > self.__XmlSizeLimit:# pragma: no cover

                # to test this feature a big file must be created (so we dont test it)
                if self.__binarycpt > (2**30) :
                    self.__binaryFilePointer.close()
                    self.NewBinaryFilename()
                    self.__binaryFilePointer = open (self.__binFileName, "wb")
                    self.__binarycpt = 0

                gsdata,gsstorage = self.globalStorage.get(str(name),(None,None))

                datatowrite = data.ravel()
                if np.array_equal(gsdata,datatowrite):
                    binaryfile = gsstorage.filenameOnly
                    seek = gsstorage.offset
                    gsstorage.usedByNInstances += 1
                    self.globalStorage[str(name)] = (gsdata,gsstorage)
                    res = gsstorage
                else:
                    res = BinaryStorage(data=data,filePointer=self.__binaryFilePointer)
                    res.filenameOnly = self.__binFileNameOnly
                    binaryfile = self.__binFileNameOnly
                    seek = self.__binarycpt
                    data.ravel().tofile(self.__binaryFilePointer)
                    self.__binarycpt += s*len(datatowrite)
                    if len(self.globalStorage) > self.maxStorageSize-1:
                        usage = [x[1].usedByNInstances for x in self.globalStorage.values()]
                        #print(usage)
                        minUsage = min(usage)
                        newGlobalStorage = {}
                        GT = len(self.globalStorage)
                        for i,d in self.globalStorage.items():
                            if d[1].usedByNInstances == minUsage:
                                GT -= 1
                                if GT < self.maxStorageSize-1:
                                    break
                                continue
                            newGlobalStorage[i]=d
                        self.globalStorage = newGlobalStorage
                    self.globalStorage[str(name)] = (datatowrite,res)


                self.filePointer.write(' <DataItem Format="Binary"'+
                ' NumberType="'+typename+'"'+
                ' Dimensions="'+dimension+'" '+
                ' Seek="'+str(seek)+'" '+
                ' Endian="Native" '+
                ' Precision="'+str(s)+'" '+

                ' Compression="Raw" >')
                self.filePointer.write(binaryfile)
                self.filePointer.write('</DataItem>\n')

                return res
            else:
                self.filePointer.write(' <DataItem Format="XML" NumberType="'+typename+'" Dimensions="'+dimension+'">')
                self.filePointer.write(" ".join(str(x) for x in data.ravel()))

                self.filePointer.write('</DataItem>\n')
                return None


from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".xdmf",XdmfWriter)
RegisterWriterClass(".xmf",XdmfWriter)

def WriteTest(tempdir,Temporal, Binary):

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,3,4]);
    print(" a  ",myMesh.GetDimensions())
    print(" b  ",myMesh.structElements.GetNumberOfElements())
    print(myMesh.GetNumberOfElements())
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

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh

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

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    import BasicTools.Containers.UnstructuredMesh as UM
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOfTriangles

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


    CRM2D = ConstantRectilinearMesh(2);

    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetXmlSizeLimit(0)
    writer.SetBinary(True)
    writer.SetMultidomain()
    writer.Open(filename=tempdir+'testdirectTwoDomains.xdmf');


    writer.Write(CRM2D, PointFields = [ np.zeros((CRM2D.GetNumberOfNodes(),) ).astype(np.float32), np.zeros((CRM2D.GetNumberOfNodes(),) ).astype(np.int) ],
                                                            PointFieldsNames = ["Test", "testint"],
                                                            CellFields= [ np.arange(CRM2D.GetNumberOfElements()*2 ).astype(np.float64)],
                                                            CellFieldsNames = [ "TestV"] );
    CRM3D = ConstantRectilinearMesh(3)
    writer.Write(CRM3D)

    writer.Close()
    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(filename = tempdir+'testdirectTwoDomains.xdmf')


    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetXmlSizeLimit(0)
    writer.SetBinary(True)
    writer.SetMultidomain()
    writer.SetTemporal()
    writer.Open(filename=tempdir+'testdirectTwoDomainsTwoSteps.xdmf');
    writer.Write(CRM2D)
    writer.Write(CRM3D)
    writer.MakeStep(Time=1.5)
    writer.Write(CRM3D)
    writer.Write(CRM2D)
    writer.Close()

    if GUI:
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(filename = tempdir+'testdirectTwoDomainsTwoSteps.xdmf')

    ####### work for PXMF 2.0 #################
    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetXmlSizeLimit(0)
    writer.SetBinary(True)
    writer.SetParafac(True)
    writer.Open(filename=tempdir+'parafac.pxdmf');
    from BasicTools.Containers.UnstructuredMeshTools import  CreateMeshFromConstantRectilinearMesh as CMFCRM
    from BasicTools.Containers.UnstructuredMeshTools import  CreateUniformMeshOfBars

    mesh1DTime = CreateUniformMeshOfBars(2,5,10)
    mesh1DTime.props['ParafacDims'] = 1
    mesh1DTime.props['ParafacDim0'] = "T"

    cmesh2D = ConstantRectilinearMesh(2)
    cmesh2D.SetDimensions([4,4])
    mesh2DParametres = CMFCRM(cmesh2D)
    mesh2DParametres.props['ParafacDims'] = 2
    mesh2DParametres.props['ParafacDim0'] = "Px"
    mesh2DParametres.props['ParafacDim1'] = "Py"

    cmesh3D = ConstantRectilinearMesh(3)
    cmesh3D.SetDimensions([8,8,8])
    mesh3DSpace = CMFCRM(cmesh3D)

    from BasicTools.FE.IntegrationsRules import LagrangeIsoParam

    IntegrationPointData = np.arange(mesh2DParametres.GetNumberOfElements()*len(LagrangeIsoParam[EN.Quadrangle_4][1]) )+0.1

    IntegrationPointDatas = [IntegrationPointData]
    writer.Write(mesh2DParametres,
                 IntegrationRule=LagrangeIsoParam,
                 IntegrationPointData=IntegrationPointDatas,
                 IntegrationPointDataNames=["IPId_0"])

    print(IntegrationPointDatas)
    writer.iptorage["IPId_0"].UpdateHeavyStorage(IntegrationPointData+10)

    writer.Write(mesh1DTime, CellFields = [np.arange(mesh1DTime.GetNumberOfElements())+0.1 ], CellFieldsNames=["IPId_0"])
    writer.Write(mesh3DSpace, CellFields = [np.arange(mesh3DSpace.GetNumberOfElements())+0.1 ], CellFieldsNames=["IPId_0"])
    writer.Close()

    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(filename = tempdir+'parafac.pxdmf')

    from BasicTools.IO.XdmfReader import XdmfReader  as XR
    f = XR(filename = tempdir+'testdirect.xdmf' )
    f.lazy = False;
    f.Read();
    print(tempdir+'testdirect.xdmf')
    f.xdmf.GetDomain(0).GetGrid(0).GetFieldsOfType("Cell")
    print(f.xdmf.GetDomain(0).GetGrid(0).attributes )

    #print("Structured Mesh in 3D")
    #import BasicTools.FE.StructuredMesh as SM
    #SM3D = SM.StructuredMesh();
    #WriteMeshToXdmf(tempdir+'StructuredMesh.xdmf',SM3D, PointFields = [ np.arange(SM3D.GetNumberOfNodes()) ],
    #                                                                   PointFieldsNames = ["Test"],
    #                                                                   CellFields= [ np.arange(SM3D.GetNumberOfElements()*3)],
    #                                                                   CellFieldsNames = [ "TestV"] );
    #from BasicTools.IO.XdmfReader import XdmfReader  as XR
    #f = XR(filename = tempdir+'StructuredMesh.xdmf' )
    #f.lazy = False;
    #f.Read();

    domain = f.xdmf.GetDomain(0)
    grid  = domain.GetGrid(0)
    grid.GetFieldsOfType("Cell")

    ### Domain Decomposition and Parafac

    res = CheckIntegrityDDM(GUI=GUI)
    if res.lower() != "ok": return res

    return 'ok'

def CheckIntegrityDDM(GUI=False):
    """ this test function can be lauched using the mpirun -n 2 ...
    to test the writer in mpi mode
    """

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()


    from BasicTools.Containers.UnstructuredMeshTools import  CreateUniformMeshOfBars

    mesh1D = CreateUniformMeshOfBars(0,5,10)
    mesh1D.props['ParafacDims'] = 1
    mesh1D.props['ParafacDim0'] = "T"

    from BasicTools.Helpers.MPIInterface import MPIInterface as MPI

    writer = XdmfWriter()
    writer.SetFileName(None)
    writer.SetXmlSizeLimit(0)
    writer.SetBinary(True)
    writer.SetMultidomain()
    writer.SetParafac(True)
    mpi = MPI()
    print("rank ", mpi.rank)
    writer.Open(filename=tempdir+'DDM_parafac.pxdmf');


    mpi = MPI()
    if mpi.Rank() == 0:
        mesh1D.props['ParafacDim0'] = "D0_P0"
        writer.Write(mesh1D, CellFields = [np.arange(mesh1D.GetNumberOfElements())+0.1 ], CellFieldsNames=["IPId_0"])

        mesh1D.nodes[:,0] += 1
        mesh1D.props['ParafacDim0'] = "D0_P1"
        writer.Write(mesh1D, CellFields = [np.arange(mesh1D.GetNumberOfElements())+0.1 ], CellFieldsNames=["IPId_0"])

    if mpi.IsParallel() == 0 :
        writer.NextDomain()

    if mpi.Rank() == mpi.size-1 :

        mesh1D.nodes[:,0] += 1
        mesh1D.props['ParafacDim0'] = "D1_P0"
        writer.Write(mesh1D, CellFields = [np.arange(mesh1D.GetNumberOfElements())+0.1 ], CellFieldsNames=["IPId_0"])

        mesh1D.nodes[:,0] += 1
        mesh1D.props['ParafacDim0'] = "D1_P1"
        writer.Write(mesh1D, CellFields = [np.arange(mesh1D.GetNumberOfElements())+0.1 ], CellFieldsNames=["IPId_0"])

    writer.Close()
    return "ok"


if __name__ == '__main__':
    #from BasicTools.Helpers.Tests import TestTempDir
    #TestTempDir.SetTempPath("/tmp/BasicTools_Test_Directory__is42mzi_safe_to_delete/")
    print(CheckIntegrity(GUI=False)) # pragma: no cover

