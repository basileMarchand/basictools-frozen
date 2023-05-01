# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Xdmf file reader
"""

import os
import sys
import numpy as np
import xml.sax

from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.IO.XdmfTools import FieldNotFound
from BasicTools.NumpyDefs import PBasicFloatType

import BasicTools.Containers.ElementNames as ElementNames


def ReadXdmf(fileName):
    """Function API for reading an xdmf file

    Parameters
    ----------
    fileName : str
        name of the file to be read

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    obj = XdmfReader(filename=fileName)
    return obj.Read()



class Xdmfbase(object):
    """ Base class for all the xdmf Releated objects """

    def ReadAttribute(self,attrs,name,default=None):
        """" Helper to read attributes"""
        if name in attrs:
            return attrs.get(name)
        else:
            if default is None :
                raise KeyError("Key :'" +name+"' not available (you can add a default value to bypass this error)")
            else:
                return default

    def TreatCDATA(self):
        """Default inplementation to read the heavy data"""
        pass
    def __AttributsToStr(self):
        """ Helper function to make easier the introspection """
        import inspect
        res = ''
        TFormat.II()
        attributes = inspect.getmembers(self, lambda a:not inspect.isroutine(a) )
        attributes = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__')  )and a[0][0].isupper() ]
        #print(attributes)
        for a in attributes:
            #print(a)
            res += TFormat.GetIndent() + str(a[0]) +' : '+ str(a[1]) +    '\n'
        TFormat.DI()
        return res

class Xdmf(Xdmfbase):
    """ Top class for the xdmf Document """
    def __init__(self):
        self.domains = []

    def GetDomain(self,nameornumber):
        """ Get a Grid (XdmfGrid) using a name or a number """
        if isinstance(nameornumber,str):
            for g in self.domains:
                if g.Name == nameornumber :
                    return g
            raise Exception("Domain '"+nameornumber +"' not Found")
        else:
            return self.domains[nameornumber]

    def __str__(self):
        res = 'Xdmf\n'
        TFormat.II()
        for d in self.domains:
            res += d.__str__()
        TFormat.DI()
        return res

class XdmfDomain(Xdmfbase):
    """A Domain. Can contain many grids"""
    def __init__(self):
        self.grids=[]
        self.Name = ''
        self.informations =[]

    def GetGrid(self, nameornumber):
        """ Get a Grid (XdmfGrid) using a name or a number """
        if isinstance(nameornumber,str):
            for g in self.grids:
                if g.Name == nameornumber : return g
            raise
        else:
            return self.grids[nameornumber]

    def ReadAttributes(self,attrs):

        self.Name = self.ReadAttribute(attrs,'Name','')

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfDomain\n'
        TFormat.II()
        for g in self.grids:
            res += g.__str__()
        TFormat.DI()
        return res



class XdmfGrid(Xdmfbase):
    """ a Grid: contains a mesh (poinst and connectivity ) and fields by (element, nodes, grid)"""
    def __init__(self):
        self.informations= []
        self.topology = None
        self.geometry = None
        self.attributes = []
        self.Name = ''
        self.GridType = 'Uniform'
        self.time = None

    def GetSupport(self):
        """Returns the support defined in the Xdmf file

        Returns
        -------
        UnstructuredMesh or ConstantRectilinearMesh
            output mesh object containing reading result
        """
        if self.geometry.Type == "ORIGIN_DXDYDZ":
            from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
            res = ConstantRectilinearMesh()
            res.SetOrigin(self.geometry.GetOrigin())
            res.SetSpacing(self.geometry.GetSpacing())
            swaper = [2,1,0]
            res.SetDimensions(self.topology.GetDimensions()[swaper])

        if self.geometry.Type == "XYZ":
            from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
            res = UnstructuredMesh()
            res.nodes = self.geometry.GetNodes().reshape((-1,3))
            res.elements = self.topology.GetConnectivity()
            res.GenerateManufacturedOriginalIDs()
            res.PrepareForOutput()
            if np.linalg.norm(res.nodes[:,2]) == 0:
                from BasicTools.Containers.UnstructuredMeshModificationTools import LowerNodesDimension
                res = LowerNodesDimension(res)

        if self.geometry.Type == "XY":
            from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
            res = UnstructuredMesh()
            res.nodes = self.geometry.GetNodes().reshape((-1,2))
            res.elements = self.topology.GetConnectivity()
            res.GenerateManufacturedOriginalIDs()
            res.PrepareForOutput()

        nodeFields = {name : self.GetPointField(name)[0] for name in  self.GetPointFieldsNames()}
        elemFields = {name : self.GetCellField(name)[0] for name in  self.GetCellFieldsNames()}

        for name,data in nodeFields.items():
            if data.shape[0] == res.GetNumberOfNodes() and data.dtype == np.int8 and np.min(data) == 0 and np.max(data)==1:
                res.nodesTags.CreateTag(name).SetIds(np.where(data==1)[0])
            else:
                res.nodeFields[name] = data

        for name, data in elemFields.items():
            if data.shape[0] == res.GetNumberOfElements() and data.dtype == np.int8 and np.min(data) == 0 and np.max(data)==1:
                res.AddElementsToTag(np.where(data)[0],tagname=name)
            else:
                res.elemFields[name] = data

        return res

    def ReadAttributes(self,attrs):
        self.Name = self.ReadAttribute(attrs,'Name',default="")
        self.GridType = self.ReadAttribute(attrs,'CollectionType',default="Uniform")

    def HasField(self,name):
        for a in self.attributes:
            if a.Name == name : return True
        return False

    def GetTime(self):
        """Returns time at which data is read

        Returns
        -------
        float
            time at which data is read
        """
        if self.time is None:
            return None
        else:
            return self.time.Value[0]

    def GetFieldsNames(self):
        return [a.Name for a in self.attributes ]

    def GetPointFieldsNames(self):
        return [a.Name for a in self.attributes if a.Center == "Node"]

    def GetCellFieldsNames(self):
        return [a.Name for a in self.attributes if a.Center == "Cell"]

    def GetGridFieldsNames(self):
        return [a.Name for a in self.attributes if a.Center == "Grid" ]

    def GetPointFields(self):
        return self.GetFieldsOfType("Node")

    def GetPointField(self,name):
        return self.GetFieldsOfType("Node",name)

    def GetCellFields(self):
        return self.GetFieldsOfType("Cell")

    def GetCellField(self,name):
        return self.GetFieldsOfType("Cell",name)

    def GetGridFields(self):
        return self.GetFieldsOfType("Grid")

    def GetGridField(self,name):
        return self.GetFieldsOfType("Grid",name)

    def GetFieldsOfType(self,ftype, name = None):
        res = []
        for a in self.attributes:
            if name is not None:
                if a.Name != name:
                    continue
            if a.Center == ftype :
                data = a.dataitems[0].GetData()
                if self.geometry.Type == "ORIGIN_DXDYDZ":
                    if a.Type == "Vector":
                        data = data.transpose(2,1,0,3)
                    else:
                        data.shape = self.topology.GetDimensions()
                        data = data.transpose(2,1,0).ravel()
                res.append(np.copy(data))
        return res

    def GetFieldData(self,name,noTranspose=False):
        for att in self.attributes:
            if att.Name == name :
                data = att.dataitems[0].GetData()
                if self.geometry.Type == "ORIGIN_DXDYDZ":
                    if att.Type == "Vector":
                        if noTranspose:
                            return data
                        else:
                            return data.transpose(2,1,0,3)
                    else:
                        if noTranspose:
                            return data
                        else:
                            return data.transpose(2,1,0)
                return data
        raise FieldNotFound(name)

    def GetFieldTermsAsColMatrix(self,fieldname, offset = 0):
        #first we check the padding max 8 zeros of padding
        padding = 0
        for i in range(8):
            padding = i
            ss = fieldname+str(offset).zfill(padding)
            #print ss
            if self.HasField(ss):
                break
            #print(padding)
        else:
            raise FieldNotFound(fieldname)
        #first we check the number of terms
        #print('padding ' + str(padding))
        cpt =0
        while(self.HasField(fieldname+ str(offset+cpt).zfill(padding) )):
            cpt +=1
        # get the firt data to get the type and the size
        d_0 = self.GetFieldData(fieldname+str(offset).zfill(padding) )

        # now we allocate the np matrix for the data
        res = np.empty([d_0.size, cpt], dtype=d_0.dtype)
        res[:,0] = d_0.reshape(d_0.size)

        for i in range(1,cpt):
            res[:,i] = self.GetFieldData(fieldname+str(offset+i).zfill(padding) ).reshape(d_0.size)

        return res

    def GetFieldTermsAsTensor(self,fieldname,sep='_',offseti=0,offsetj=0):
        from itertools import product
        #first we check the padding max 8 zeros of padding
        paddingi = 0
        paddingj = 0
        for i,j in product(range(8),range(8)):
            paddingi = i
            paddingj = j
            ss = fieldname + str(offseti).zfill(paddingi) + sep + str(offsetj).zfill(paddingj)
            if self.HasField(ss):
                break
        else:
            raise FieldNotFound(fieldname+"*"+ str(offseti)+ sep + "*" +str(offsetj))


        #first we check the number of terms
        cpti =0
        while(self.HasField(fieldname+str(offseti+cpti).zfill(paddingi) + sep + str(offsetj).zfill(paddingj) )):
            cpti +=1

        cptj =0
        while(self.HasField(fieldname+str(offseti).zfill(paddingi)  + sep + str(offsetj+cptj).zfill(paddingj) )):
            cptj +=1

        # get the firt data to get the type and the size
        d_0_0 = self.GetFieldData(fieldname+str(offseti).zfill(paddingi)  + sep + str(offsetj).zfill(paddingj) )

        # now we allocate the np matrix for the data
        res = np.empty([d_0_0.size, cpti, cptj], dtype=d_0_0.dtype)
        for i in range(0,cpti):
            for j in range(0,cptj):
                #print(fieldname + str(offseti+i).zfill(paddingi)  + sep + str(offsetj+j).zfill(paddingj))
                res[:,i,j] = self.GetFieldData(fieldname + str(offseti+i).zfill(paddingi)  + sep + str(offsetj+j).zfill(paddingj) ).reshape(d_0_0.size)
        return res

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfGrid\n'
        TFormat.II()
        for info in self.informations:
            res += info.__str__()
        res += TFormat.GetIndent() +'Name : '+ self.Name +  '\n'
        res += self.geometry.__str__()
        res += self.topology.__str__()
        for att in self.attributes:
            res += att.__str__()
        TFormat.DI()
        return res

class XdmfInformation(Xdmfbase):
    """ class for extra information in the xmdf file"""
    def __init__(self):
        self.Name = ''
        self.Value = ''

    def ReadAttributes(self,attrs):
        self.Name = self.ReadAttribute(attrs,'Name')
        self.Value = self.ReadAttribute(attrs,'Value')

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfInformation\n'
        res += self._Xdmfbase__AttributsToStr()
        return res

class XdmfTopology(Xdmfbase):
    """ XdmfTopology class: stores the connectivity of the Grid"""
    def __init__(self):
        self.dataitems= []
        self.Dimensions = None
        self.Type = None

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfTopology\n'
        res += self._Xdmfbase__AttributsToStr()
        return res

    def ReadAttributes(self,attrs):

        self.Type = self.ReadAttribute(attrs,'Type', default=-1)
        if self.Type == -1:
            self.Type = self.ReadAttribute(attrs,'TopologyType')

#        if self.Type != "Mixed":
        try:
            self.Dimensions = np.array(self.ReadAttribute(attrs,'Dimensions').split(), dtype='int')[::-1]
        except :
            pass
#             d = self.dataitems[0].ReadAttribute(attrs,'Dimensions').split()
#             print(d)
#             self.Dimensions = np.array(d, dtype='int')[::-1]
#             print(self.Dimensions)


    def GetConnectivity(self):
        from BasicTools.IO.XdmfTools import XdmfNumberToEN
        from BasicTools.IO.XdmfTools import XdmfNameToEN
        from BasicTools.Containers.UnstructuredMesh import ElementsContainer,AllElements
        res = AllElements()
        from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh


        if self.Type == "Mixed":
            rawdata = self.dataitems[0].GetData()
            ts = len(rawdata)
            cpt = 0
            elcpt = 0
            while cpt < ts:
                #print(rawdata)
                #print(XdmfNumberToEN)
                EN = XdmfNumberToEN[rawdata[cpt]]
                if rawdata[cpt] == 0x2 or rawdata[cpt] == 0x1 :
                    cpt += 1
                cpt += 1
                elements = res.GetElementsOfType(EN)
                connectivity =rawdata[cpt:elements.GetNumberOfNodesPerElement()+cpt]
                elements.AddNewElement(connectivity,elcpt)
                elcpt +=1
                cpt += elements.GetNumberOfNodesPerElement()
            return res
        elif self.Type == "3DCoRectMesh" :
            raise
        elif self.Type == "3DSMesh" :
            m = ConstantRectilinearMesh()
            m.SetDimensions(self.Dimensions)

            #m.GenerateFullConnectivity()
            return m.elements
        else :
            en = XdmfNameToEN[self.Type]
            data = ElementsContainer(en)
            data.connectivity = self.dataitems[0].GetData().ravel()
            size = (len(data.connectivity)//ElementNames.numberOfNodes[en],ElementNames.numberOfNodes[en])
            data.connectivity.shape = size
            data.cpt = size[0]
            res[en] = data
            return res



    def GetDimensions(self):
        if self.Dimensions is None:
            return self.dataitems[0].Dimensions

        return self.Dimensions[::-1]

class XdmfGeometry(Xdmfbase):
    """XdmfGeometry class:  stores the point positions """

    def __init__(self):
        self.dataitems= []
        self.Type = None

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfGeometry\n'
        res += self._Xdmfbase__AttributsToStr()
        return res

    def ReadAttributes(self,attrs):
        if 'Type' in attrs:
            self.Type = self.ReadAttribute(attrs,'Type')
        else:
            self.Type = self.ReadAttribute(attrs,'GeometryType')

    def GetOrigin(self):
        if self.Type == "ORIGIN_DXDYDZ":
            #self.Read()
            return self.dataitems[0].GetData()[::-1]
        else:
            raise Exception# pragma: no cover

    def GetSpacing(self):
        if self.Type == "ORIGIN_DXDYDZ":
            #self.Read()
            return self.dataitems[1].GetData()[::-1]
        else:
            raise Exception# pragma: no cover

    def GetNodes(self):
        if self.Type == "XYZ" or self.Type == "XY":
            return self.dataitems[0].GetData()
        else:
            raise Exception# pragma: no cover



class XdmfAttribute(Xdmfbase):
    """  XdmfAttribute class: to store the data over the grids """

    def __init__(self):
        self.dataitems= []
        self.Name = ''
        self.Type =  ''
        self.Center = ''
        self.CDATA = ''

    def ReadAttributes(self,attrs):
        self.Name = self.ReadAttribute(attrs,'Name')
        try :
            self.Type = self.ReadAttribute(attrs,'Type')
        except:
            self.Type = self.ReadAttribute(attrs,'AttributeType')
        self.Center = self.ReadAttribute(attrs,'Center')

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfAttribute\n'
        TFormat.II()
        res +=TFormat.GetIndent() + 'Name : '+ self.Name +  '\n'
        res +=TFormat.GetIndent() + 'Type : '+ self.Type +  '\n'
        res +=TFormat.GetIndent() + 'Center : '+ self.Center +  '\n'
        for d in self.dataitems:
            res += d.__str__()
        TFormat.DI()
        return res

class XdmfTime(Xdmfbase):
    """  XdmfTime class: to store the data over the grids """

    def __init__(self):
        self.Value = None

    def ReadAttributes(self,attrs):
        self.Value = np.array(self.ReadAttribute(attrs,'Value').split(), dtype=PBasicFloatType)

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfTime'
        TFormat.II()
        res +=TFormat.GetIndent() + 'Value : '+ str(self.Value) +  '\n'
        TFormat.DI()
        return res

class XdmfDataItem(Xdmfbase):
    """ XdmfDataItem class : class to manage the reading of the data Heavy and light """

    def __init__(self):
        self.Type = None
        self.Dimensions= []
        self.Precision = None
        self.Format = None
        self.Data = []
        self.CDATA = ''
        self.Seek = 0
        self.Endian = "Native"

    def ReadAttributes(self,attrs, path):
        self.path = path
        self.Dimensions = np.array(self.ReadAttribute(attrs,'Dimensions').split(), dtype='int')
        if 'DataType' in attrs:
            self.Type = self.ReadAttribute(attrs,'DataType')
        elif 'NumberType' in attrs:
            self.Type = self.ReadAttribute(attrs,'NumberType')
        else:
            self.Type = 'float'

        if((self.Type.lower() == 'float' or self.Type.lower() == 'int' )and 'Precision' in attrs):
            self.Precision = int(self.ReadAttribute(attrs,'Precision'))
        self.Format = self.ReadAttribute(attrs,'Format')
        self.Seek = int(self.ReadAttribute(attrs,'Seek',0))
        self.Endian = self.ReadAttribute(attrs,'Endian',"Native")

    def __str__(self):
        res = TFormat.GetIndent() +'XdmfDataItem \n'
        TFormat.II()
        res += TFormat.GetIndent() + 'Type : '+ self.Type +  '\n'
        res += TFormat.GetIndent() + 'Dimensions : '+ str(self.Dimensions) +  '\n'
        if self.Type.lower() == 'float' :
            res += TFormat.GetIndent() + 'Precision : '+ str(self.Precision) +  '\n'
        elif  self.Type.lower() == 'int' :
            res += TFormat.GetIndent() + 'Precision : '+ str(self.Precision) +  '\n'

        res += TFormat.GetIndent() + 'Format : '+ str(self.Format) +  '\n'
        res += TFormat.GetIndent() + 'Data : \n  '+  TFormat.GetIndent()  +str(self.Data) +  '\n'
        res += TFormat.GetIndent() + 'CDATA : \n  '+  TFormat.GetIndent()  +str(self.CDATA) +  '\n'
        TFormat.DI()
        return res

    def TreatCDATA(self):
        if len(self.CDATA) == 0:
            return


        if self.Format  == 'XML':
            if self.Type.lower() =='float':
                if self.Precision == 4:
                    numpytype = 'float32'
                else:
                    numpytype = 'float_'
            elif self.Type.lower() == 'int':
                numpytype = 'int_'
            elif self.Type.lower() == 'char':
                numpytype = 'int8'

            self.Data= np.array(self.CDATA.split(), dtype=numpytype)
            self.Data = self.Data.reshape(self.Dimensions)
            self.CDATA = ''
        elif self.Format  == 'HDF':

            filename,dataSetPath  = str(self.CDATA).lstrip().rstrip().split(":")
            #print(dataSetPath)
            from h5py import File as __File
            f = __File(os.path.join(self.path, filename),'r')
            #print(f[dataSetPath])
            self.Data =  np.array(f[dataSetPath])
            self.Data = self.Data.reshape(self.Dimensions)
            #print(self.Data)
            self.CDATA = ''

        elif self.Format  == 'Binary':
            if self.Type.lower() == 'float':
                if self.Precision == 4 :
                    numpytype = 'float32'
                else:
                    numpytype = 'float_'
            elif self.Type.lower() =='int':
                if self.Precision == 4:
                    numpytype = np.int32
                elif self.Precision == 8 :
                    numpytype = np.int64
                else:
                    raise Exception(f"Dont know how to treat this type of Precision: {self.Precision}" )
            elif self.Type.lower() == 'char':
                numpytype = 'int8'

            binfilename  = str(self.CDATA).lstrip().rstrip()
            binfile = open (os.path.join(self.path, binfilename ), "rb")
            binfile.seek(self.Seek)

            self.Data = np.fromfile(binfile, dtype=numpytype, count=np.prod(self.Dimensions))
            if self.Endian == "Native" :
                pass
            elif self.Endian == "Big" and sys.byteorder == "little":
                self.Data.byteswap(inplace=True)
            elif self.Endian == "Little" and sys.byteorder == "big":
                self.Data.byteswap(inplace=True)

            self.Data.shape = self.Dimensions
            binfile.close()
            self.CDATA = ''
        else :
            raise Exception("Heavy data in format '" + self.Format + "' not suported yet") # pragma: no cover

    def GetData(self):
        self.TreatCDATA()
        return self.Data


class XdmfReader(xml.sax.ContentHandler):
    """Xdmf Reader class
    """
    def __init__(self,filename=''):
        super(XdmfReader,self).__init__()
        self.xdmf = Xdmf()
        self.pile = []
        self.path = ''
        self.filename = ''
        self.SetFileName(filename)
        self.readed = False
        self.lazy = True
        self.canHandleTemporal = True
        self.time = np.array([])
        self.timeToRead = -1
        self.encoding = None

    def Reset(self):
        """Resets the reader
        """
        self.xdmf = Xdmf()
        self.pile = []
        self.readed = False
        self.timeToRead = -1
        self.time = np.array([])

    def SetFileName(self, filename):
        """Sets the name of file to read

        Parameters
        ----------
        filename : str
            file name to set
        """
        #if same filename no need to read the file again
        if  self.filename == filename: return

        self.readed = False
        self.filename = filename
        self.path = os.path.dirname(filename)

    def ReadMetaData(self):
        """Function that performs the reading of the metadata of a xdmf file:
        sets the time attribute of the reader
        """
        self.lazy = True
        times = []
        self.Read()
        for i, grid in enumerate(self.xdmf.domains[0].grids):
            t  = grid.GetTime()
            if t == None:
                continue
            times.append(t)
        self.time = np.array(times)

    def GetAvailableTimes(self):
        """Returns the available times at which data can be read

        Returns
        -------
        np.ndarray
            available times at which data can be read
        """
        return self.time

    def SetTimeToRead(self, time=None, timeIndex=None):
        """Sets the time at which the data is read

        Parameters
        ----------
        time : float, optional
            time at which the data is read, by default None
        timeIndex : int, optional
            time index at which the data is read, by default None

        Returns
        -------
        int
            time index at which the data is read
        """
        if time is not None:
            self.timeToRead = time
        else:
            if timeIndex is not None:
                self.timeToRead = self.time[timeIndex]
            else:
                self.timeToRead = -1

    def Read(self, fileName = None):
        """Function that performs the reading of the data defined in an ut file

        Parameters
        ----------
        fileName : str, optional
            name of the file to read, by default None

        Returns
        -------
        UnstructuredMesh or ConstantRectilinearMesh
            output mesh object containing reading result
        """

        if self.timeToRead == -1.:
            timeIndex = len(self.time)-1
        else:
            if len(self.time) == 0:
                timeIndex = -1
            else:
                timeIndex = np.argmin(abs(self.time - self.timeToRead ))

        if fileName is not None:
            self.SetFileName(fileName)
            self.readed = False

        if len(self.filename) == 0 :
            raise Exception('Need a filename ')

        # read only one time the filel
        if self.readed: return
        self.Reset()

        thefile = open(self.filename,"r")
        # to deactivate the DTD validation
        parser = xml.sax.make_parser()
        parser.setContentHandler(self)
        parser.setFeature(xml.sax.handler.feature_external_ges, False)
        parser.parse(thefile)
        thefile.close()
        return self.xdmf.GetDomain(-1).GetGrid(timeIndex).GetSupport()

    # this a a overloaded function (must start with lower case)      !!!!!
    def startElement(self, name, attrs):

        if name == "Xdmf":
            self.pile.append(self.xdmf)
            return

        father = self.pile[-1]

        if name == "Domain":
            res = XdmfDomain()
            res.ReadAttributes(attrs)
            father.domains.append(res)
        elif name == 'Grid':
            res = XdmfGrid()
            res.ReadAttributes(attrs)
            # for the moment we use a flat representation (no GridType collection)
            if "GridType" in attrs and attrs.get("GridType").lower() == "collection" :
                res = father
            else :
                father.grids.append(res)
        elif name == 'Information':
            res = XdmfInformation()
            res.ReadAttributes(attrs)
            father.informations.append(res)
        elif name == 'Topology':
            res = XdmfTopology()
            res.ReadAttributes(attrs)
            father.topology = res
        elif name == 'Geometry':
            res = XdmfGeometry()
            res.ReadAttributes(attrs)
            father.geometry = res
        elif name == 'DataItem':
            res = XdmfDataItem()
            res.ReadAttributes(attrs, self.path)
            father.dataitems.append(res)
        elif name == 'Attribute':
            res = XdmfAttribute()
            res.ReadAttributes(attrs)
            father.attributes.append(res)
        elif name == "Time":
            res = XdmfTime()
            res.ReadAttributes(attrs)
            father.time = res

        else:
            raise Exception("Unkown tag :  '"+ name +"' Sorry!") # pragma: no cover

        self.pile.append(res)

    # this a a overloaded function (must start with lower case)      !!!!!
    def characters(self, content):
        #print("*** " + content + " +++ ")
        father = self.pile[-1]
        if isinstance(father,XdmfDataItem):
            father.CDATA += content  # Note: '+=', not '='

    # this a a overloaded function (must start with lower case)      !!!!!
    def endElement(self, name):
        if self.lazy :
            self.pile.pop()
        else:
            self.pile.pop().TreatCDATA()

    def __str__(self):
        res = ''
        for d in self.xdmf.domains:
            res = d.__str__() + "\n"
        return res

def GetTensorRepOfField(domain,fieldname):
    """ Get The tensor representation of a field in a xdmf domain

    Get The tensor representation of a field in a xdmf domain, for the moment
    works for Cannonic and Train Tensor formats
    """
    import BasicTools.TensorTools.Formats as st
    # we check the nature of the field (info in the parent domain)
    fieldtype = 'CP'
    for info in domain.informations:
        if info.Name == (fieldname+"_Type"):
            fieldtype = info.Value

    if fieldtype == "CP":
        res = st.CanonicTensor()
        for g in  domain.grids:
            G = st.FullTensor(nRanks=2)
            G.SetRanksNames(["-".join([i.Value for i in g.informations if not i.Name == 'Dims' ]) ,'alpha1'])
            G.array = g.GetFieldTermsAsColMatrix(fieldname+"_")
            res.AddSubTensor(G)
    elif fieldtype == "TT":
        res = st.TensorTrain()
        l = len(domain.grids)
        for cpt,g in  zip(range(0,l),domain.grids):
            if cpt == 0:
                G = st.FullTensor(nRanks=2)
                G.SetRanksNames(["-".join([i.Value for i in g.informations if not i.Name == 'Dims']) ,'alpha1'])
                G.array = g.GetFieldTermsAsColMatrix(fieldname+"_")
            elif cpt == l-1:
                G = st.FullTensor(nRanks=2)
                G.SetRanksNames(["-".join([i.Value for i in g.informations if not i.Name == 'Dims']) ,'alpha'+str(cpt)])
                G.array = g.GetFieldTermsAsColMatrix(fieldname+"_")
            else:
                G = st.FullTensor(nRanks=3)
                G.SetRanksNames(["-".join([i.Value for i in g.informations if not i.Name == 'Dims']) ,'alpha'+str(cpt+1),'alpha'+str(cpt)])
                G.array = g.GetFieldTermsAsTensor(fieldname+"_")
            res.AddSubTensor(G)
    return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".xdmf",XdmfReader)
RegisterReaderClass(".pxdmf",XdmfReader)

def CheckIntegrity():


    # All the try pasrt are to execute the code for error handling.
    # The code in the try(s) must raise exceptions
    # TODO Verified if the exceptions are the good ones

    # Get TestData directory
    from  BasicTools.TestData import GetTestDataPath
    TestDataPath = GetTestDataPath()

    # Create a Reader
    res = XdmfReader()
    try :
        res.Read()
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    res = XdmfReader(filename = TestDataPath + "Unstructured.xmf" )
    res.lazy = False
    res.Read()

    res = XdmfReader(filename = TestDataPath + "Unstructured.xmf" )
    # read only the xml part
    res.Read()


    rep = res.__str__()
    #print(rep)
    #Test xdmf part***********************
    a = res.xdmf.__str__()

    try:
        res.xdmf.GetDomain("Imaginary Domain")
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    res.xdmf.GetDomain("Dom 1")

    domain = res.xdmf.GetDomain(0)
    #Test domain part***********************
    try :
        domain.GetGrid("imaginary grid")
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    domain.GetGrid(0)
    grid  = domain.GetGrid("Grid")


    #Test Grid part**************************************************************


    names = grid.GetFieldsNames()
    if names[0] != 'RTData': raise Exception()

    grid.HasField(names[0])
    grid.HasField('toto')

    try:
        data = res.xdmf.GetDomain("Dom 1").GetGrid("Grid").GetFieldData('ImaginaryField')
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    grid.GetFieldTermsAsColMatrix('term_')

    try:
        grid.GetFieldTermsAsColMatrix('ImaginaryField_')
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass


    grid.GetFieldTermsAsTensor('TensorField_',offsetj=1)

    try:
        grid.GetFieldTermsAsTensor('ImaginaryField_',offsetj=1)
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    grid.GetFieldData('IntField')

    try:
        grid.GetFieldData('UnknownField')
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass


    data = res.xdmf.domains[0].GetGrid("Grid").GetFieldData('RTData')
    if  data[49] != 260.0: raise

    geo = res.xdmf.domains[0].GetGrid("Grid").geometry.dataitems[0].GetData()
    #print(geo)
    if geo[0,2] != -1: raise

    topo = res.xdmf.domains[0].GetGrid("Grid").topology.dataitems[0].GetData()
    #print(topo)
    if topo[2] != 1: raise

    ######################### Structured #########################
    res = XdmfReader(filename = TestDataPath + "Structured.xmf" )
    # read only the xml part
    res.Read()
    domain = res.xdmf.GetDomain(0)
    #Test domain part***********************

    grid  = domain.GetGrid(0)
    #print(grid.GetFieldsOfType("Node"))
    #print("--")
    #print(grid.GetFieldsOfType("Cell"))
    grid.geometry.GetOrigin()
    grid.geometry.GetSpacing()
    ##################################
    Example1()
    Example2()

    return 'OK'



def Example1():
    import BasicTools.TestData as test
    # Create a Reader
    reader = XdmfReader(filename = test.GetTestDataPath() + "Unstructured.xmf")
    # Do the reading (only the xml part, to read all the data set lazy to False)
    #res.lazy = False
    reader.Read()

    # Get the domaine "Dom 1"
    dom = reader.xdmf.GetDomain("Dom 1")

    # Get the first Grid
    grid = dom.GetGrid(0)
    grid.topology.GetDimensions()
    names = grid.GetPointFieldsNames()
    allFields = grid.GetPointFields()
    if len(allFields):
        print(grid.GetPointField(names[0]))

    names = grid.GetCellFieldsNames()
    allFields = grid.GetCellFields()
    if len(allFields):
        print(grid.GetCellField(names[0]))

    names = grid.GetGridFieldsNames()
    allFields = grid.GetGridFields()
    if len(allFields):
        print(grid.GetGridField(names[0]))

    #Get one field (or term)
    dataField1= grid.GetFieldData('RTData')

    dataField1= grid.GetFieldsOfType('Node')
    #print(dataField1.shape)
    #Get all the term as a matrix
    dataField2=  grid.GetFieldTermsAsColMatrix('term_')
    #print(dataField2.shape)
    #Get all the term as a matrix
    dataField3=  grid.GetFieldTermsAsTensor('TensorField_',offsetj=1)
    #print(dataField3.shape)


def Example2():
    import BasicTools.TestData as test

    reader = XdmfReader(filename = test.GetTestDataPath() + "TensorTestData.xmf")

    # Do the reading (only the xml part, to read all the data set lazy to False)
    #res.lazy = False
    reader.Read()


    # Get the domaine "Dom 1"
    dom = reader.xdmf.GetDomain(0)

    CT = GetTensorRepOfField(dom,'CanonicField')

    TT = GetTensorRepOfField(dom,'TTField')
    TT.CheckIntegrity(verbose=False)
    #print(TT)

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
