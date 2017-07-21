# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 09:47:58 2016

@author: d584808
"""


from BasicTools.Helpers.TextFormatHelper import TFormat
from BasicTools.IO.XdmfTools import FieldNotFound
import numpy as np


class Xdmfbase(object):
    """ Base class for all the xdmf Releated objects """

    def ReadAttribute(self,attrs,name,default=None):
        """" Helper to read attributes"""
        if name in attrs:
            return attrs.get(name)
        else:
            if( default==None) :
                raise KeyError("Key :'" +name+"' not avilable (you can add a default value to bypass this error)")
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
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
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
                if g.Name == nameornumber : return g;
            raise Exception("Domain '"+nameornumber +"' not Found");
        else:
            return self.domains[nameornumber]

    def __str__(self):
        res = 'Xdmf\n'
        TFormat.II()
        for d in self.domains:
            res += d.__str__();
        TFormat.DI()
        return res

class XdmfDomain(Xdmfbase):
    """A Domain. Can contain many grids"""
    def __init__(self):
        self.grids=[];
        self.Name = '';
        self.informations =[];

    def GetGrid(self,nameornumber):
        """ Get a Grid (XdmfGrid) using a name or a number """
        if isinstance(nameornumber,str):
            for g in self.grids:
                if g.Name == nameornumber : return g;
            raise
        else:
            return self.grids[nameornumber]

    def ReadAttributes(self,attrs):

        self.Name = self.ReadAttribute(attrs,'Name','')

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfDomain\n'
        TFormat.II()
        for g in self.grids:
            res += g.__str__();
        TFormat.DI()
        return res



class XdmfGrid(Xdmfbase):
    """ a Grid: contains a mesh (poinst and connectivity ) and fields by (element, nodes, grid)"""
    def __init__(self):
        self.informations= []
        self.topology = None;
        self.geometry = None;
        self.attributes = []
        self.Name ='';

    def GetSupport(self):
        if self.geometry.Type == "ORIGIN_DXDYDZ":
            from BasicTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh
            res = ConstantRectilinearMesh()
            res.SetOrigin(self.geometry.GetOrigin())
            res.SetSpacing(self.geometry.GetSpacing())
            res.SetDimensions(self.topology.GetDimensions())
            return res


    def ReadAttributes(self,attrs):
        self.Name = self.ReadAttribute(attrs,'Name')
        self.GridType = self.ReadAttribute(attrs,'CollectionType',default="Uniform")

    def HasField(self,name):
        for a in self.attributes:
            if a.Name == name : return True;
        return False

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

    def GetCellFields(self):
        return self.GetFieldsOfType("Cell")

    def GetGridFields(self):
        return self.GetFieldsOfType("Grid")

    def GetFieldsOfType(self,ftype):
        import numpy as np
        res = []
        for a in self.attributes:
            if a.Center == ftype :
                data = a.dataitems[0].GetData()
                if self.geometry.Type == "ORIGIN_DXDYDZ":
                    if a.Type == "Vector":
                        data = data.transpose(2,1,0,3)
                    else:
                        data = data.T
                res.append(np.copy(data))
        return res ;


    def GetFieldData(self,name):
        for a in self.attributes:
            if a.Name == name :
                data = a.dataitems[0].GetData();
                print(data.shape)
                if self.geometry.Type == "ORIGIN_DXDYDZ":
                    if a.Type == "Vector":
                        return data.transpose(2,1,0,3)
                    else:
                        return data.transpose(2,1,0)
                return data
        raise FieldNotFound(name)

    def GetFieldTermsAsColMatrix(self,fieldname, offset = 0):
        import numpy as __np
        #first we check the padding max 8 zeros of padding
        padding = 0
        for i in range(8):
            padding = i
            ss = fieldname+str(offset).zfill(padding)
            #print ss
            if( self.HasField(ss)):
                break
            #print(padding)
        else:
            raise FieldNotFound(fieldname)
        #first we check the number of terms
        #print('padding ' + str(padding))
        cpt =0;
        while(self.HasField(fieldname+ str(offset+cpt).zfill(padding) )):
            cpt +=1;
        # get the firt data to get the type and the size
        d_0 = self.GetFieldData(fieldname+str(offset).zfill(padding) )

        # now we allocate the np matrix for the data
        res = __np.empty([d_0.size, cpt], dtype=d_0.dtype);
        res[:,0] = d_0.reshape(d_0.size);

        for i in range(1,cpt):
            res[:,i] = self.GetFieldData(fieldname+str(offset+i).zfill(padding) ).reshape(d_0.size)

        return res
    def GetFieldTermsAsTensor(self,fieldname,sep='_',offseti=0,offsetj=0):
        import numpy as __np
        from itertools import product
        #first we check the padding max 8 zeros of padding
        paddingi = 0
        paddingj = 0
        for i,j in product(range(8),range(8)):
            paddingi = i
            paddingj = j
            ss = fieldname + str(offseti).zfill(paddingi) + sep + str(offsetj).zfill(paddingj)
            if( self.HasField(ss)):
                break
        else:
            raise FieldNotFound(fieldname+"*"+ str(offseti)+ sep + "*" +str(offsetj))


        #first we check the number of terms
        cpti =0;
        while(self.HasField(fieldname+str(offseti+cpti).zfill(paddingi) + sep + str(offsetj).zfill(paddingj) )):
           cpti +=1;

        cptj =0;
        while(self.HasField(fieldname+str(offseti).zfill(paddingi)  + sep + str(offsetj+cptj).zfill(paddingj) )):
            cptj +=1;

        # get the firt data to get the type and the size
        d_0_0 = self.GetFieldData(fieldname+str(offseti).zfill(paddingi)  + sep + str(offsetj).zfill(paddingj) )

        # now we allocate the np matrix for the data
        res = __np.empty([d_0_0.size, cpti, cptj], dtype=d_0_0.dtype);
        for i in range(0,cpti):
            for j in range(0,cptj):
                #print(fieldname + str(offseti+i).zfill(paddingi)  + sep + str(offsetj+j).zfill(paddingj))
                res[:,i,j] = self.GetFieldData(fieldname + str(offseti+i).zfill(paddingi)  + sep + str(offsetj+j).zfill(paddingj) ).reshape(d_0_0.size);
        return res

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfGrid\n'
        TFormat.II()
        for i in self.informations:
            res += i.__str__();
        res += TFormat.GetIndent() +'Name : '+ self.Name +  '\n'
        res += self.geometry.__str__();
        res += self.topology.__str__();
        for a in self.attributes:
            res += a.__str__();
        TFormat.DI()
        return res

class XdmfInformation(Xdmfbase):
    """ class for extra information in the xmdf file"""
    def __init__(self):
        self.Name = '';
        self.Value = '';

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
        import numpy as np

        self.Type = self.ReadAttribute(attrs,'Type', default=-1)
        if self.Type is -1:
            self.Type = self.ReadAttribute(attrs,'TopologyType')

        if self.Type != "Mixed":
            self.Dimensions = np.array(self.ReadAttribute(attrs,'Dimensions').split(), dtype='int')[::-1]

    def GetConnectivity(self):

        if self.Type != "3DCoRectMesh":
            return self.dataitems[0].GetData()
        else:
            raise Exception# pragma: no cover


    def GetDimensions(self):
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
        self.Type = self.ReadAttribute(attrs,'Type')

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
        if self.Type == "XYZ":
            return self.dataitems[0].GetData()
        else:
            raise Exception# pragma: no cover



class XdmfAttribute(Xdmfbase):
    """  XdmfAttribute class: to store the data over the grids """

    def __init__(self):
        self.dataitems= []
        self.Name = '';
        self.Type =  '';
        self.Center = '';
        self.CDATA = '';

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
            res += d.__str__();
        TFormat.DI();
        return res

class XdmfTime(Xdmfbase):
    """  XdmfTime class: to store the data over the grids """

    def __init__(self):
        self.Value = None;

    def ReadAttributes(self,attrs):
        self.Value = np.array(self.ReadAttribute(attrs,'Value').split(), dtype=np.float)

    def __str__(self):
        res = TFormat.GetIndent() + 'XdmfTime'
        TFormat.II()
        res +=TFormat.GetIndent() + 'Value : '+ str(self.Value) +  '\n'
        TFormat.DI();
        return res

class XdmfDataItem(Xdmfbase):
    """ XdmfDataItem class : class to manage the reading of the data Heavy and light """

    def __init__(self):
        self.Type = None
        self.Dimensions= [];
        self.Precision = None;
        self.Format = None
        self.Data = [];
        self.CDATA = '';
        self.Seek = 0;

    def ReadAttributes(self,attrs, path):
        import numpy as __np
        self.path = path
        self.Dimensions = __np.array(self.ReadAttribute(attrs,'Dimensions').split(), dtype='int')
        try :
            self.Type = self.ReadAttribute(attrs,'DataType')
        except:
            self.Type = self.ReadAttribute(attrs,'NumberType')
        if(self.Type.lower() == 'float' and 'Precision' in attrs):
            self.Precision = int(self.ReadAttribute(attrs,'Precision'))
        self.Format = self.ReadAttribute(attrs,'Format')
        self.Seek = int(self.ReadAttribute(attrs,'Seek',0))

    def __str__(self):
        res = TFormat.GetIndent() +'XdmfDataItem \n'
        TFormat.II();
        res += TFormat.GetIndent() + 'Type : '+ self.Type +  '\n'
        res += TFormat.GetIndent() + 'Dimensions : '+ str(self.Dimensions) +  '\n'
        if(self.Type.lower() == 'float'):
            res += TFormat.GetIndent() + 'Precision : '+ str(self.Precision) +  '\n'
        res += TFormat.GetIndent() + 'Format : '+ str(self.Format) +  '\n'
        res += TFormat.GetIndent() + 'Data : \n  '+  TFormat.GetIndent()  +str(self.Data) +  '\n'
        res += TFormat.GetIndent() + 'CDATA : \n  '+  TFormat.GetIndent()  +str(self.CDATA) +  '\n'
        TFormat.DI();
        return res

    def TreatCDATA(self):
        import numpy as np
        import os as os
        if len(self.CDATA) == 0:
            return ;


        if( self.Format  == 'XML'):

          if(self.Type.lower() =='float'):
              if(self.Precision == 4):
                  numpytype = 'float32'
              else:
                  numpytype = 'float_'
          elif(self.Type.lower() =='int'):
              numpytype = 'int_'

          self.Data= np.array(self.CDATA.split(), dtype=numpytype);
          self.Data = self.Data.reshape(self.Dimensions)
          self.CDATA = '';
        elif ( self.Format  == 'HDF'):

            filename,dataSetPath  = str(self.CDATA).lstrip().rstrip().split(":")
            #print(dataSetPath)
            from h5py import File as __File
            f = __File(os.path.join(self.path, filename),'r')
            #print(f[dataSetPath])
            self.Data =  np.array(f[dataSetPath])
            #print(self.Data)
            self.CDATA = '';

        elif ( self.Format  == 'Binary'):
            if(self.Type.lower() =='float'):
                if(self.Precision == 4):
                  numpytype = 'float32'
                else:
                  numpytype = 'float_'
            elif(self.Type.lower() =='int'):
                numpytype = 'int_'

            binfilename  = str(self.CDATA).lstrip().rstrip()
            binfile = open (os.path.join(self.path, binfilename ), "rb")
            binfile.seek(self.Seek)

            self.Data = np.fromfile(binfile, dtype=numpytype, count=np.prod(self.Dimensions))
            self.Data.shape = self.Dimensions
            #print(self.Data.shape)
            binfile.close()
            self.CDATA = '';
        else :
            raise Exception("Heavy data in format '" + self.Format + "' not suported yet") # pragma: no cover

    def GetData(self):
        self.TreatCDATA()
        return self.Data


import xml.sax
class XdmfReader(xml.sax.ContentHandler):

    def __init__(self,filename=''):
        self.xdmf = Xdmf();
        self.pile = [];
        self.path = '';
        self.filename = '';
        self.SetFileName(filename);
        self.readed = False;
        self.lazy = True;

    def Reset(self):
        self.xdmf = Xdmf();
        self.pile = [];
        self.readed = False;

    def SetFileName(self,filename):
        #if same filename no need to read the file again
        if  self.filename == filename: return
        import os as __os

        self.readed = False;
        self.filename = filename;
        self.path = __os.path.dirname(filename)


    def Read(self):
        if len(self.filename) == 0 :
            raise Exception('Need a filename ')

        # read only one time the file
        if( self.readed ): return;
        self.Reset();

        thefile = open(self.filename,"r")
        # to deactivate the DTD validation
        parser = xml.sax.make_parser()
        parser.setContentHandler(self)
        parser.setFeature(xml.sax.handler.feature_external_ges, False)
        parser.parse(thefile)
        thefile.close();

    # this a a overloaded function (must start with lower case)      !!!!!
    def startElement(self, name, attrs):

        if name == "Xdmf":
            self.pile.append(self.xdmf)
            return

        father = self.pile[-1];

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
            father.topology = res;
        elif name == 'Geometry':
            res = XdmfGeometry()
            res.ReadAttributes(attrs)
            father.geometry = res;
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

        else:
            raise Exception("Unkown tag :  '"+ name +"' Sorry!") # pragma: no cover

        self.pile.append(res)

    # this a a overloaded function (must start with lower case)      !!!!!
    def characters(self, content):
        #print("*** " + content + " +++ ")
        father = self.pile[-1];
        if isinstance(father,XdmfDataItem):
            father.CDATA += content  # Note: '+=', not '='

    # this a a overloaded function (must start with lower case)      !!!!!
    def endElement(self, name):
        if self.lazy :
            self.pile.pop();
        else:
            self.pile.pop().TreatCDATA();

    def __str__(self):
        res = ''
        for d in self.xdmf.domains:
            res = d.__str__() + "\n";
        return res

def GetTensorRepOfField(domain,fieldname):
    """ Get The tensor representation of a field in a xdmf domain

    Get The tensor representation of a field in a xdmf domain, for the moment
    works for Cannonic and Train Tensor formats
    """
    import BasicTools.T.Formats as st
    # we check the nature of the field (info in the parent domain)
    fieldtype = 'CP'
    for info in domain.informations:
        if info.Name == (fieldname+"_Type"):
            fieldtype = info.Value

    if fieldtype == "CP":
        res = st.CanonicTensor()
        for g in  domain.grids:
            G = st.FullTensor(nRanks=2);
            G.SetRanksNames(["-".join([i.Value for i in g.informations if not(i.Name == 'Dims')]) ,'alpha1'])
            G.array = g.GetFieldTermsAsColMatrix(fieldname+"_")
            res.AddSubTensor(G)
    elif fieldtype == "TT":
        res = st.TensorTrain()
        l = len(domain.grids)
        for cpt,g in  zip(range(0,l),domain.grids):
            if cpt == 0:
                G = st.FullTensor(nRanks=2);
                G.SetRanksNames(["-".join([i.Value for i in g.informations if not(i.Name == 'Dims')]) ,'alpha1'])
                G.array = g.GetFieldTermsAsColMatrix(fieldname+"_")
            elif (cpt == l-1):
                G = st.FullTensor(nRanks=2);
                G.SetRanksNames(["-".join([i.Value for i in g.informations if not(i.Name == 'Dims')]) ,'alpha'+str(cpt)])
                G.array = g.GetFieldTermsAsColMatrix(fieldname+"_")
            else:
                G = st.FullTensor(nRanks=3);
                G.SetRanksNames(["-".join([i.Value for i in g.informations if not(i.Name == 'Dims')]) ,'alpha'+str(cpt+1),'alpha'+str(cpt)])
                G.array = g.GetFieldTermsAsTensor(fieldname+"_")
            res.AddSubTensor(G)
    return res


def CheckIntegrity():


    # All the try pasrt are to execute the code for error handling.
    # The code in the try(s) must raise exceptions
    # TODO Verified if the exceptions are the good ones

    # Get TestData directory
    import os
    TestDataPath = os.path.dirname(os.path.abspath(__file__))+os.sep + '..' + os.sep +'TestData' + os.sep

    # Create a Reader
    res = XdmfReader()
    try :
        res.Read();
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    res = XdmfReader(filename = TestDataPath + "Unstructured.xmf" )
    res.lazy = False;
    res.Read();

    res = XdmfReader(filename = TestDataPath + "Unstructured.xmf" )
    # read only the xml part
    res.Read();


    rep = res.__str__()
    #print(rep)
    #Test xdmf part***********************
    a = res.xdmf.__str__();

    try:
        res.xdmf.GetDomain("Imaginary Domain");
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    res.xdmf.GetDomain("Dom 1");

    domain = res.xdmf.GetDomain(0)
    #Test domain part***********************
    try :
        domain.GetGrid("imaginary grid");
        return 'Not OK' # pragma: no cover
    except Exception as e :
        pass

    domain.GetGrid(0);
    grid  = domain.GetGrid("Grid")


    #Test Grid part**************************************************************


    names = grid.GetFieldsNames()
    if(names[0] != 'RTData'): raise Exception();

    grid.HasField(names[0])
    grid.HasField('toto')

    try:
        data = res.xdmf.GetDomain("Dom 1").GetGrid("Grid").GetFieldData('ImaginaryField');
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


    data = res.xdmf.domains[0].GetGrid("Grid").GetFieldData('RTData');
    if( data[49] != 260.0 ): raise ;

    geo = res.xdmf.domains[0].GetGrid("Grid").geometry.dataitems[0].GetData()
    #print(geo)
    if( geo[0,2] != -1 ): raise ;

    topo = res.xdmf.domains[0].GetGrid("Grid").topology.dataitems[0].GetData()
    #print(topo)
    if( topo[2] != 1 ): raise ;

    ######################### Structured #########################
    res = XdmfReader(filename = TestDataPath + "Structured.xmf" )
    # read only the xml part
    res.Read();
    domain = res.xdmf.GetDomain(0)
    #Test domain part***********************

    grid  = domain.GetGrid(0)
    print(grid.GetFieldsOfType("Node"))
    print("--")
    print(grid.GetFieldsOfType("Cell"))
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
    #res.lazy = False;
    reader.Read();

    # Get the domaine "Dom 1"
    dom = reader.xdmf.GetDomain("Dom 1");

    # Get the first Grid
    grid = dom.GetGrid(0)
    grid.topology.GetDimensions()
    names = grid.GetPointFieldsNames()
    names = grid.GetCellFieldsNames()
    names = grid.GetGridFieldsNames()
    allFields = grid.GetPointFields()
    allFields = grid.GetCellFields()
    allFields = grid.GetGridFields()
    #Get one field (or term)
    dataField1= grid.GetFieldData('RTData');

    dataField1= grid.GetFieldsOfType('Node');
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
    #res.lazy = False;
    reader.Read();


    # Get the domaine "Dom 1"
    dom = reader.xdmf.GetDomain(0);

    CT = GetTensorRepOfField(dom,'CanonicField')

    TT = GetTensorRepOfField(dom,'TTField')
    TT.CheckIntegrity(verbose=False)
    #print(TT)

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover
    a = XdmfReader("/data/fbordeu/conforme/Bar.xmf")
    a.Read()
    grid = a.xdmf.GetDomain(0).GetGrid(-1)
    #print(grid.GridType)
    print(grid.geometry.GetNodes()[0,:])
    print(grid.topology.GetConnectivity())

