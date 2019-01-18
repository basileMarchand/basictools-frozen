# -*- coding: utf-8 -*-
import numpy as np


import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tag
from BasicTools.Containers.MeshBase import Tags

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

AllElements = object()

class ElementsContainer(BaseOutputObject):
    def __init__(self,elementType):
        super(ElementsContainer,self).__init__()
        self.elementType = elementType
        self.connectivity = np.empty((0,0),dtype=np.int)
        self.globaloffset   = 0
        self.tags = Tags()
        self.cpt = 0;

        self.originalIds = np.empty((0,),dtype=np.int)
        self.originalOffset = 0

    def GetNumberOfElements(self):
        return self.cpt
        #return self.connectivity.shape[0]

    def Merge(self,other,offset=None):
        other.tighten()

        self.Reserve(self.cpt+other.cpt)

        if offset is None:
            offset = 0

        self.connectivity[self.cpt:,:] = other.connectivity+offset
        self.originalIds[self.cpt:] = -1*np.arange(other.cpt)

        for tag in other.tags:
            self.GetTag(tag.name).AddToTag(tag.GetIds() + self.cpt)

        self.cpt += other.cpt

    def AddNewElement(self,conn,originalid):
        if self.cpt >= self.connectivity.shape[0]:
            self.Reserve(2*self.cpt+1)

        self.connectivity[self.cpt,:] =conn
        self.originalIds[self.cpt] =originalid
        self.cpt +=1

        #if(self.connectivity.shape[1] == 0):
        #    self.connectivity = np.empty((0,len(conn)),dtype=np.int)
        #self.connectivity = np.vstack((self.connectivity,np.array(conn,dtype=int) ))
        #self.originalIds = np.append(self.originalIds,originalid)

        return self.cpt

    def GetNumberOfNodesPerElement(self):
        if  self.connectivity.shape[1] : return self.connectivity.shape[1]
        return ElementNames.numberOfNodes[self.elementType]

    def GetNodesIdFor(self,ids):
        return np.unique(self.connectivity[ids,:])

    def GetTag(self, tagName):
        return self.tags.CreateTag(tagName,False)

    def Reserve(self,nbElements):
        if (nbElements != self.connectivity.shape[0]):
            self.connectivity =  np.resize(self.connectivity, (nbElements,self.GetNumberOfNodesPerElement()))
            self.originalIds =  np.resize(self.originalIds, (nbElements,))

            #self.connectivity = np.empty((nbElements,self.GetNumberOfNodesPerElement()),dtype=np.int)
            #self.originalIds = np.empty((nbElements,),dtype=np.int)
    def Allocate(self,nbElements):
        self.Reserve(nbElements)
        self.cpt = nbElements

    def tighten(self):
        self.Reserve(self.cpt)

    # you must compute the globaloffset first to make this function work
    def AddElementToTag(self,elemNumber,tagname):
        if elemNumber-self.globaloffset <  self.GetNumberOfElements():
            self.tags.CreateTag(tagname,False).AddToTag(elemNumber-self.globaloffset)
            return 1
        else:
            return 0

    def DeleteElementsById(self,ids):
        mask = np.ones(self.GetNumberOfElements(),dtype=np.bool)
        mask[ids] = False
        from BasicTools.Containers.UnstructuredMeshTools import ExtractElementsByMask
        return  ExtractElementsByMask(self,mask)

    def __str__(self):
        res  = "ElementsContainer \n"
        res += "  elementType    : {}\n".format(self.elementType)
        res += "  Number Of Elements : {}\n".format(self.GetNumberOfElements())
        res += "  Tags          : " + str([x.name for x in self.tags]) + "\n"
        return res

class AllElements(object):
    def __init__(self):
        super(AllElements,self).__init__()
        self.storage = {}

    ## to make all the key always in order
    ## the number of different types of elements is reduced so, I don't think
    ## this is gonna add alot of overhead to the library
    def keys(self):
        return sorted(self.storage.keys())

    def __iter__(self):
        return iter(self.keys())

    def items(self):
        return  sorted(self.storage.items())

    #send basis functions calls to the storage dictionary
    def __setitem__(self, key, value):
            self.storage[key] = value

    def __len__(self):
        return len(self.storage)

    def __contains__(self, k):
        return k in self.storage

    def __getitem__(self,key):
        return self.storage[key]

    def __delitem__(self,key):
        del self.storage[key]


class UnstructuredMesh(MeshBase):

    def IsUnstructured(self):
        return True

    def __init__(self):
        super(UnstructuredMesh,self).__init__()
        self.nodes = np.empty((0,3),dtype=np.double)
        self.originalIDNodes = np.empty((0,),dtype=np.int)
        self.elements = AllElements();
        self.boundingMin = [0,0,0];
        self.boundingMax = [0,0,0];

    def GetNumberOfNodes(self):
        return self.nodes.shape[0]

    def GetDimensionality(self):
        return self.nodes.shape[1]

    def GetNumberOfElements(self,dim = None):
        n = 0
        for elemname, data in self.elements.items():
            if dim == None:
                n += data.GetNumberOfElements()
            else:
                if ElementNames.dimension[elemname] == dim:
                    n += data.GetNumberOfElements()

        return n

    def MergeElements(self,other,force=False):
        if (self.nodes is not other.nodes) and (not force) :
            raise(RuntimeError("the two meshes does not share the same nodes fiels (potentially dangerous)"))

        for name,data in other.elements.items():
            self.GetElementsOfType(name).Merge(data)

    def ComputeGlobalOffset(self):
        cpt = 0
        for type, data in self.elements.items():
            data.globaloffset = cpt
            n = data.GetNumberOfElements()
            cpt = cpt + n

    def GetElementsOfType(self,typename):
        if not typename in self.elements:
            self.elements[typename] = ElementsContainer(typename)
        return self.elements[typename]

    def ComputeBoundingBox(self):
        self.boundingMin = np.amin(self.nodes, axis=0);
        self.boundingMax = np.amax(self.nodes, axis=0);

    def AddNodeToTagUsingOriginalId(self,oid,tagname):

        w = np.where(self.originalIDNodes == oid)
        if len(w[0]) > 0 :
            self.GetNodalTag(tagname).AddToTag(w[0])
        else:
            raise Exception("No node with id " + str(oid)) #pragma: no cover


    def AddElementToTagUsingOriginalId(self,oid,tagname):
        for ntype, data in self.elements.items():
            w = np.where(data.originalIds[:data.cpt] == oid)
            if len(w[0]) > 0 :
                data.tags.CreateTag(tagname,False).AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid)) #pragma: no cover


    # you must compute the globaloffset first to make this function work
    def AddElementToTag(self,elemNumber,tagname):
        for ntype, data in self.elements.items():
            if data.AddElementToTag(elemNumber,tagname):
                return
        raise Exception("No element found") #pragma: no cover

    def DeleteElemTags(self, tagNames):
        for ntype, data in self.elements.items():
            data.tags.DeleteTags(tagNames)


    def GetPosOfNode(self, i ):
        return self.nodes[i,:]

    def GetPosOfNodes(self):
        return self.nodes

    def GetNamesOfElemTags(self):
        res = set()
        for ntype, data in self.elements.items():
            for tag in data.tags:
                res.add(tag.name)

        return list(res)
#    def GetElementTag(self,tagname):
#        ne = self.GetNumberOfElements()
#        res = np.zeros((ne,1),dtype=np.int)
#        for ntype, elem in self.elements.iteritems():
#            if elem.tags.has_key(tagname):
#                print(elem.tags[tagname].id.shape)
#                res[elem.globaloffset+elem.tags[tagname].id] = 1;
#        return res

    def GetElementsOriginalIDs(self):
        res = np.empty(self.GetNumberOfElements(),dtype=np.int)
        cpt = 0
        for ntype, data in self.elements.items():
            res[0:data.GetNumberOfElements()] = data.originalIds+data.originalOffset
            cpt += data.GetNumberOfElements()
        return res

    def GetElementsInTag(self,tagname,useOriginalId=False) :
        self.ComputeGlobalOffset();
        ne = self.GetNumberOfElements()
        res = np.zeros((ne,),dtype=np.int)
        cpt =0
        for ntype, elem in self.elements.items():
            if tagname in elem.tags:

                tag = elem.tags[tagname].GetIds()
                if useOriginalId:
                    res[cpt:cpt+len(tag) ] = elem.originalIds[tag];
                else:
                    res[cpt:cpt+len(elem.tags[tagname].GetIds()) ] = elem.globaloffset+tag;
                cpt +=  len(tag)
        return res[0:cpt]

    def PrepareForOutput(self):
       self.ComputeGlobalOffset()
       for ntype, data in self.elements.items():
             data.tighten()
             for tag in data.tags:
                   tag.Tighten()

    def GenerateManufacturedOriginalIDs(self):
       self.originalIDNodes = np.arange(self.GetNumberOfNodes())
       counter = 0
       for key, value in self.elements.items():
           value.originalIds = np.arange(counter,counter+value.GetNumberOfElements())
           counter += value.GetNumberOfElements()


    def __str__(self):
        res  = "UnstructuredMesh \n"
        res += "  Number Of Nodes    : {}\n".format(self.GetNumberOfNodes())
        res += "  Number Of Elements : {}".format(self.GetNumberOfElements())
        for name,data in self.elements.items():
            if data.GetNumberOfElements():
                res += " ({}:{})".format(name,data.GetNumberOfElements())
        res += "\n"
        res += "  Node Tags          : " + str(self.nodesTags) + "\n"
        res += "  Cell Tags          : " + str([x for x in self.GetNamesOfElemTags()])+ "\n"
        if len(self.nodeFields.keys()):
            res += "  nodeFields         : " + str(self.nodeFields.keys()) + "\n"
        if len(self.elemFields.keys()):
            res += "  elemFields         : " + str(self.elemFields.keys()) + "\n"
        return res

def CheckIntegrity():
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOfTriangles
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshFromConstantRectilinearMesh

    res = CreateMeshOfTriangles([[0,0,0],[1,2,3],[0,1,0]], [[0,1,2]])

    elements = res.GetElementsOfType(ElementNames.Triangle_3)

    elements = res.GetElementsOfType(ElementNames.Bar_2)
    elements.AddNewElement([1,2],1)

    elements.GetNumberOfNodesPerElement()

    print(res.IsUnstructured())

    res.ComputeGlobalOffset()
    res.AddElementToTag(1,"SecondElement")

    if res.GetNumberOfElements() != 2: raise Exception()
    res.ComputeGlobalOffset()

    print(res.GetDimensionality())
    res.ComputeBoundingBox()
    print(res.boundingMin)
    print(res.boundingMax)

    res.nodesTags.CreateTag('toto')
    print(res.GetNodalTag('toto'))
    print(res.GetNodalTag('toto2'))

    res.AddElementToTagUsingOriginalId(1,"bars")

    if res.GetPosOfNodes()[1,1] != 2: raise Exception()

    print(res.PrepareForOutput())
    print(res.GetElementsInTag("bars"))
    print(res.GetElementsInTag("bars",useOriginalId=True))

    print(res.GetElementsOfType(ElementNames.Bar_2).GetTag("toto"))
    print(res.GetPosOfNode(0))
    print(res)
    res.DeleteElemTags(["SecondElement"])
    print(res)

    resII = CreateMeshOfTriangles([[0,0,0],[1,2,3],[0,1,0]], [[0,1,2]])
    resII.AddNodeToTagUsingOriginalId(0,"First Point")
    resII.GenerateManufacturedOriginalIDs()

    try:
        resII.MergeElements(res)
        raise#pragma: no cover
    except:
        pass
    resII.MergeElements(res,force=True)
    print("----")
    print(resII)
    print(resII.GetElementsOfType(ElementNames.Triangle_3).DeleteElementsById([0]))
    print(resII.GetNumberOfElements(dim=2))
    del resII.elements[ElementNames.Triangle_3]
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
