# -*- coding: utf-8 -*-
import numpy as np

from BasicTools.FE.MeshBase import MeshBase
from BasicTools.FE.MeshBase import Tag
from BasicTools.FE.MeshBase import Tags
import BasicTools.FE.ElementNames as ElementNames
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject


class ElementsContainer(BaseOutputObject):
    def __init__(self,elementType):
        self.elementType = elementType
        self.connectivity = np.empty((0,0),dtype=np.int)
        self.globaloffset   = 0
        self.originalIds = np.empty((0,),dtype=np.int)
        self.tags = Tags()
        self.cpt = 0;

    def GetNumberOfElements(self):
        return self.cpt
        #return self.connectivity.shape[0]

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

    def __str__(self):
        res  = "ElementsContainer \n"
        res += "  elementType    : {}\n".format(self.elementType)
        res += "  Number Of Elements : {}\n".format(self.GetNumberOfElements())
        res += "  Tags          : " + str([x.name for x in self.tags]) + "\n"
        return res

class UnstructuredMesh(MeshBase):

    def IsUnstructured(self):
        return True

    def __init__(self):
        super(UnstructuredMesh,self).__init__()
        self.nodes = np.empty((0,3),dtype=np.double)
        self.originalIDNodes = np.empty((0,),dtype=np.int)
        self.elements = {}
        self.boundingMin = [0,0,0];
        self.boundingMax = [0,0,0];
    
    def GetNumberOfNodes(self):
        return self.nodes.shape[0]

    def GetDimensionality(self):
        return self.nodes.shape[1]

    def GetNumberOfElements(self):
        n = 0
        for type, data in self.elements.iteritems():
            n += data.GetNumberOfElements()
        return n

    def ComputeGlobalOffset(self):
        cpt = 0
        for type, data in self.elements.iteritems():
            data.globaloffset = cpt
            n = data.GetNumberOfElements()
            cpt = cpt + n

    def GetElementsOfType(self,typename):
        if self.elements.has_key(typename):
            return self.elements[typename]
        else:
            els = ElementsContainer(typename)
            self.elements[typename] = els
            return els

    def ComputeBoundingBox(self):
        self.boundingMin = np.amin(self.nodes, axis=0);
        self.boundingMax = np.amax(self.nodes, axis=0);

    def AddElementToTagUsingOriginalId(self,oid,tagname):
        for ntype, data in self.elements.iteritems():
            w = np.where(data.originalIds == oid)
            if len(w[0]) > 0 :
                data.tags.CreateTag(tagname,False).AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid)) #pragma: no cover


    # you must compute the globaloffset first to make this function work
    def AddElementToTag(self,elemNumber,tagname):
        for ntype, data in self.elements.iteritems():
            if data.AddElementToTag(elemNumber,tagname):
                return
        raise Exception("No element found") #pragma: no cover

    def DeleteElemTags(self, tagNames):
        for ntype, data in self.elements.iteritems():
            data.tags.DeleteTags(tagNames)


    def GetPosOfNode(self, i ):
        return self.nodes[i,:]

    def GetPosOfNodes(self):
        return self.nodes

    def GetNamesOfElemTags(self):
        res = set()
        for ntype, data in self.elements.iteritems():
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

    def GetElementsInTag(self,tagname,useOriginalId=False) :
        self.ComputeGlobalOffset();
        ne = self.GetNumberOfElements()
        res = np.zeros((ne,),dtype=np.int)
        cpt =0
        for ntype, elem in self.elements.iteritems():
            if elem.tags.has_key(tagname):

                tag = elem.tags[tagname].GetIds()
                if useOriginalId:
                    res[cpt:cpt+len(tag) ] = elem.originalIds[tag];
                else:
                    res[cpt:cpt+len(elem.tags[tagname].GetIds()) ] = elem.globaloffset+tag;
                cpt +=  len(tag)
        return res[0:cpt]

    def PrepareForOutput(self):
       self.ComputeGlobalOffset()
       for ntype, data in self.elements.iteritems():
             data.tighten()
             for tag in data.tags:
                   tag.Tighten()



    def GenerateManufacturedOriginalIDs(self):
       counter = 0
       for key, value in self.elements.iteritems():
         for i in xrange(value.GetNumberOfElements()):
           value.originalIds[i] = counter
           counter += 1

    def __str__(self):
        res  = "UnstructuredMesh \n"
        res += "  Number Of Nodes    : {}\n".format(self.GetNumberOfNodes())
        res += "  Number Of Elements : {}\n".format(self.GetNumberOfElements())
        res += "  Node Tags          : " + str(self.nodesTags) + "\n"
        res += "  Cell Tags          : " + str([x for x in self.GetNamesOfElemTags()])+ "\n"
        return res

def CheckIntegrity():
    from BasicTools.FE.UnstructuredMeshTools import CreateMeshOfTriangles
    from BasicTools.FE.UnstructuredMeshTools import CreateMeshFromConstantRectilinearMesh

    res = CreateMeshOfTriangles([[0,0,0],[1,2,3]], [[0,2,3]])

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
    print(res.GetNodalTag('toto2') )

    res.AddElementToTagUsingOriginalId(1,"bars")

    if res.GetPosOfNodes()[1,1] != 2: raise Exception()

    print(res.PrepareForOutput())
    print(res.GetElementsInTag("bars"))
    print(res.GetElementsInTag("bars",useOriginalId=True))

    print(res.GetElementsOfType("bars").GetTag("toto"))
    print(res.GetPosOfNode(0))
    print(res)
    res.DeleteElemTags(["SecondElement"])
    print(res)
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
