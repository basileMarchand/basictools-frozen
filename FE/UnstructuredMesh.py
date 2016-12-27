# -*- coding: utf-8 -*-
from OTTools.FE.MeshBase import MeshBase
from OTTools.FE.MeshBase import Tag
import numpy as np
import OTTools.FE.ElementNames as ElementNames

class ElementsContainer():
    def __init__(self,typename):
        self.elementType = typename
        self.connectivity = np.empty((0,0),dtype=np.int)
        self.globaloffset   = 0
        self.originalIds = np.empty((0,),dtype=np.int)
        self.tags = {}
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
        if not self.tags.has_key(tagName):
           self.tags[tagName] = Tag(tagName)
        return self.tags[tagName]
    
    def Reserve(self,nbElements):
        if (nbElements != self.connectivity.shape[0]):
            self.connectivity =  np.resize(self.connectivity, (nbElements,self.GetNumberOfNodesPerElement()))
            self.originalIds =  np.resize(self.originalIds, (nbElements,))
            
            #self.connectivity = np.empty((nbElements,self.GetNumberOfNodesPerElement()),dtype=np.int)
            #self.originalIds = np.empty((nbElements,),dtype=np.int)

    def tighten(self):
        self.Reserve(self.cpt)
        
        

class UnstructuredMesh(MeshBase):

    def IsUnstructured(self):
        return True

    def __init__(self):
        super(UnstructuredMesh,self).__init__()
        self.nodes = np.empty((0,3),dtype=np.double)
        self.originalIDNodes =[]
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
                if not data.tags.has_key(tagname):
                    data.tags[tagname] = Tag(tagname)
                data.tags[tagname].AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid)) #pragma: no cover

    def GetPosOfNodes(self):
        return self.nodes

    def GetNamesOfCellTags(self):
        res = set()
        for ntype, data in self.elements.iteritems():
            for tagname in data.tags:
                res.add(tagname)

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
                
                tag = elem.tags[tagname].id
                if useOriginalId:
                    res[cpt:cpt+len(tag) ] = elem.originalIds[tag];
                else:
                    res[cpt:cpt+len(elem.tags[tagname].id) ] = elem.globaloffset+tag;
                cpt +=  len(tag)
        return res[0:cpt]

    def PrepareForOutput(self):
       self.ComputeGlobalOffset()
       for ntype, data in self.elements.iteritems():
             for name, tag in data.tags.iteritems():
                   tag.tighten()

       

    """def GenerateManufacturedOriginalIDs(self):
       counter = 1
       for key, value in self.elements.iteritems():
         for i in xrange(value.originalIds.shape[0]):
           value.originalIds[i] = counter
           counter += 1"""

    def __str__(self):
        res  = "UnstructuredMesh \n"
        res += "  Number Of Nodes    : {}\n".format(self.GetNumberOfNodes())
        res += "  Number Of Elements : {}\n".format(self.GetNumberOfElements())
        res += "  Node Tags          : " + str([x for x in self.nodesTags]) + "\n"
        res += "  Cell Tags          : " + str([x for x in self.GetNamesOfCellTags()])+ "\n"
        return res

def CheckIntegrity():
    from OTTools.FE.UnstructuredMeshTools import CreateMeshOfTriangles
    from OTTools.FE.UnstructuredMeshTools import CreateMeshFromConstantRectilinearMesh
    
    res = CreateMeshOfTriangles([[0,0,0],[1,2,3]], [[0,2,3]])

    elements = res.GetElementsOfType(ElementNames.Triangle_3)

    elements = res.GetElementsOfType(ElementNames.Bar_2)
    elements.AddNewElement([1,2],1)
    elements.GetNumberOfNodesPerElement()

    print(res.IsUnstructured())

    if res.GetNumberOfElements() != 2: raise Exception()
    res.ComputeGlobalOffset()

    print(res.GetDimensionality())
    res.ComputeBoundingBox()
    print(res.boundingMin)
    print(res.boundingMax)

    res.nodesTags['toto'] = Tag('toto')
    print(res.GetNodalTag('toto'))
    print(res.GetNodalTag('toto2') )

    res.AddElementToTagUsingOriginalId(1,"bars")

    if res.GetPosOfNodes()[1,1] != 2: raise Exception()

    print(res.PrepareForOutput())
    print(res.GetElementsInTag("bars"))
    print(res.GetElementsInTag("bars",useOriginalId=True))

    print(res.GetElementsOfType("bars").GetTag("toto"))
    print(res)

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
