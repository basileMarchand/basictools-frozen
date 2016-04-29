# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import MeshBase
from OTTools.FE.MeshBase import Tag

import numpy as np

class ElementsContainer():
    def __init__(self,typename):
        self.elementType = typename
        self.connectivity = np.empty((0,0),dtype=np.int)
        self.globaloffset   = 0
        self.originalIds = np.empty((0,),dtype=np.int)
        self.tags = {}
        
        
    def GetNumberOfElements(self):
        return self.connectivity.shape[0]
    
    def AddNewElement(self,conn,originalid):
        if(self.connectivity.shape[1] == 0):
            self.connectivity = np.empty((0,len(conn)),dtype=np.int)    
        self.connectivity = np.vstack((self.connectivity,conn))
        self.originalIds = np.append(self.originalIds,originalid)
        
    def GetNumberOfNodesPerElement(self):
        return self.connectivity.shape[1]



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
        
    def GetNodalTag(self, tagname):
        if not self.nodesTags.has_key(tagname):
            self.nodesTags[tagname] = Tag(tagname)
        return self.nodesTags[tagname]
        
    def AddElementToTagUsingOriginalId(self,oid,tagname):
        for ntype, data in self.elements.iteritems(): 
            w = np.where(data.originalIds == oid)
            if len(w[0]) > 0 :
                if not data.tags.has_key(tagname):
                    data.tags[tagname] = Tag(tagname)
                data.tags[tagname].AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid))
        
    def GetPosOfNode(self):
        return self.nodes
        
    def GetNamesOfCellTags(self):
        res = set()
        for ntype, data in self.elements.iteritems():
            for tagname in data.tags:
                res.add(tagname)
                
        return res
#    def GetElementTag(self,tagname):
#        ne = self.GetNumberOfElements()
#        res = np.zeros((ne,1),dtype=np.int)
#        for ntype, elem in self.elements.iteritems():
#            if elem.tags.has_key(tagname):
#                print(elem.tags[tagname].id.shape)
#                res[elem.globaloffset+elem.tags[tagname].id] = 1;                 
#        return res
        
    def GetElementsInTag(self,tagname,useOriginalId=False) :
        ne = self.GetNumberOfElements()
        res = np.zeros((ne,),dtype=np.int)
        cpt =0
        for ntype, elem in self.elements.iteritems():
            if elem.tags.has_key(tagname):
                if useOriginalId:
                    res[cpt:cpt+len(elem.tags[tagname].id) ] = elem.originalIds[elem.tags[tagname].id];                 
                else:
                    res[cpt:cpt+len(elem.tags[tagname].id) ] = elem.globaloffset+elem.tags[tagname].id;                 
                cpt +=  len(elem.tags[tagname].id)
        return res[0:cpt]

    def PrepareForOutput(self):
       self.ComputeGlobalOffset()

        
    def __str__(self):
        res  = "UnstructuredMesh \n" 
        res += "  Number Of Nodes : {}\n".format(self.GetNumberOfNodes()) 
        res += "  Number Of Elements : {}\n".format(self.GetNumberOfElements()) 
        res += "  Node Tags : " + str([x for x in self.nodesTags]) + "\n"
        res += "  Cell Tags : " + str([x for x in self.GetNamesOfCellTags()])+ "\n"
        return res