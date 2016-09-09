# -*- coding: utf-8 -*-
from OTTools.FE.MeshBase import MeshBase
from OTTools.FE.MeshBase import Tag
import OTTools.FE.ElementNames as EN
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
        self.connectivity = np.vstack((self.connectivity,np.array(conn,dtype=int) ))
        self.originalIds = np.append(self.originalIds,originalid)
        return self.connectivity.shape[0]

    def GetNumberOfNodesPerElement(self):
        return self.connectivity.shape[1]

    def GetTag(self, tagName):
        if not self.tags.has_key(tagName):
           self.tags[tagName] = Tag(tagName)
        return self.tags[tagName]

def CreateMeshOfTriangles(points,tris):
    res = UnstructuredMesh()
    res
    res.nodes = np.array(points, dtype=np.double)
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int)

    elements = res.GetElementsOfType(EN.Triangle_3)
    elements.connectivity = np.array(tris,dtype=np.int)
    elements.originalIds = np.arange(0,elements.GetNumberOfElements(),dtype=np.int)
    return res

def CreateMeshFromConstantRectilinearMesh(CRM):
    res = UnstructuredMesh()
    res.nodes = CRM.GetPosOfNodes();
    res.originalIDNodes = np.arange(0,res.GetNumberOfNodes(),dtype=np.int);
    #CRM.originalIDNodes

    nbelements = CRM.GetNumberOfElements()

    res.PrintDebug(nbelements)

    if(CRM.GetDimensionality() == 3):
        elements = res.GetElementsOfType(EN.Hexaedron_8)
        elements.connectivity = np.zeros((nbelements,8),dtype=np.int)
    else:
        elements = res.GetElementsOfType(EN.Quadrangle_4)
        elements.connectivity = np.zeros((nbelements,4),dtype=np.int)

    for i in xrange(nbelements):
        elements.connectivity[i,:] = CRM.GetConnectivityForElement(i)
    elements.originalIds = np.arange(0,elements.GetNumberOfElements(),dtype=np.int)

    return res

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
            raise Exception("No element with id " + str(oid)) #pragma: no cover

    def GetPosOfNodes(self):
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
    res = CreateMeshOfTriangles([[0,0,0],[1,2,3]], [[0,2,3]])

    elements = res.GetElementsOfType(EN.Triangle_3)

    elements = res.GetElementsOfType(EN.Bar_2)
    elements.AddNewElement([1,2],1)
    elements.GetNumberOfNodesPerElement()

    print(res.IsUnstructured())

    if res.GetNumberOfElements() != 2: raise Exception()
    res.ComputeGlobalOffset()

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

    ###########################
    from OTTools.FE.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2,2]);
    myMesh.SetSpacing([1, 1, 1]);
    print(CreateMeshFromConstantRectilinearMesh(myMesh))

    ###########################
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([2,2]);
    myMesh.SetSpacing([1, 1]);
    print(CreateMeshFromConstantRectilinearMesh(myMesh))

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
