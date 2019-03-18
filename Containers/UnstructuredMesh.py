# -*- coding: utf-8 -*-
import numpy as np


import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tag
from BasicTools.Containers.MeshBase import Tags

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject,froze_it


AllElements = object()

@froze_it
class ElementsContainer(BaseOutputObject):
    """
    Class to hold a list of element of the same type

    elementType : a string form BasicTools.Containers.ElementNames
    connectivity : the connectivity matrix starting form 0
    tags : the tags holder class

    originalIds : the id or number from the previous mesh/file
    originalOffset : the offset from the previous mesh/file
    the user cans use this data to find the mapping from the inintial mesh/file
    to the currect mesh (self)


    self.globaloffset  : this value is calculate automaticaly by the mesh
    self.cpt : an internal counter to do efficient add of elements one by one

    the user is responsible to call self.tighten() to compact the connectivity
    matrix after the population ( calls AddNewElement(...) or allocate(...))

    """
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
        """
        return the number of elements in this container
        """
        return self.cpt
        #return self.connectivity.shape[0]

    def Merge(self,other,offset=None):
        """
        Merge the elements from the other container into this.

        Non elimination of double elements is done

        if an offset is supplied the connectivity of the other container is
        shifted by the value of the offset during the merge
        """
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
        """
        append a new element to the connectivity

        inputs:
        conn : connectivity of the added element
        originalid : the original id of the added element

        return the total number of elements in the container
        """
        if self.cpt >= self.connectivity.shape[0]:
            self.Reserve(2*self.cpt+1)

        self.connectivity[self.cpt,:] = conn
        self.originalIds[self.cpt] = originalid
        self.cpt +=1

        return self.cpt

    def GetNumberOfNodesPerElement(self):
        """
        return the number of nodes per element for the elements in this container
        """
        if  self.connectivity.shape[1] : return self.connectivity.shape[1]
        return ElementNames.numberOfNodes[self.elementType]

    def GetNodesIdFor(self,ids):
        """
        return the nodes used by the list of elements

        input:
            ids : list of ids of element to treat (always a local id list)
        """
        return np.unique(self.connectivity[ids,:])

    def GetTag(self, tagName):
        """
        return the tag based is a name
        if the tag does not exist a new tag is created
        """
        return self.tags.CreateTag(tagName,False)

    def Reserve(self,nbElements):
        """
        Reserve the storage for nbElements

        the user is responsible to call self.tighten() to compact the connectivity
        matrix after the population

        """
        if (nbElements != self.connectivity.shape[0]):
            self.connectivity =  np.resize(self.connectivity, (nbElements,self.GetNumberOfNodesPerElement()))
            self.originalIds =  np.resize(self.originalIds, (nbElements,))

            #self.connectivity = np.empty((nbElements,self.GetNumberOfNodesPerElement()),dtype=np.int)
            #self.originalIds = np.empty((nbElements,),dtype=np.int)
    def Allocate(self,nbElements):
        """
        Allocate the storage for nbElements

        the user is responsible of filling the connectivity and the originalid
        with valid values

        """
        self.Reserve(nbElements)
        self.cpt = nbElements

    def tighten(self):
        """
        to compact the storage an free non used space
        """
        self.Reserve(self.cpt)


    def AddElementToTag(self,globalElemNumber,tagname):
        """
        Add an element to a tag using a global element number
        The user must compute the globaloffset first to make this function work

        """
        if globalElemNumber -self.globaloffset <  self.GetNumberOfElements():
            self.tags.CreateTag(tagname,False).AddToTag(globalElemNumber-self.globaloffset)
            return 1
        else:
            return 0

    def __str__(self):
        res  = "ElementsContainer \n"
        res += "  elementType    : {}\n".format(self.elementType)
        res += "  Number Of Elements : {}\n".format(self.GetNumberOfElements())
        res += "  Tags          : " + str([x.name for x in self.tags]) + "\n"
        return res

@froze_it
class AllElements(object):
    """
    Class to store a list of element containers
    This class is a sorted by keys dictioniary to keep all the always in order

    note:
      FB: the number of different types of elements is low, I don't think
      this is gonna add alot of overhead to the library
    """

    def __init__(self):
        super(AllElements,self).__init__()
        self.storage = {}

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

    def GetElementsOfType(self,typename):
        if not typename in self:
            self[typename] = ElementsContainer(typename)
        return self[typename]

class UnstructuredMesh(MeshBase):
    """
    class to store a UnstructuredMesh:
        self.nodes : the points positions
        self.orignilaIdNodes : the ids of the previous mesh/file
        self.elements : the list of all the elememnt in the mesh
        self.boundingMin/Max : the bounding box of the mesh (use ComputeBoundingBox
         to compute it)
    """
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
        """
        return the total number of nodes in the mesh
        """
        return self.nodes.shape[0]

    def GetDimensionality(self):
        """
        return the dimensionality 2/3
        """
        return self.nodes.shape[1]

    def GetNumberOfElements(self,dim = None):
        """
        Compute and return the total number of elements in the mesh
        """
        n = 0
        for elemname, data in self.elements.items():
            if dim == None:
                n += data.GetNumberOfElements()
            else:
                if ElementNames.dimension[elemname] == dim:
                    n += data.GetNumberOfElements()
        return n



    def MergeElements(self,other,force=False):
        """
        Merge the element for a second mesh into this
        the nodes array must be the same (not only equal)

        the user can force the merge if needed (force variable)
        """
        if (self.nodes is not other.nodes) and (not force) :
            raise(RuntimeError("the two meshes does not share the same nodes fiels (potentially dangerous)"))

        for name,data in other.elements.items():
            self.GetElementsOfType(name).Merge(data)

    def ComputeGlobalOffset(self):
        """
        Recompute the Global Offset,
        This is necessary for some operation.
        Recomendation : Call it after changing the topology
        """

        cpt = 0
        for type, data in self.elements.items():
            data.globaloffset = cpt
            n = data.GetNumberOfElements()
            cpt = cpt + n

    def GetElementsOfType(self,typename):
        """
        return the element container for the typename element
        """
        return self.elements.GetElementsOfType(typename)

    def ComputeBoundingBox(self):
        """
        to recumpute the bounding box
        """
        self.boundingMin = np.amin(self.nodes, axis=0);
        self.boundingMax = np.amax(self.nodes, axis=0);

    def AddNodeToTagUsingOriginalId(self,oid,tagname):
        """
        add a node (using the original id ) to a tag (tagname)
        """
        w = np.where(self.originalIDNodes == oid)
        if len(w[0]) > 0 :
            self.GetNodalTag(tagname).AddToTag(w[0])
        else:
            raise Exception("No node with id " + str(oid)) #pragma: no cover


    def AddElementToTagUsingOriginalId(self,oid,tagname):
        """
        add a element (using the originalid) to a tag (tagname)
        """
        for ntype, data in self.elements.items():
            w = np.where(data.originalIds[:data.cpt] == oid)
            if len(w[0]) > 0 :
                data.tags.CreateTag(tagname,False).AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid)) #pragma: no cover


    def AddElementToTag(self,globalElemNumber,tagname):
        """
        add a element (using the global element number) to a tag (tagname)
        # you must compute the globaloffset first to make this function work
        """
        for ntype, data in self.elements.items():
            if data.AddElementToTag(globalElemNumber,tagname):
                return
        raise Exception("No element found") #pragma: no cover

    def DeleteElemTags(self, tagNames):
        """
        delete element tags
        """
        #check not a string but a list like
        assert not isinstance(tagNames, str)
        for ntype, data in self.elements.items():
            data.tags.DeleteTags(tagNames)

    def GetPosOfNode(self, i ):
        """
        return the position of the point i
        """
        return self.nodes[i,:]

    def GetPosOfNodes(self):
        """
        return the position of all the nodes
        """
        return self.nodes

    def GetNamesOfElemTags(self):
        """
        return a list containing all the element tags present in the mehs
        """
        res = set()
        for ntype, data in self.elements.items():
            for tag in data.tags:
                res.add(tag.name)

        return list(res)

    def GetElementsOriginalIDs(self):
        """
        return a single list with all the originalid concatenated
        """
        res = np.empty(self.GetNumberOfElements(),dtype=np.int)
        cpt = 0
        for ntype, data in self.elements.items():
            res[0:data.GetNumberOfElements()] = data.originalIds+data.originalOffset
            cpt += data.GetNumberOfElements()
        return res

    def GetElementsInTag(self,tagname,useOriginalId=False) :
        """
        return a list with the ids of the elements in a tag
        """
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
        """
        function to free the extra memory used during a incremental creation of a mesh
        and final treatement (offset computation)
        """
        self.ComputeGlobalOffset()
        for ntype, data in self.elements.items():
             data.tighten()
             for tag in data.tags:
                   tag.Tighten()

    def GenerateManufacturedOriginalIDs(self):
        """
        function to generate a valid originalid data
        """
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
            res += "  nodeFields         : " + str(list(self.nodeFields.keys())) + "\n"
        if len(self.elemFields.keys()):
            res += "  elemFields         : " + str(list(self.elemFields.keys())) + "\n"
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
    print(resII.GetNumberOfElements(dim=2))
    del resII.elements[ElementNames.Triangle_3]
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
