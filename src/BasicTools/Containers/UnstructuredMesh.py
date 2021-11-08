# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.MeshBase import MeshBase
from BasicTools.Containers.MeshBase import Tags
from BasicTools.Containers.Filters import ElementFilter

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject, froze_it

@froze_it
class ElementsContainer(BaseOutputObject):
    """
    Class to hold a list of element of the same type

    * elementType : a string form BasicTools.Containers.ElementNames
    * connectivity : the connectivity matrix starting form 0
    * tags : the tags holder class
    * originalIds : the id or number from the previous mesh/file
    * originalOffset : the offset from the previous mesh/file

    The user can use this data to find the mapping from the inintial mesh/file
    to the currect mesh (self).

    * self.globaloffset  : this value is calculate automaticaly by the mesh
    * self.cpt : an internal counter to do efficient add of elements one by one

    The user is responsible to call self.tighten() to compact the connectivity
    matrix after the population ( calls AddNewElement(...) or allocate(...))
    """
    def __init__(self,elementType):
        super(ElementsContainer,self).__init__()
        self.elementType = elementType
        self.connectivity = np.empty((0,ElementNames.numberOfNodes[elementType]),dtype=np.int)
        self.globaloffset   = 0
        self.tags = Tags()
        self.cpt = 0

        self.originalIds = np.empty((0,),dtype=np.int)
        self.originalOffset = 0
        self.mutable = True

    def __eq__(self, other):

        if self.elementType != other.elementType:
            return False

        self.tighten()
        other.tighten()

        if not np.array_equal(self.connectivity, other.connectivity):
            return False

        if not np.array_equal(self.originalIds, other.originalIds):
            return False

        if self.tags != other.tags:
            return False

        if self.originalOffset != other.originalOffset:
            return False

        return True

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
        if other.cpt == 0:
            return

        self.Reserve(self.cpt+other.cpt)

        if offset is None:
            offset = 0

        self.connectivity[self.cpt:,:] = other.connectivity+offset
        self.originalIds[self.cpt:] = -1*np.arange(other.cpt)

        for tag in other.tags:
            self.GetTag(tag.name).AddToTag(tag.GetIds() + self.cpt)

        self.cpt += other.cpt

    def AddNewElements(self,conn,originalids=None):
        """
        append a new element to the connectivity

        inputs:
        conn : connectivity of the added element
        originalid : the original id of the added element

        return the total number of elements in the container
        """
        onoe = self.GetNumberOfElements()
        self.Allocate(onoe+conn.shape[0])

        self.connectivity[onoe:onoe+conn.shape[0],:] = conn

        if originalids is None:
            self.originalIds[onoe:onoe+conn.shape[0]] = -1
        else:
            self.originalIds[onoe:onoe+conn.shape[0]] = originalids

        return self.cpt

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
        if  self.connectivity.shape[1]:
            return self.connectivity.shape[1]
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
        if nbElements != self.connectivity.shape[0]:
            self.connectivity =  np.resize(self.connectivity, (nbElements,self.GetNumberOfNodesPerElement()))
            self.originalIds =  np.resize(self.originalIds, (nbElements,))

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
        self.tags.Tighten()


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
        res  = "    ElementsContainer, "
        res += "  Type : ({},{}), ".format(self.elementType,self.GetNumberOfElements())
        res += "  Tags : " + " ".join([ ("("+x.name+":"+str(len(x)) +")") for x in self.tags]) + "\n"
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

    def __eq__(self, other):
        if len(self.storage) != len(other.storage):
            return False

        for i in self:
            if i in other:
                if self[i] != other[i]:
                    return False
            else:
                return False
        return True

    def keys(self):
        return sorted(self.storage.keys())

    def values(self):
        return [ values for key,values in sorted(self.storage.items()) ]

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

    def GetTagsNames(self):
        res = set()
        for data in self.values():
            res.update(data.tags.keys() )

        return list(res)

    def __str__(self):
        res = ""
        for data in self.storage.values():
            res += str(data)
        return res

@froze_it
class UnstructuredMesh(MeshBase):
    """
    Class storing an unstructured (i.e. general) mesh:

    * self.nodes : the points positions
    * self.orignilaIdNodes : the ids of the previous mesh/file
    * self.elements : the list of all the elememnt in the mesh
    * self.boundingMin/Max : the bounding box of the mesh (use ComputeBoundingBox to compute it)

    The manual construction of this class must always end with a call to the
    function.
    """
    def IsUnstructured(self):
        return True

    def __init__(self):
        super(UnstructuredMesh, self).__init__()
        self.nodes = np.empty((0,3),dtype=np.double)
        self.originalIDNodes = np.empty((0,),dtype=np.int)
        self.elements = AllElements()
        self.boundingMin = np.array([0.,0,0])
        self.boundingMax = np.array([0.,0,0])

    def __copy__(self):
        res = UnstructuredMesh()
        res._assign(self)
        res.nodes = self.nodes
        res.originalIDNodes = self.originalIDNodes
        res.elements = self.elements
        res.boundingMin = self.boundingMin
        res.boundingMax = self.boundingMax
        return res

    def __eq__(self, other):
        if not super(UnstructuredMesh,self).__eq__(other):
            return False

        if not np.array_equal(self.nodes,other.nodes):
            return False

        if not np.array_equal(self.originalIDNodes, other.originalIDNodes):
            return False

        if self.elements != other.elements:
            return False

        return True

    def ConvertDataForNativeTreatment(self):
        self.nodes = np.asarray(self.nodes, dtype=float, order="C")
        for data in self.elements.values():
            data.connectivity = np.asarray(data.connectivity, dtype=int, order="C")

    def GetNumberOfNodes(self):
        """
        return the total number of nodes in the mesh
        """
        return self.nodes.shape[0]

    def SetNodes(self, arrayLike, originalIDNodes=None, generateOriginalIDs = False):
        """
        Set nodes and original Ids in the correct internal data
        """
        self.nodes = np.require(arrayLike,dtype=float,requirements=['C','A'])
        if originalIDNodes is not None:
            self.originalIDNodes = np.require(originalIDNodes,dtype=int,requirements=['C','A'])
        elif generateOriginalIDs:
            self.originalIDNodes = np.arange(self.GetNumberOfNodes())

    def GetDimensionality(self):
        """
        return the dimensionality 2/3
        """
        return self.nodes.shape[1]

    def GetNumberOfElements(self, dim = None):
        """
        Compute and return the total number of elements in the mesh
        """
        numberOfElements = 0
        if dim == None:
            for elemname, data in self.elements.items():
                numberOfElements += data.GetNumberOfElements()
        else:
            for elemname, data in self.elements.items():
                if ElementNames.dimension[elemname] == dim:
                    numberOfElements += data.GetNumberOfElements()
        return numberOfElements

    def MergeElements(self,other,force=False):
        """
        Merge the element for a second mesh into this
        the nodes array must be the same (not only equal)

        the user can force the merge if needed (force variable)
        """
        if (self.nodes is not other.nodes) and (not force) :
            raise(RuntimeError("the two meshes does not share the same nodes fields (potentially dangerous)"))

        for name,data in other.elements.items():
            self.GetElementsOfType(name).Merge(data)

    def ComputeGlobalOffset(self):
        """
        Recompute the Global Offset,
        This is necessary for some operation.
        Recomendation : Call it after changing the topology
        """

        cpt = 0
        for data in self.elements.values():
            data.globaloffset = cpt
            cpt += data.GetNumberOfElements()

    def ComputeBoundingBox(self):
        """
        to recumpute the bounding box
        """
        self.boundingMin = np.amin(self.nodes, axis=0)
        self.boundingMax = np.amax(self.nodes, axis=0)

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
        for data in self.elements.values():
            w = np.where(data.originalIds[:data.cpt] == oid)
            if len(w[0]) > 0 :
                data.tags.CreateTag(tagname,False).AddToTag(w[0])
                break
        else:
            raise Exception("No element with id " + str(oid)) #pragma: no cover

    def AddElementsToTag(self,globalElemNumbers,tagname):
        """
        add elements (using the global element number) to a tag (tagname)
        # you must compute the globaloffset first to make this function work
        """
        elementNotTreated = np.unique(globalElemNumbers)
        cpt = 0
        for data in self.elements.values():
            cpt2 = data.GetNumberOfElements() +cpt
            f = elementNotTreated < cpt2
            dataToTreat = elementNotTreated[f]
            if len(dataToTreat):
                data.tags.CreateTag(tagname,False).AddToTag(dataToTreat-data.globaloffset)
            elementNotTreated = elementNotTreated[np.logical_not(f)]
            cpt += data.GetNumberOfElements()

        if len(elementNotTreated):
            raise Exception("No element found") #pragma: no cover


    def AddElementToTag(self,globalElemNumber,tagname):
        """
        add a element (using the global element number) to a tag (tagname)
        # you must compute the globaloffset first to make this function work
        """
        for data in self.elements.values():
            if data.AddElementToTag(globalElemNumber,tagname):
                return
        raise Exception("No element found") #pragma: no cover

    def DeleteElemTags(self, tagNames):
        """
        delete element tags
        """
        #check not a string but a list like
        assert not isinstance(tagNames, str)
        for data in self.elements.values():
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

    def GetElementsOriginalIDs(self,dim = None):
        """
        return a single list with all the originalid concatenated
        """
        res = np.empty(self.GetNumberOfElements(dim=dim),dtype=np.int)
        cpt = 0
        from BasicTools.Containers.Filters import ElementFilter
        for name,data,ids in ElementFilter(self,dimensionality = dim):
            res[0+cpt:len(ids)+cpt] = data.originalIds[ids]+data.originalOffset
            cpt += len(ids)
        return res

    def SetElementsOriginalIDs(self,originalIDs):
        """
        Set from a single list all the originalid
        """
        cpt = 0
        for data in self.elements.values():
            data.originalIds = originalIDs[cpt:data.GetNumberOfElements()+cpt]
            cpt += data.GetNumberOfElements()


    def GetElementsInTag(self,tagname,useOriginalId=False) :
        """
        return a list with the ids of the elements in a tag
        The user must compute the globaloffset first to make this function work
        """
        res = np.zeros((self.GetNumberOfElements(),),dtype=np.int)
        cpt =0
        for data in self.elements.values():
            if tagname in data.tags:
                tag = data.tags[tagname].GetIds()
                if useOriginalId:
                    res[cpt:cpt+len(tag) ] = data.originalIds[tag]
                else:
                    res[cpt:cpt+len(tag) ] = data.globaloffset+tag
                cpt +=  len(tag)
        return res[0:cpt]

    def PrepareForOutput(self):
        """
        function to free the extra memory used during a incremental creation of a mesh
        and final treatement (offset computation)
        """
        self.nodesTags.Tighten()
        for data in self.elements.values():
            data.tighten()
        self.ComputeGlobalOffset()
        self.VerifyIntegrity()

    def __str__(self):
        res  = "UnstructuredMesh \n"
        res += "  Number Of Nodes    : {} \n".format(self.GetNumberOfNodes())
        res += "    Tags : " + " ".join( ["("+x.name+":"+str(len(x))+")" for x in  self.nodesTags ]) + "\n"

        res += "  Number Of Elements : {} \n".format(self.GetNumberOfElements())
        for data in self.elements.values():
            res += str(data)
        if len(self.nodeFields.keys()) > 0:
            res += "  nodeFields         : " + str(list(self.nodeFields.keys())) + "\n"
        if len(self.elemFields.keys()) > 0:
            res += "  elemFields         : " + str(list(self.elemFields.keys())) + "\n"
        return res

    def VerifyIntegrity(self):
        #verification nodes an originalIdNodes are compatible
        if len(self.nodes.shape) !=2:
            raise Exception("Error in the shape of nodes")

        if self.nodes.flags['C_CONTIGUOUS'] is False:
            raise Exception("Error in the order of nodes")

        if len(self.originalIDNodes.shape) !=1:
            print(self.originalIDNodes.shape)
            raise Exception("Error in the shape of originalIDNodes")

        if self.originalIDNodes.flags['C_CONTIGUOUS'] is False:
            raise Exception("Error in the order of originalIDNodes")

        if self.originalIDNodes.shape[0] != self.nodes.shape[0]:
            print(self.originalIDNodes.shape)
            print(self.nodes.shape)
            raise Exception("nodes and originalIDNodes are incompatible")

        nbnodes = self.nodes.shape[0]

        #verification of min max in nodes tags
        for elemtype,data in self.nodesTags.items():
            ids = data.GetIds()
            if len(ids) == 0:
                continue
            if ids[0] < 0:
                raise Exception("Ids of '"+ elemtype +"' tag out of bound (<0)")
            if ids[-1] >= nbnodes:
                print(ids)
                print(nbnodes)
                raise Exception("Ids of '"+ elemtype +"' tag out of bound > nbnodes")

        # verification of elements
        for elemtype, data in self.elements.items():
            if data.connectivity.flags['C_CONTIGUOUS'] is False:
                raise Exception("Error :  connectivity not C_CONTIGUOUS")

            if len(data.connectivity.shape) != 2:
                raise Exception("Wrong shape of connetivity of elements '"+ elemtype +"'")

            if data.connectivity.shape is False:
                raise Exception("Error in the order of connecitivty")

            if data.originalIds.shape[0] != data.connectivity.shape[0]:
                print(data.originalIds.shape[0])
                print(data.connectivity.shape[0])
                raise Exception("connectivity and originalIds are incompatible '"+elemtype+"'")

            if data.originalIds.flags['C_CONTIGUOUS'] is False:
                raise Exception("Error in the order of originalIds")

            if ElementNames.numberOfNodes[elemtype] != data.connectivity.shape[1]:
                print(elemtype)
                print(ElementNames.numberOfNodes[elemtype])
                print(data.connectivity.shape[1])
                raise Exception("Incompatible number of columns of the connectivity array")

            if len(data.connectivity) and  np.amin(data.connectivity) < 0:
                raise Exception("connectivity of '"+elemtype+"' out of bound (<0)")

            if len(data.connectivity) and np.amax(data.connectivity) >= nbnodes:
                print(data.connectivity)
                print(nbnodes)
                raise Exception("connecitivty of '"+elemtype+"' out of bound > nbnodes")

            nbelements = data.connectivity.shape[0]


            #verification of min max in nodes tags
            for tagName,tagData in data.tags.items():
                ids = tagData.GetIds()
                if len(ids) == 0:
                    continue
                if ids[0] < 0:
                    print(nbelements)
                    print(ids)
                    raise Exception("Ids of '"+ tagName +"' tag out of bound (<0)")
                if ids[-1] >= nbelements:
                    print(nbelements)
                    print(ids)
                    raise Exception("Ids of '"+ tagName +"' tag out of bound > nbelements")

def CheckIntegrity():
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,2,3],[0,1,0]], [[0,1,2]])
    res.ConvertDataForNativeTreatment()

    elements = res.GetElementsOfType(ElementNames.Triangle_3)

    elements = res.GetElementsOfType(ElementNames.Bar_2)
    elements.AddNewElement([1,2],1)

    elements.GetNumberOfNodesPerElement()

    print(res.IsUnstructured())

    res.ComputeGlobalOffset()
    res.AddElementToTag(1,"SecondElement")

    if res.GetNumberOfElements() != 2:
        raise Exception()
    res.ComputeGlobalOffset()

    print(res.GetDimensionality())
    res.ComputeBoundingBox()
    print(res.boundingMin)
    print(res.boundingMax)

    res.nodesTags.CreateTag('toto')
    print(res.GetNodalTag('toto'))
    print(res.GetNodalTag('toto2'))

    res.AddElementToTagUsingOriginalId(1,"bars")

    if res.GetPosOfNodes()[1,1] != 2:
        raise Exception()

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
        raise #pragma: no cover
    except:
        pass

    resII.MergeElements(res,force=True)
    print("----")
    print(resII)
    print(resII.GetNumberOfElements(dim=2))
    print(resII.elements.GetTagsNames())
    del resII.elements[ElementNames.Triangle_3]
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
