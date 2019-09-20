# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject  as BOO
import BasicTools.Containers.ElementNames as EN
"""
Filter classes of elements

"""
class Filter(BOO):
    """
    Base class to construct node and element filters

    """
    def __init__(self,mesh = None, zones = None, tags = None, zone = None, tag = None):
        """
        Constructor

        :param UnstructuredMesh mesh: the mesh to filter
        :param list(Callable) zones: a list zone to treat
        :param Callable zone: a zone to treat
        :param List(str) tags: the list of tag names to treat
        :param str tag: a tag name to add to the list

        """
        super(Filter,self).__init__()

        self.tags = list()
        if tags is not None:
            self.SetTags(tags)
        if tag is not None:
            self.AddTag(tag)

        self.zones = list()
        if zones is not None:
            self.SetZones(zones)
        if zone is not None:
            self.AddZone(zone)

        self.mesh = mesh

    def SetTags(self,tagNames):
        """
        Set the tag list name to treate

        :param list(str) tagNames: list of tag names

        """
        self.tags = list(tagNames)

    def AddTag(self,tagName):
        """
        Add a tagname to the list of tag to treat

        :param str tagName: tag name to add
        """
        self.tags.append(tagName)

    def SetZones(self,zonesList):
        """
        Set the zone list to treate

        :param list(callable) zonesList: list of zone

        """
        self.zones = zonesList

    def AddZone(self,zone):
        """
        Add a zone to the list of zone to be treated by the filter

        :param callable zone: a callable object it takes one argument with the
        points positions
        """
        self.zones.append(zone)

    def _CheckTags_(self,tags,numberOfObjects):
        """
        Internal fuction to compute the ids to be treated based on the tags

        :param Tags tags: The tags container
        :param int numberOfObjects: the total number of object (number of points,
           or number of element in the current element container)

        :return list(int) or None:

        """
        if len(self.tags) == 0:
            return None

        res = np.zeros(0,dtype=int)
        for tag in self.tags:
            if tag in tags:
                res = np.union1d(res,tags[tag].GetIds())
        return res

    def intersect1D(self,first,second):
        """
        Function to generate an intersection of two vectors (like np.intersect1d)
        but with the particularity of treat the case where the inputs can be None


        :param list(int)/np.array(int)/None first:  first vector of indices
        :param list(int)/np.array(int)/None second:  second vector of indices
        :return: The intersection of two list
        :rtype: Union[list(int),None]

        """
        if first is None:
            if second is None:
                return None
            else:
                return second
        else:
            if second is None:
                return first
            else:
                return np.intersect1d(first,second,assume_unique=True)

class NodeFilter(Filter):
    """

       Specialized class for node filtering zone and tag

    """
    def __init__(self,mesh, etags = None, etag = None, **kwargs):
        """
        Constructor: all the parammeter are passed to the base class

        """
        self.etags = list()
        if etags is not None:
            self.SetETags(etags)
        if etag is not None:
            self.AddETag(etag)

        super(NodeFilter,self).__init__(mesh,**kwargs)

    def SetETags(self,tagNames):
        """
        Set the tag list name to treate

        :param list(str) tagNames: list of tag names

        """
        self.etags = list(tagNames)

    def AddETag(self,tagName):
        """
        Add a tagname to the list of tag to treat

        :param str tagName: tag name to add
        """
        self.etags.append(tagName)

    def _CheckZones_(self,pos,numberOfObjects):
        """
        Internal function to compute the ids to be treated based on the zones

        :param np.array((n,3)) pos: np.array with the positions to be treated
        :param int numberOfObjects: total number of points
        :return list(int) or None: list of nodes to treat or None for all nodes
        """
        if len(self.zones) == 0:
            return None

        res = np.zeros(numberOfObjects,dtype=bool)

        for zone in self.zones:
            np.logical_or(res, zone(self.mesh.nodes)<=0 ,out=res)

        return np.where(res)[0]

    def GetIdsToTreat(self):
        """
        Main function of this class

        :return list(int): the filtered ids

        """
        if len(self.etags):
            ff =ElementFilter(self.mesh,tags=self.etags)
            class OP():
                 def __init__(self):
                     self.set = set()

                 def __call__(self,name,data,ids):
                     self.set.update(data.GetNodesIdFor(ids))
            op = OP()
            ff.ApplyOnElements(op)
            if len(op.set):
               resE = list(op.set)
            else:
                resE = None
        else:
            resE = None

        res  = self._CheckTags_(self.mesh.nodesTags,self.mesh.GetNumberOfNodes())
        res2 = self._CheckZones_(self.mesh.nodes,self.mesh.GetNumberOfNodes() )
        res3 = self.intersect1D(res,res2)
        res3 = self.intersect1D(res3,resE)
        if res3 is None:
            return range(self.mesh.GetNumberOfNodes())
        else:
            return res3

    def ApplyOnNodes(self,op):
        """
        Function to apply filter using an operator

        :param callable op: An instance of a callable object, the object can have
        the PreCondition function and/or the Postcondition function. Theses
        functions are called (if exist) (with the mesh as the first argument)
        before and after the main call ( op(mesh,nodes,ids) )


        """

        pc = getattr(op,"PreCondition",None)

        if callable(pc):
            pc(self.mesh)

        op(self.mesh,self.mesh.nodes,self.GetIdsToTreat() )

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)


class ElementFilter(Filter):
    """
       Specialized class for element filtering by dimensionality, zone and tag

       for the zones three types of treatment are possible:
           Based on the center of the element: self.zoneTreatment = "center"
           if all nodes of the element are in the zone ; self.zoneTreatment = "allnodes"
           if at least one nodes of the element is in the zone ; self.zoneTreatment = "leastonenode"

    """

    def __init__(self,mesh,dimensionality=None,**kwargs):
        """
        constructor : the dimensionality parameter can be used to select a type
        of element based in the dimension (tris are 2D, tetras are 3D...). All
        the rest of the argument are passed to the base class

        :param UnstructuredMesh mesh: the mesh to filter
        :param int dimensionality: the dimensionality filter, [-3 -2 -1 0 1 2 3 or None]
            the - sign is for the complementary part (-2 = all non 2D elements)


        """
        super(ElementFilter,self).__init__(mesh,**kwargs)
        self.dimensionality = dimensionality

        self.zoneTreatment = "center" # "center", "allnodes", "leastonenode"

    def SetDimensionality(self,dim):
        """
        Set the dimensionality of the elements to be treated

        :param int or None dimensionality: the dimensionality filter, [-3 -2 -1 0 1 2 3 or None]
            the - sign is for the complementary part (-2 = all non 2D elements).
            Set to None to reset (not to apply dimensionality as a criteria)
        """
        self.dimensionality = dim

    def _CheckDimensionality_(self,elements):
        """
        Internal function check if a type of elemtn must be treated based on
        the dimensionality criteria

        :param ElementsContainer elements: Elements to treat
        :return: True, if the elements must be treated
                              False, otherwise
                              None, dimensionality criteria not active
        :rtype  bool/None:
        """

        if self.dimensionality is None:
            return None
        else:
            eldim = EN.dimension[elements.elementType]
            if self.dimensionality  >= 0:
                if eldim != self.dimensionality:
                    return False
            else:
                if eldim == -self.dimensionality:
                    return False
        return True


    def _CheckZones_(self, elements):
        """
        Internal function to compute the ids of the elements to be treated based
        on the zones.

        :param ElementsContainer elements: Elements to treat
        :return list(int) or None: list of element to treat
                                   None, if all the elements must be treated

        """
        if len(self.zones) == 0:
            return None

        if self.zoneTreatment == "center":
            from BasicTools.Containers.MeshTools import GetElementsCenters
            centers = GetElementsCenters(nodes=self.mesh.nodes,elements=elements)

        numberOfObjects = elements.GetNumberOfElements()

        res = np.zeros(numberOfObjects,dtype=bool)

        for zone in self.zones :
            if self.zoneTreatment == "center":
                res2 = zone(centers)<=0
            elif self.zoneTreatment == "allnodes":
                z = zone(self.mesh.nodes)
                res2 = np.sum(z[elements.connectivity],axis=1) == elements.GetNumberOfNodesPerElement()
            elif self.zoneTreatment == "leastonenode":
                z = zone(self.mesh.nodes)
                res2 = np.sum(z[elements.connectivity],axis=1) > 0
            else:#pragma: no cover
                raise(Exception("zoneTreatment unknown"))

            np.logical_or(res, res2 ,out=res)

        return np.where(res)[0]


    def GetIdsToTreat(self,elements):
        """
        Main function of this class.

        :param ElementsContainer elements: Elements to treat
        :return: the filtered ids for the elements
        :rtype: list(int)
        """

        if self._CheckDimensionality_(elements) == False:
            return []

        res  = self._CheckTags_(elements.tags,elements.GetNumberOfElements())
        res2 = self._CheckZones_(elements)
        res3 = self.intersect1D(res,res2)
        if res3 is None:
            return range(elements.GetNumberOfElements())
        else:
            return res3

    def __iter__(self):
        """
        Iteration interface to ease the use of the filter

        :example:

            myFilter = ElementFilter(myMesh,dimensionality=2)

            for name,elements,ids in myFilter:

                print("This function is called only for 2D elements")

                print("Number of 2D element of type " + str(name)+ " is : "  + str(len(ids) )

        """
        for name,data in self.mesh.elements.items():
            ids = self.GetIdsToTreat(data)
            if len(ids) == 0: continue
            yield name, data, ids

    def ApplyOnElements(self,op):
        """
        Function to apply the filter  using an operator

        :param callable op: An instance of a callable object, the object can have
            the PreCondition function and/or the Postcondition function. Theses
            functions are called (if exist) (with the mesh as the first argument)
            before and after the main call ( op(name,elements,ids) )
        """
        pc = getattr(op,"PreCondition",None)

        if callable(pc):
            pc(self.mesh)

        for name,elements,ids in self:
            op(name,elements,ids)

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)


def CheckIntegrity( GUI=False):
    """
    .. literalinclude:: ../../Containers/Filters.py
       :pyobject: CheckIntegrity

    """
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    nx = 11; ny = 12; nz = 13;
    mesh = CreateCube(dimensions=[nx,ny,nz],origin=[0,0,0.], spacing=[1./(nx-1),1./(ny-1), 10./(nz-1)], ofTetras=True )
    print(mesh)



    class NOP():
        def __init__(self):
            self.cpt = 0

        def PreCondition(self,mesh):
            self.cpt = 0

        def __call__(self,mesh,node,ids):
            self.cpt += len(ids)

        def PostCondition(self,mesh):
            print("The counter is at {}".format(self.cpt) )

    ff = NodeFilter(mesh)
    ff.ApplyOnNodes(NOP())


    ff = NodeFilter(mesh,tags=["x0y0z0","x0y0z1",],tag="x1y0z0")
    ff.AddTag("x1y0z1")


    op = NOP()
    print(ff)
    ff.ApplyOnNodes(op)

    if op.cpt != 4:#pragma: no cover
        raise(Exception("Error finding the point"))

    ff.AddZone(lambda p: -1)

    op = NOP()
    print(ff)
    ff.ApplyOnNodes(op)
    if op.cpt != 4:#pragma: no cover
        raise(Exception("Error finding the point"))

    # example of counting the number of element in the eTag X0
    cpt = 0
    ff = ElementFilter(mesh)
    ff.AddTag("X0")
    for name,data,ids in ff:
        cpt += len(ids)

    print("Number of Element in tag X0 {}".format(cpt))

    ff = ElementFilter(mesh, zone = lambda p: (p[:,2]-mesh.boundingMin[2]-0.001),zones=[])
    for name,data,ids in ff:
        data.tags.CreateTag("ZZ0").SetIds(ids)


    # example of counting the number of element in the eTag ZZ0

    class OP():
        def __init__(self):
            self.cpt = -1

        def PreCondition(self,mesh):
            self.cpt = 0

        def __call__(self,name,data,ids):
            self.cpt += len(ids)
            print(name)

        def PostCondition(self,mesh):
            print("The counter is at {}".format(self.cpt) )


    ff.zoneTreatment = "allnodes"
    ff.ApplyOnElements(op)
    ff.zoneTreatment = "leastonenode"
    ff.ApplyOnElements(op)


    cpt = 0
    ff = ElementFilter(mesh)
    ff.ApplyOnElements(op)



    ff.AddTag("ZZ0")
    op = OP()
    ff.SetDimensionality(-2)
    ff.ApplyOnElements(op)

    ff.SetDimensionality(2)

    ff.ApplyOnElements(op)

    if op.cpt != (2*(nx-1)*(ny-1)) : # pragma: no cover
        raise(Exception("Error in the number of elements in the tag = " + str(op.cpt)+ " must be " + str((2*(nx-1)*(ny-1)))))

    print("Number of Element in tag ZZ0 {}".format(op.cpt))
    print(mesh)

    if GUI :# pragma: no cover
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

    return "ok"
if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))# pragma: no cover
