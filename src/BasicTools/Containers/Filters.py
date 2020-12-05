# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject  as BOO
import BasicTools.Containers.ElementNames as EN
"""
Filter classes of elements

"""
class Filter(BOO):
    """
    Base class to construct node and element filters

    :param UnstructuredMesh mesh: the mesh to filter
    :param list(Callable) zones: a list zone to treat
    :param Callable zone: a zone to treat
    :param List(str) tags: the list of tag names to treat
    :param str tag: a tag name to add to the list
    """
    def __init__(self,mesh = None, zones = None, tags = None, zone = None, tag = None):
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

        :rtype: list(int) or None
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
    def __init__(self,mesh=None, etags = None, etag = None, **kwargs):
        """
        Constructor: all the parammeter are passed to the base class
        """
        self.etags = list()
        if etags is not None:
            self.SetETags(etags)
        if etag is not None:
            self.AddETag(etag)

        super(NodeFilter,self).__init__(mesh=mesh,**kwargs)

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

        :return: the filtered ids
        :rtype: list(int)
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

class FilterOP(BOO):
    """
      Specialized class to compute the operation over filters
    """
    def __init__(self,mesh=None,filters=None):
        super(FilterOP,self).__init__()

        if filters is not None:
            self.filters = filters
        else:
            self.filters = []

        self.mesh = mesh

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, m):
        self._mesh = m
        for f in self.filters:
            f.mesh = m

    def Complementary(self):
        for name,data in self.mesh.elements.items():
            ids = self.GetIdsToTreat(data)
            if len(ids) == data.GetNumberOfElements(): continue
            mask = np.ones(data.GetNumberOfElements(),dtype=bool)
            mask[ids] = False
            yield name, data, np.where(mask)[0]

    def __iter__(self):
        """
        Iteration interface to ease the use of the filter

        :example:

            myFilter = DerivedElementFilterName(myMesh)
            myFilter.filters.append(myOtherFilter1)
            myFilter.filters.append(myOtherFilter2)

            for name,elements,ids in myFilter:

                print("This function is called on the union of the 2 filters")

                print("Number of element of type " + str(name)+ " is : "  + str(len(ids) )
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

class UnionElementFilter(FilterOP):
    """
      Specialized class to compute the union of filter (add)
    """
    def __init__(self,mesh=None,filters=None):
        super(UnionElementFilter,self).__init__(mesh=mesh,filters=filters)

    def GetIdsToTreat(self, data):
        ids = set()
        for ff in self.filters:
            ids.update(ff.GetIdsToTreat(data))
        return list(ids)

class IntersectionElementFilter(FilterOP):
    """
      Specialized class to compute the intersection of filters
    """
    def __init__(self,mesh=None,filters=None):
        super(IntersectionElementFilter,self).__init__(mesh=mesh,filters=filters)

    def GetIdsToTreat(self, data):
        ids = None
        for ff in self.filters:
            if ids is None:
                ids = ff.GetIdsToTreat(data)
            else:
                ids = np.intersect1d(ids,ff.GetIdsToTreat(data) )
        return list(ids)


class ElementFilter(Filter):
    """
       Specialized class for element filtering by dimensionality, zone and tag

       for the zones three types of treatment are possible:
           Based on the center of the element: self.zoneTreatment = "center"
           if all nodes of the element are in the zone ; self.zoneTreatment = "allnodes"
           if at least one nodes of the element is in the zone ; self.zoneTreatment = "leastonenode"

    """

    def __init__(self,mesh=None,dimensionality=None,**kwargs):
        """
        constructor : the dimensionality parameter can be used to select a type
        of element based in the dimension (tris are 2D, tetras are 3D...). All
        the rest of the argument are passed to the base class

        :param UnstructuredMesh mesh: the mesh to filter
        :param int dimensionality: the dimensionality filter, [-3 -2 -1 0 1 2 3 or None]
            the - sign is for the complementary part (-2 = all non 2D elements)


        """
        super(ElementFilter,self).__init__(mesh=mesh,**kwargs)
        self.dimensionality = dimensionality

        self.zoneTreatment = "center" # "center", "allnodes", "leastonenode"

    def __str__(self):
        res = "ElementFilter\n"
        res += "  dimensionality: "+ str(self.dimensionality) + " \n"
        res += "  tags          : "+ str(self.tags) + " \n"
        res += "  zones         : "+ str(self.zones) + " \n"
        res += "  zoneTreatment : "+ str(self.zoneTreatment) + " \n"
        return res

    def SetZoneTreatment(self,zt):
        if zt in ["center", "allnodes", "leastonenode"]:
            self.zoneTreatment = zt
        else:
            raise(Exception(f"Zone treatment not valide ({zt}), possible options are : center, allnodes, leastonenode"))

    def SetDimensionality(self,dim):
        """
        Set the dimensionality of the elements to be treated

        :param dimensionality: the dimensionality filter, [-3 -2 -1 0 1 2 3 or None]
            the - sign is for the complementary part (-2 = all non 2D elements).
            Set to None to reset (not to apply dimensionality as a criteria)
        :type dimensionality: int or None
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
        :rtype: Bool/None
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
        :return list of element to treat (None, if all the elements must be treated)
        :rtype: linst(int) or None
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
                z = zone(self.mesh.nodes)<=0
                res2 = np.sum(z[elements.connectivity],axis=1) == elements.GetNumberOfNodesPerElement()
            elif self.zoneTreatment == "leastonenode":
                z = zone(self.mesh.nodes)<=0
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
        self.zonesField = res2
        res3 = self.intersect1D(res,res2)
        if res3 is None:
            return range(elements.GetNumberOfElements())
        else:
            return res3

    def Complementary(self):
        for name,data  in self.mesh.elements.items():
            ids = self.GetIdsToTreat(data)
            if len(ids) == data.GetNumberOfElements(): continue
            mask = np.ones(data.GetNumberOfElements(),dtype=bool)
            mask[ids] = False
            yield name, data, np.where(mask)[0]

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

def ElementFilterToImplicitField(ff, pseudoDistance=2):
    """ Function to generate a iso zero levelset on the mesh to represent
    the shape of the filter. This discretise iso zero on the mesh cant always
    hold a 'perfect' representation of the filter, so a modified iso zero is
    created. An aditional parameter pseudoDistance can be encreased to create
    a pseudo distance. This field is created using the connectivity and not
    the real distances.

    arguments:
        ff : a ElementFilter
        pseudoDistance : distance (in number of elements) to compute pseudo
           distance

    return
        iso zeros of the ElementFilter Geometry (numpy array)
    """
    def UpdateInsideOutsideNodes(mesh, phi, insideNodes=None,outsideNodes=None, iso=0.0):
        """ function to build masks (insideNodes, outsideNodes) with the information
        about if a particular nodes is connected (	through an element) to the
        inside phi < iso or to the ouside phi > iso
        """
        for name, data  in mesh.elements.items():
            if mesh.GetDimensionality() == EN.dimension[name]:
                phis = phi[data.connectivity]
                if insideNodes is not None:
                    elmask = np.any(phis<iso,axis=1)
                    insideNodes[ data.connectivity[elmask] ] = True
                if outsideNodes is not None:
                    elmask = np.any(phis>iso,axis=1)
                    outsideNodes[ data.connectivity[elmask] ] = True

    mesh = ff.mesh
    phi = np.zeros(mesh.GetNumberOfNodes())
    insideNodes = np.zeros(mesh.GetNumberOfNodes(),dtype=bool)
    dim = 0
    for name, data, ids  in ff:
        insideNodes[data.connectivity[ids,:]]  = True
        dim = max(dim,EN.dimension[name])

    outsideNodes = np.zeros(mesh.GetNumberOfNodes(),dtype=bool)
    for name, data, ids  in ff.Complementary():
        outsideNodes[data.connectivity[ids,:]]  = True

    phi[insideNodes] = -1
    phi[outsideNodes] = 1
    phi[np.logical_and(insideNodes,outsideNodes)] = 0

    And = np.logical_and

    if dim == mesh.GetDimensionality():
        # correction of point values
        # if a point have only zeros or negative values on neighbors point then the point is set to inside
        # if a point have only zeros or positive values on neighbors point then the point is set to outside

        insideNodes.fill(False)
        outsideNodes.fill(False)

        UpdateInsideOutsideNodes(mesh,phi,insideNodes,outsideNodes, 0)
        mask = phi == 0
        phi[And(mask, And(np.logical_not(outsideNodes),insideNodes ))] = -1/2
        phi[And(mask, And(np.logical_not(insideNodes), outsideNodes))] = 1/2


    insideWork = True
    outsideWork = True
    for i in range(1, pseudoDistance):
        if insideWork:
            insideNodes.fill(False)
            UpdateInsideOutsideNodes(mesh, phi, insideNodes, None, float(i) )
            mask = phi == i
            finalmaks = And(mask, np.logical_not(insideNodes) )
            if np.any(finalmaks):
                phi[finalmaks] = i+1
            else:
                insideWork = False

        if outsideWork:
            outsideNodes.fill(False)
            UpdateInsideOutsideNodes(mesh, phi, None, outsideNodes, float(-i) )
            mask = phi == -i
            finalmaks = And(mask, np.logical_not(outsideNodes) )
            if np.any(finalmaks):
                phi[finalmaks] = -(i+1)
            else:
                outsideWork = False

        if not outsideWork and not insideWork:
            break
    return phi

def CheckIntegrity( GUI=False):
    """
    .. literalinclude:: ../../src/BasicTools/Containers/Filters.py
       :pyobject: CheckIntegrity
    """
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
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

    ## to check if a filter can be used 2 times
    for name,data,ids in ff:
        data.tags.CreateTag("ZZ0O").SetIds(ids)


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


    op = OP()

    ff.SetZoneTreatment("allnodes")
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

    if GUI: # pragma: no cover
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

    ff = ElementFilter(mesh,tag="Some")
    mesh.elements[EN.Tetrahedron_4].tags.CreateTag("Some").SetIds(list(range(4000)))


    from BasicTools.Helpers.Timer import Timer
    a = Timer("ElementFilterToImplicitField").Start()
    phi = ElementFilterToImplicitField(ff, pseudoDistance=50)
    a.Stop()
    print(a)

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    WriteMeshToXdmf("test.xdmf",mesh, PointFields=[phi],PointFieldsNames=["Phi"] )

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=False)) # pragma: no cover
