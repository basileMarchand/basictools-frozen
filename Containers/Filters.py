# -*- coding: utf-8 -*-

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject  as BOO
import BasicTools.Containers.ElementNames as EN
"""
Filter classes of elements

"""
class Filter(BOO):
    def __init__(self,mesh =None, zone=None, tag=None):
        super(Filter,self).__init__()
        self.tag = tag
        self.mesh = mesh
        self.zone = zone

    def SetTag(self,data):
        self.tag = data

    def SetZone(self,data):
        self.zone = data


class NodeFilter(Filter):
    def __init__(self):
        super(NodeFilter,self).__init__()

    def CheckTag(self):

        if self.tag is None:
            return range(self.GetNumberOfNodes())

        if self.tag in self.mesh.nodesTags:
            return self.mesh.nodestags[self.tag].GetIds()

        raise (Exception("Tag Not present in mesh nodes"))

    def CheckZone(self):
        if self.zone in None:
            return range(self.mesh.GetNumberOfNodes())

        return np.where(self.zone(self.mesh.nodes)<=0)[0]

    def GetIdsToTreat(self):
        res = None
        if self.tag is not None:
            res = self.CheckTag()

        if self.zone is not None:
            res2 = self.CheckZone()
            if res is not None:
                res = np.intersect1D(res,res2)
            else:
                res = res2

        return res

    def ApplyOnNodes(self,op):


        pc = getattr(op,"PreCondition",None)
        if callable(pc):
            pc(self.mesh)

        op(self.mesh,self.mesh.nodes,self.GetIdsToTreat() )

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)


class ElementFilter(Filter):
    def __init__(self,mesh,zone = None, tag = None, dimensionality = None):
        super(ElementFilter,self).__init__(mesh, zone=zone, tag=tag)
        self.dimensionality = dimensionality # possible are [-3 -2 -1 0 1 2 3 ]
        # "the - sign is for the complementary -2 = all non 2D elements"

        self.zoneTreatement = "center" # "center", "allnodes", "somenodes"

    def SetDimensionality(self,data):
        self.dimensionality = data

    def CheckDimensionality(self,elemens):

        if self.dimensionality is None:
            return True
        else:
            if self.dimensionality  > 0:
                eldim = EN.dimension[elemens.elementType]
                if eldim != self.dimensionality:
                    return False
            else:
                if eldim == -self.dimensionality:
                    return False
        return True

    def CheckTag(self,elements):
        if self.tag is None:
            return range(elements.GetNumberOfElements())

        if self.tag in elements.tags:
            return elements.tags[self.tag].GetIds()

        return []


    def CheckZone(self, elements):
        if self.zone is None:
            return range(elements.GetNumberOfElements())

        if self.zoneTreatement == "center":
            from BasicTools.Containers.MeshTools import GetElementsCenters
            centers = GetElementsCenters(nodes=self.mesh.nodes,elements=elements)
            return np.where(self.zone(centers)<=0)[0]

        elif self.zoneTreatement == "allnodes":
            return np.sum(self.mesh.nodestags[self.tag].GetIds(),axis=1) == elements.GetNumberOfNodes()

        elif self.zoneTreatement == "somenodes":
            return np.sum(self.mesh.nodestags[self.tag].GetIds(),axis=1) > 0
        else:
            raise(Exception)

    def GetIdsToTreat(self,elements):
        if not self.CheckDimensionality(elements):
            return []

        res = None

        if self.tag is None:
            res = range(elements.GetNumberOfElements())
        else:
            res = self.CheckTag(elements)

        if self.zone is not None:
            res2 = self.CheckZone(elements)
            if res is not None:
                res = np.intersect1d(res,res2)
            else:
                res = res2

        return res

    def __iter__(self):
        for name,data in self.mesh.elements.items():
            ids = self.GetIdsToTreat(data)
            if len(ids) == 0: continue
            yield name,data, ids

    def ApplyOnElements(self,op):
        pc = getattr(op,"PreCondition",None)
        if callable(pc):
            pc(self.mesh)

        for name,elements,ids in self:
            op(name,elements,ids)

        pc = getattr(op,"PostCondition",None)
        if callable(pc):
            pc(self.mesh)


def CheckIntegrity( GUI=False):
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    nx = 11; ny = 12; nz = 13;
    mesh = CreateCube(dimensions=[nx,ny,nz],origin=[0,0,0.], spacing=[1./(nx-1),1./(ny-1), 10./(nz-1)], ofTetras=True )
    print(mesh)

    # example of counting the number of element in the eTag X0
    cpt = 0
    ff = ElementFilter(mesh)
    ff.SetTag("X0")
    for name,data,ids in ff:
        cpt += len(ids)

    print("Number of Element in tag X0 {}".format(cpt))

    ff = ElementFilter(mesh, zone = lambda p: (p[:,2]-mesh.boundingMin[2]-0.001))
    for name,data,ids in ff:
        data.tags.CreateTag("ZZ0").SetIds(ids)

    # example of counting the number of element in the eTag X0
    cpt = 0
    ff = ElementFilter(mesh)
    ff.SetTag("ZZ0")
    ff.SetDimensionality(2)

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
    ff.ApplyOnElements(op)

    if op.cpt != (2*(nx-1)*(ny-1)) :
        raise(Exception("Error in the number of elements in the tag = " + str(op.cpt) ))
    print("Number of Element in tag ZZ0 {}".format(op.cpt))
    print(mesh)

    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

    return "ok"
if __name__ == '__main__':

    print(CheckIntegrity( GUI=False))

