# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers import Filters

class DofNumberingDict(BaseOutputObject ):
    def __init__(self):
        """Class to generate the numbering.
        The internal implementation uses Dict.
        """
        super(DofNumberingDict,self).__init__()
        self.numbering = dict()
        self.size = 0
        self._doftopointLeft = None
        self._doftopointRight= None
        self._doftocellLeft =  None
        self._doftocellRight = None
        self.almanac = dict()
        self.newAlmanac = dict()

        self.fromConnectivity = True

    def __getitem__(self,key):
        if key == "size":
           # print("Please use the new API of DofNumbering : DofNumbering.size")
            return self.size
        if key == "fromConnectivity":
           # print("Please use the new API of DofNumbering : DofNumbering.fromConnectivity")
            return self.fromConnectivity
        if key == "almanac":
           # print("Please use the new API of DofNumbering : DofNumbering.almanac")
            return self.almanac

        if key == "doftopointLeft":
           # print("Please use the new API of DofNumbering : DofNumbering.doftopointLeft")
            return self.doftopointLeft
        if key == "doftopointRight":
           # print("Please use the new API of DofNumbering : DofNumbering.doftopointRight")
            return self.doftopointRight

        if key == "doftocellLeft":
           # print("Please use the new API of DofNumbering : DofNumbering.doftopointLeft")
            return self.doftocellLeft
        if key == "doftocellRight":
           # print("Please use the new API of DofNumbering : DofNumbering.doftopointRight")
            return self.doftocellRight

        return self.numbering[key]

    def get(self,key,default=None):
        if key in self.numbering:
            return self.numbering[key]
        else:
            return default

    def __contains__(self, k):
        return k in self.numbering

    def GetSize(self):
        return self.size

    def ComputeNumberingFromConnectivity(self,mesh,space):

        self.size = mesh.GetNumberOfNodes()
        self._doftopointLeft = range(self.size)
        self._doftopointRight = range(self.size)
        self._doftocellLeft =  []
        self._doftocellRight = []
        self.fromConnectivity = True

        for name,data in mesh.elements.items():
            self.numbering[name] = data.connectivity

        almanac = {}
        for i in range(self.size):
            almanac[('P', i, None)] = i
        self.almanac = almanac

        return self

    def ComputeNumberingGeneral(self, mesh, space, elementFilter=None, discontinuous=False):
        self.fromConnectivity = False
        almanac = self.almanac

        if elementFilter is None:
            elementFilter = Filters.ElementFilter(mesh)
        else:
            elementFilter.SetMesh(mesh)

        cpt = self.size
        self.PrintDebug("bulk ")
        usedDim = 0
        for name, data, elementsIds in elementFilter:
            usedDim = max(usedDim,EN.dimension[name] )
            res = self.GetHashFor(data, space[name], elementsIds, discontinuous)

            if name in self.numbering:
                dofs = self.numbering[name]
            else:
                dofs = np.zeros((data.GetNumberOfElements(),space[name].GetNumberOfShapeFunctions()), dtype=PBasicIndexType) -1

            self.PrintDebug(name + " Done")
            for i in range(len(res)):
                localRes = res[i]
                localDofs = dofs[:,i]
                for j, elementId in enumerate(elementsIds):
                    d = almanac.setdefault(localRes[j],cpt)
                    cpt += (d == cpt)
                    localDofs[elementId] = d

            self.numbering[name] = dofs
            self.PrintDebug(name + " Done Done")
        self.PrintDebug("bulk Done")
        self.PrintDebug("complementary ")
        from BasicTools.Containers.Filters import IntersectionElementFilter, ElementFilter, ComplementaryObject
        outside = IntersectionElementFilter(mesh=mesh, filters =[ElementFilter(dimensionality=usedDim-1), ComplementaryObject( filters = [elementFilter])] )

        for name,data,elementsIds in outside:
            res = self.GetHashFor(data,space[name],elementsIds,discontinuous)
            if name in self.numbering:
                dofs = self.numbering[name]
            else:
                dofs = self.numbering.setdefault(name,np.zeros((data.GetNumberOfElements(),space[name].GetNumberOfShapeFunctions()), dtype=PBasicIndexType) -1)
            for i in range(len(res)):
                localRes = res[i]
                localDofs = dofs[:,i]
                for j,elementId in enumerate(elementsIds):
                    if localDofs[elementId] >= 0:
                        continue
                    localDofs[elementId] = almanac.get(localRes[j],-1)
            self.numbering[name] = dofs
        self.PrintDebug("complementary Done")
        self.size = cpt
        #-------------------------------------------------------------------------
        self.mesh = mesh
        # we keep a reference to the mesh because we need it to compute the
        return self

    def GetDofOfPoint(self,pid):
        return self.almanac[("P",pid,None)]

    @property
    def doftopointLeft(self):
        if self._doftopointLeft is None:
            self.computeDofToPoint()
        return self._doftopointLeft

    @property
    def doftopointRight(self):
        if self._doftopointRight is None:
            self.computeDofToPoint()
        return self._doftopointRight

    def computeDofToPoint(self):
        extractorLeftSide = np.empty(self.size,dtype=PBasicIndexType)
        extractorRightSide = np.empty(self.size,dtype=PBasicIndexType)

        cpt = 0
        # if k[0] is 'P' then k[1] is the node number
        for k,v in self.almanac.items():
            if k[0] == 'P':
                extractorLeftSide[cpt] = k[1]
                extractorRightSide[cpt] = v
                cpt += 1

        self._doftopointLeft = np.resize(extractorLeftSide, (cpt,))
        self._doftopointRight = np.resize(extractorRightSide, (cpt,))

    @property
    def doftocellLeft(self):
        if self._doftocellLeft is None:
            self.computeDofToCell()
        return self._doftocellLeft

    @property
    def doftocellRight(self):
        if self._doftocellRight is None:
            self.computeDofToCell()
        return self._doftocellRight


    def computeDofToCell(self):
        mesh = self.mesh
        extractorLeftSide = np.empty(self.size,dtype=PBasicIndexType)
        extractorRightSide = np.empty(self.size,dtype=PBasicIndexType)

        cpt = 0
        # if k[0] is the elementName then k[1] is the connectivity
        # we generate the same almanac with the number of each element
        elemDic = {}
        for name, data in mesh.elements.items():
            elemDic[name] = {}
            elemDic2 = elemDic[name]
            sortedConnectivity = np.sort(data.connectivity,axis=1)

            for i in range(data.GetNumberOfElements()):
                elemDic2[tuple(sortedConnectivity[i,:])] = i

        for k,v in self.almanac.items():
            #if not k[0] in {'P',"F","F2","G"} :
            #we need the global number of the element (not the local to the element container)
            if k[0] in elemDic.keys():
                localDict = elemDic[k[0]]
                if k[1] in localDict.keys():
                    extractorLeftSide[cpt] = mesh.elements[k[0]].globaloffset + localDict[k[1]]
                    extractorRightSide[cpt] = v
                    cpt += 1

        self._doftocellLeft = np.resize(extractorLeftSide, (cpt,))
        self._doftocellRight = np.resize(extractorRightSide, (cpt,))

    def GetHashFor(self, data, sp, elementIds, discontinuous):

        numberOfShapeFunctions = sp.GetNumberOfShapeFunctions()
        res = []
        name = data.elementType

        elementsIdsConnectivity  = data.connectivity[elementIds,:]

        for j in range(numberOfShapeFunctions):
            on,idxI,idxII = sp.dofAttachments[j]
            if on == "P":
                T = "P"
                shapeFunctionConnectivity = elementsIdsConnectivity [:,idxI]

                if discontinuous :
                    res.append( [ (T+str(elementIds),x,idxII) for i,x in zip(elementIds,shapeFunctionConnectivity)  ]  )
                else:
                    res.append( [ (T,x,idxII) for x in shapeFunctionConnectivity  ]  )
            elif on == "C":
                sortedConnectivity = np.sort(elementsIdsConnectivity,axis=1)
                T = name
                if idxII is not None:
                    raise
                res.append([(T,tuple(sortedConnectivity[i,:]),idxI) for i in range(len(elementIds)) ] )
            elif on == "F2":
                edge = EN.faces2[name][idxI]
                T = edge[0]
                nn = np.sort(elementsIdsConnectivity[:,edge[1]],axis=1)
                if discontinuous :
                    res.append( [ (T+str(elementIds),tuple(x) ,0) for i,x in zip(elementIds,nn)  ]  )
                else:
                    res.append( [ (T,tuple(x),0) for x in nn  ]  )
            elif on == "F":
                edge = EN.faces[name][idxI]
                T = edge[0]
                nn = np.sort(elementsIdsConnectivity[:,edge[1]],axis=1)
                if discontinuous :
                    res.append( [ (T,tuple(x),i) for i,x in zip(elementIds,nn)  ]  )
                else:
                    res.append( [ (T,tuple(x),0) for x in nn  ]  )
            elif on == "G":
                """G is for global """
                key = (on,0,idxII)
                res.append( [ key for x in elementIds ]  )
            elif on == "IP":
                res.append( [ (on,tuple(localCoon),idxI) for localCoon,i in zip(elementsIdsConnectivity,elementIds) ]  )
            else:
                print(on)
                raise
        return res

def CheckIntegrity(GUI=False):
    import BasicTools.FE.DofNumbering  as DN
    return DN.CheckIntegrityUsingAlgo("DictBase",GUI)

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
