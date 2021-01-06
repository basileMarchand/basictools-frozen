# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers import Filters

class DofNumbering(BaseOutputObject ):
    def __init__(self):
        super(DofNumbering,self).__init__()
        self.numbering = dict()
        self.size = 0
        self._doftopointLeft = None
        self._doftopointRight= None
        self._doftocellLeft =  None
        self._doftocellRight = None
        self.almanac = dict()
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
        self._doftopointLeft = slice(self.size)
        self._doftopointRight = slice(self.size)
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

    def ComputeNumberingGeneral(self,mesh,space,elementFilter=None,discontinuous=False):
        self.fromConnectivity = False
        almanac = self.almanac

        if elementFilter is None:
            elementFilter = Filters.ElementFilter(mesh)

        cpt = self.size
        self.PrintDebug("bulk ")
        useddim = 0
        cctt = 0
        for name,data, elids in elementFilter:
            cctt  += len(elids)
            useddim = max(useddim,EN.dimension[name] )
            res = self.GetHashFor(data,space[name],elids,discontinuous)
            if name in self.numbering:
                dofs = self.numbering[name]
            else:
                dofs = np.zeros((data.GetNumberOfElements(),space[name].GetNumberOfShapeFunctions()), dtype=np.int_) -1

            self.PrintDebug(name + " Done")
            for i in range(len(res)):
                lres = res[i]
                ldofs = dofs[:,i]
                for j,elid in enumerate(elids):
                    d = almanac.setdefault(lres[j],cpt)
                    cpt += (d == cpt)
                    ldofs[elid] = d

            self.numbering[name] = dofs
            self.PrintDebug(name + " Done Done")
        self.PrintDebug("bulk Done")
        self.PrintDebug("complementary ")
        from BasicTools.Containers.Filters import IntersectionElementFilter, ElementFilter, ComplementaryObject
        outside = IntersectionElementFilter(mesh=mesh, filters =[ElementFilter(dimensionality=useddim-1), ComplementaryObject( filters = [elementFilter])] )
        cctt = 0
        for name,data,elids in outside:
            cctt += len(elids)
            res = self.GetHashFor(data,space[name],elids,discontinuous)
            if name in self.numbering:
                dofs = self.numbering[name]
            else:
                dofs = self.numbering.setdefault(name,np.zeros((data.GetNumberOfElements(),space[name].GetNumberOfShapeFunctions()), dtype=np.int_) -1)
            for i in range(len(res)):
                lres = res[i]
                ldofs = dofs[:,i]
                for j,elid in enumerate(elids):
                    if ldofs[elid] >= 0:
                        continue
                    ldofs[elid] = almanac.get(lres[j],-1)
            self.numbering[name] = dofs
        self.PrintDebug("complementary Done")
        self.size = cpt
        #-------------------------------------------------------------------------
        self.mesh = mesh
        # we keep a reference to the mesh because we need it to compute the
        return self


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
        extractorLeftSide = np.empty(self.size,dtype=np.int)
        extractorRightSide = np.empty(self.size,dtype=np.int)

        tmpcpt = 0
        # if k[0] is 'P' then k[1] is the node number
        for k,v in self.almanac.items():
            if k[0] == 'P':
                extractorLeftSide[tmpcpt] = k[1]
                extractorRightSide[tmpcpt] = v
                tmpcpt += 1

        self._doftopointLeft = np.resize(extractorLeftSide, (tmpcpt,))
        self._doftopointRight = np.resize(extractorRightSide, (tmpcpt,))

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
        extractorLeftSide = np.empty(self.size,dtype=np.int)
        extractorRightSide = np.empty(self.size,dtype=np.int)

        tmpcpt = 0
        # if k[0] is the elementname then k[1] is the connecivity
        # we generate the same almanac with the number of each element
        elemDic = {}
        for name,data in mesh.elements.items():
            elemDic[name] = {}
            elemDic2 = elemDic[name]
            sortedconnectivity = np.sort(data.connectivity,axis=1)

            for i in range(data.GetNumberOfElements()):
                elemDic2[tuple(sortedconnectivity[i,:])] = i

        for k,v in self.almanac.items():
            #if not k[0] in {'P',"F","F2","G"} :
            #we need the global number of the element (not the local to the element container)
            if k[0] in elemDic.keys():
                localdic = elemDic[k[0]]
                if k[1] in localdic.keys():
                    extractorLeftSide[tmpcpt] = mesh.elements[k[0]].globaloffset + localdic[k[1]]
                    extractorRightSide[tmpcpt] = v
                    tmpcpt += 1

        self._doftocellLeft = np.resize(extractorLeftSide, (tmpcpt,))
        self._doftocellRight = np.resize(extractorRightSide, (tmpcpt,))

    def GetHashFor(self,data,sp,elids,discontinuous):

        numberOfShapeFunctions = sp.GetNumberOfShapeFunctions()
        res = []
        name = data.elementType

        elidsConnectivity  = data.connectivity[elids,:]

        for j in range(numberOfShapeFunctions):
            on,idxI,idxII = sp.dofAttachments[j]
            if on == "P":
                T = "P"
                shapeFunctionConnectivity = elidsConnectivity [:,idxI]

                if discontinuous :
                    res.append( [ (T+str(elids),x,idxII) for i,x in zip(elids,shapeFunctionConnectivity)  ]  )
                else:
                    res.append( [ (T,x,idxII) for x in shapeFunctionConnectivity  ]  )
            elif on == "C":
                sortedConnectivity = np.sort(elidsConnectivity,axis=1)
                T = name
                if idxII is not None:
                    raise
                res.append([(T,tuple(sortedConnectivity[i,:]),idxI) for i in range(len(elids)) ] )
            elif on == "F2":
                edge = EN.faces2[name][idxI]
                T = edge[0]
                nn = np.sort(elidsConnectivity[:,edge[1]],axis=1)
                if discontinuous :
                    res.append( [ (T+str(elids),tuple(x) ,0) for i,x in zip(elids,nn)  ]  )
                else:
                    res.append( [ (T,tuple(x),0) for x in nn  ]  )
            elif on == "F":
                edge = EN.faces[name][idxI]
                T = edge[0]
                nn = np.sort(elidsConnectivity[:,edge[1]],axis=1)
                if discontinuous :
                    res.append( [ (T,tuple(x),i) for i,x in zip(elids,nn)  ]  )
                else:
                    res.append( [ (T,tuple(x),0) for x in nn  ]  )
            elif on == "G":
                """G is for global """
                key = (on,0,idxII)
                res.append( [ key for x in elids ]  )
            elif on == "IP":
                res.append( [ (on,tuple(lcoon),idxI) for lcoon,i in zip(elidsConnectivity,elids) ]  )
            else:
                print(on)
                raise
        return res

def ComputeDofNumbering(mesh,Space,dofs=None,sign=1,fromConnectivity=False,elementFilter=None,discontinuous=False):

    if sign != 1:
        raise(NotImplementedError())

    if fromConnectivity:
        res = DofNumbering()
        if dofs is not None or sign != 1 or elementFilter is not None or discontinuous:
           raise(Exception("cant take dofs, sign, discontinuous or elementFilter different from the default values"))
        res.ComputeNumberingFromConnectivity(mesh,Space)
        return res
    else:
        if dofs is None:
            res = DofNumbering()
        else:
            res = dofs

        res.ComputeNumberingGeneral(mesh=mesh, space=Space, elementFilter=elementFilter, discontinuous=discontinuous )
        return res

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    #res2 = CreateCube([100.,100.,100.],[-1.0,-1.0,-1.0],[2./46, 2./46,2./46])
    res2 = CreateCube([2.,2.,2.],[-1.0,-1.0,-1.0],[2./46, 2./46,2./46])
    #NodesPermutation(res2,np.array([0, 4, 6,2 ,1 ,5 ,7, 3]))

    print(res2)

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

    print("Numbering ")
    import time
    st = time.time()
    #numbering = ComputeDofNumbering(res2,LagrangeSpaceGeo,tag="X0",sign=-1)
    numbering = ComputeDofNumbering(res2,LagrangeSpaceGeo)
    numbering = ComputeDofNumbering(res2,LagrangeSpaceGeo,fromConnectivity=True)
    print(time.time()-st)

    print(numbering)
    for name,data in res2.elements.items():
        if name in numbering:
            print(name)
            print(data.connectivity)
            print(numbering[name])


    import BasicTools.Containers.ElementNames as EN
    print(res2.GetElementsOfType(EN.Quadrangle_4).tags["X1"].GetIds() )

    from BasicTools.FE.Spaces.FESpaces import ConstantSpaceGlobal

    t_numbering = ComputeDofNumbering(res2,ConstantSpaceGlobal)
    print(t_numbering)

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP2

    t_numbering = ComputeDofNumbering(res2,LagrangeSpaceP2)
    print(t_numbering)
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
