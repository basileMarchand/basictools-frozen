# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.NumpyDefs import PBasicIndexType
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers import Filters

class DofNumberingNumpy(BaseOutputObject ):
    def __init__(self):
        super(DofNumberingNumpy,self).__init__()
        self.numbering = dict()
        self.size = 0
        self._doftopointLeft = None
        self._doftopointRight= None
        self._doftocellLeft =  None
        self._doftocellRight = None
        self._almanac = None

        self.fromConnectivity = True

    def __getitem__(self,key):
        if key == "size":
           # print("Please use the new API of DofNumbering : DofNumbering.size")
            return self.size
        if key == "fromConnectivity":
           # print("Please use the new API of DofNumbering : DofNumbering.fromConnectivity")
            return self.fromConnectivity
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
        self._doftopointLeft  = range(self.size)
        self._doftopointRight = range(self.size)
        self._doftocellLeft   = []
        self._doftocellRight  = []
        self.fromConnectivity = True

        for name,data in mesh.elements.items():
            self.numbering[name] = data.connectivity

        self._almanac = dict()

        data = np.arange(self.size)
        self._almanac[('P','P')] = (data,data,data)
        self.totalNumberOfPoints = mesh.GetNumberOfNodes()
        self.computeDofToPoint()
        return self

    def GetKeyFromNameAndIdxI(self,on,name,idxI):
        if on == 'C':
            on = ('C',name)
        elif on == 'F2':
            on = ('C',EN.faces2[name][idxI][0])
        elif on == 'F':
            on  = ('C',EN.faces[name][idxI][0])
        elif on == "IP":
            on =  ("IP",name)
        elif on == "G":
            on = ('G','G')
        elif on == 'P':
            on = ('P','P')

        else:
            print(on)
            raise
        return on

    def ComputeNumberingGeneral(self,mesh,space,elementFilter=None,discontinuous=False):

        if discontinuous:
            raise(Exception("ERROR!! discontinuous=True, for the moment this is not available"))

        if self._almanac is not None:
            raise(Exception("Cant pass 2 time in ComputeNumberingGeneral"))

        self.fromConnectivity = False

        if elementFilter is None:
            elementFilter = Filters.ElementFilter(mesh)

        self.totalNumberOfPoints = mesh.GetNumberOfNodes()

        def GetNumberingColsSize(k):
            if k[0] == 'C':
                nn = EN.numberOfNodes[k[1]]
            elif k[0] == "IP":
                nn = EN.numberOfNodes[k[1]] +1
            elif k[0] == "G":
                nn = 1
            elif k[0] == "P":
                nn = 1
#            elif k == "G":
#                nn = 1
#            elif k == 'P':
#                nn = 1
            else:
                print(k)
                raise
                nn = EN.numberOfNodes[k]
            return nn

        sizes = dict()

        self.PrintVerbose("Numbering counting ")
        #count the total number of shape funtions per element
        for name,data, elids in elementFilter:
            sp = space[name]
            nsf= sp.GetNumberOfShapeFunctions()
            for sf in range(nsf):
                on,idxI,idxII = sp.dofAttachments[sf]
                key = self.GetKeyFromNameAndIdxI(on,name,idxI)
                sizes[key] = sizes.get(key,0) + len(elids)

        self.PrintDebug("Numbering memory allocation")
        ## allocation of matrices to store all the dofs
        storage = {}
        for k,v in sizes.items():
            nn = GetNumberingColsSize(k)
            storage[k] = np.zeros((v,nn),dtype=PBasicIndexType)-1
            sizes[k] = 0

        self.PrintDebug("Numbering generation dofs keys")
        for elemName, data , ids in elementFilter:
            self.PrintDebug("working on " + elemName)
            #BOO.ResetStartTime()
            sp = space[elemName]
            nsf= sp.GetNumberOfShapeFunctions()
            elidsConnectivity  = data.connectivity[ids,:]
            for sf in range(nsf):
                nn = GetNumberingColsSize(k)
                key, idxI, idxII = self.GetHashFor(data,sp,ids,sf,False,elidsConnectivity=elidsConnectivity)
                cpt = sizes[key]
                storage[key][cpt:cpt+len(ids),:] = idxI
                sizes[key] = cpt + len(ids)

        cpt = self.size
        tempAlmanac  = dict()

        self.PrintDebug("Numbering generation of uniques")
        # recover the unique dofs and generate the numbering
        for k,v in storage.items():
            unique, indices, inverse = np.unique(np.sort(v,axis=1),return_index=True,return_inverse=True,axis=0)
            newdofs = np.arange(len(indices)) + cpt
            tempAlmanac[k] = (unique,newdofs,inverse)
            cpt += len(indices)
            sizes[k] = 0

        self.size = cpt

        self.PrintDebug("Numbering the bulk")

        #push the new dofs to the almanac
        maxuseddim = 0

        extractorLeftSide = np.empty(self.size,dtype=PBasicIndexType)
        extractorRightSide = np.empty(self.size,dtype=PBasicIndexType)
        extractorcpt = 0

        elementcpt = 0

        for elemName, data in mesh.elements.items():
            ids = elementFilter.GetIdsToTreat(data)
            if len(ids) == 0:
                elementcpt += data.GetNumberOfElements()
                continue
            maxuseddim = max(maxuseddim,EN.dimension[elemName] )

            sp = space[elemName]
            nsf= sp.GetNumberOfShapeFunctions()
            self.numbering[elemName] = np.zeros((data.GetNumberOfElements(),nsf),dtype=PBasicIndexType)-1
            done = False
            for sf in range(nsf):
                nn = GetNumberingColsSize(k)
                on,idxI,idxII = sp.dofAttachments[sf]
                key = self.GetKeyFromNameAndIdxI(on,elemName,idxI)
                cpt = sizes[key]
                (unique,newdofs,inverse) = tempAlmanac[key]
                dofs = newdofs[inverse][cpt:cpt+len(ids)]
                self.numbering[elemName][ids,sf] = dofs
                sizes[key] = cpt + len(ids)

                # the dofs are attached to the the elements
                #print(elemName,key, done)
                if key[1] == elemName and not done:
                    done = True
                    extractorLeftSide[extractorcpt:extractorcpt+len(ids)] = ids
                    extractorLeftSide[extractorcpt:extractorcpt+len(ids)] += elementcpt
                    extractorRightSide[extractorcpt:extractorcpt+len(ids)] = dofs
                    extractorcpt += len(ids)
#                    raise

            elementcpt += data.GetNumberOfElements()

        self._doftocellLeft = np.resize(extractorLeftSide, (extractorcpt,))
        self._doftocellRight = np.resize(extractorRightSide, (extractorcpt,))

        from BasicTools.Containers.Filters import IntersectionElementFilter, ElementFilter, ComplementaryObject,UnionElementFilter
        # we work on the elements of dim < maxuseddim not present in the original filter
        complementary = ComplementaryObject(mesh=mesh, filters = [elementFilter])
        filters = [ElementFilter(mesh=mesh,dimensionality=i) for i in range(maxuseddim) ]
        #filters.append(complementary)
        outside = IntersectionElementFilter(mesh=mesh, filters=[UnionElementFilter(mesh=mesh,filters=filters) ,complementary ]  )

        #push numbering for elements touching the already numbered elements
        # for example the 2D elements on the surface of 3D elements
        self.PrintVerbose("Numbering the complementary")
        for elemName, data, ids in outside:
            sp = space[elemName]
            nsf= sp.GetNumberOfShapeFunctions()
            if elemName not in self.numbering:
                self.numbering[elemName] = np.zeros((data.GetNumberOfElements(),nsf),dtype=PBasicIndexType)
            elidsConnectivity  = data.connectivity[ids,:]
            for sf in range(nsf):
                nn = GetNumberingColsSize(k)
                on,idxI,idxII = sp.dofAttachments[sf]
                key = self.GetKeyFromNameAndIdxI(on,elemName,idxI)
                if key not in tempAlmanac:
                    continue
                (unique,newdofs,inverse) = tempAlmanac[key]
                name, idxI, idxII = self.GetHashFor(data,sp,ids,sf,False,elidsConnectivity=elidsConnectivity)
                v = np.vstack((unique,idxI))
                uniqueII, indices, inverse = np.unique(np.sort(v,axis=1),return_index=True,return_inverse=True,axis=0)
                newnewdofs = np.hstack((newdofs,np.zeros(len(idxI),dtype=PBasicIndexType)-1 ))[inverse][len(unique):]
                self.numbering[elemName][ids,sf] = newnewdofs
        self.PrintVerbose("Numbering Done")
        self._almanac = tempAlmanac
        self.computeDofToPoint()

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
        #print(self._almanac.keys())
        self.pointDofs = np.zeros(self.totalNumberOfPoints,dtype=PBasicIndexType)-1
        if ('P','P') not in self._almanac:
            return

        unique, newdofs, inverse = self._almanac[('P','P')]

        self._doftopointLeft = unique.flatten()
        self._doftopointRight = newdofs.flatten()

        self.pointDofs[self._doftopointLeft] = self._doftopointRight

    def GetDofOfPoint(self,pointid):
        return self.pointDofs[pointid]

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
        raise
        extractorLeftSide = np.empty(self.size,dtype=PBasicIndexType)
        extractorRightSide = np.empty(self.size,dtype=PBasicIndexType)

        tmpcpt = 0
        for elem,data in self.mesh.elements.items():
            unique, newdofs, inverse = self._almanac[elem]
            extractorLeftSide[tmpcpt:tmpcpt+len(inverse)] = inverse
            extractorRightSide[tmpcpt:tmpcpt+len(inverse)] = newdofs[inverse]
            tmpcpt += len(newdofs)
        # if k[0] is the elementname then k[1] is the connecivity
        # we generate the same almanac with the number of each element
        elemDic = {}
        for name,data in self.mesh.elements.items():
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
                    extractorLeftSide[tmpcpt] = self.mesh.elements[k[0]].globaloffset + localdic[k[1]]
                    extractorRightSide[tmpcpt] = v
                    tmpcpt += 1

        self._doftocellLeft = np.resize(extractorLeftSide, (tmpcpt,))
        self._doftocellRight = np.resize(extractorRightSide, (tmpcpt,))

    def GetHashFor(self,data,sp,elids,sf,discontinuous=False,elidsConnectivity=None):
        res = []
        name = data.elementType
        #self.PrintDebug("2.1")
        #print(type(elids))
        if elidsConnectivity is None:
            elidsConnectivity  = data.connectivity[elids,:]

        #self.PrintDebug("2.2")
        on,idxI,idxII = sp.dofAttachments[sf]

        key = self.GetKeyFromNameAndIdxI(on,name,idxI)

        if on == "P":
            shapeFunctionConnectivity = elidsConnectivity [:,idxI]
            return key, shapeFunctionConnectivity[:,None], idxII

        if on == "C":
            return key, elidsConnectivity, idxI

        if on == "F":
            edge = EN.faces[name][idxI]
            return key, elidsConnectivity[:,edge[1]],0

        if on == "F2":
            edge = EN.faces2[name][idxI]
            return key, elidsConnectivity[:,edge[1]],0

        if on == "G":
            """G is for global """
            return key, np.zeros((len(elids),1),dtype=PBasicIndexType), idxII

        if on == "IP" :
            return key, np.hstack((elidsConnectivity,-idxI*np.ones(len(elids),dtype=PBasicIndexType)[:,None] ) ),idxI

        return res

def CheckIntegrity(GUI=False):
    import BasicTools.FE.DofNumbering  as DN
    return DN.CheckIntegrityUsingAlgo("NumpyBase",GUI)


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover

