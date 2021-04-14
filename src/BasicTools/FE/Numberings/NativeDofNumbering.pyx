#distutils: language = c++
#cython: language_level = 3
cimport numpy as cnp
cnp.import_array()
from libcpp cimport bool

from libcpp.string cimport string
from eigency.core cimport *



cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM
cimport BasicTools.Containers.NativeFilters as cNF
cimport BasicTools.FE.Spaces.NativeSpace as cNS

import BasicTools.Containers.NativeUnstructuredMesh as NUM
import BasicTools.FE.Spaces.NativeSpace as NS

cdef extern from "FE/DofNumberingCpp.h" :
    cdef cppclass DofNumberingCpp:
        DofNumberingCpp() except +
        int GetSize()
        int GetFromConnectivity()
        void SetFromConnectivity(bool )
        void ComputeNumberingFromConnectivity(cNUM.NativeUnstructuredMesh* )
        void ComputeNumberingGeneral(cNUM.NativeUnstructuredMesh*, cNS.Space*, cNF.ElementFilterBase*)

        PlainObjectBase & GetNumberingFor(const string & elemtype)
        bool HasNumberingFor(const string & elemtype)
        string ToStr()

cdef class NativeDofNumbering:
    cdef DofNumberingCpp dn_cpp

    #def get(self,key,default=None):
    #    return self.numbering.get(key,default)

    def __getitem__(self,key):
        if key == "size":
           # print("Please use the new API of DofNumbering : DofNumbering.size")
            return self.size
        if key == "fromConnectivity":
           # print("Please use the new API of DofNumbering : DofNumbering.fromConnectivity")
            return self.fromConnectivity
        #print("NativeDofNumbering[]  "+key)
        #print("HasNumberingFor", self.dn_cpp.HasNumberingFor(key.encode()))
        #print("HasNumberingFor", self.dn_cpp.HasNumberingFor(key.encode()))
        if self.dn_cpp.HasNumberingFor(key.encode()):
            try:
                return  ndarray_view(self.dn_cpp.GetNumberingFor(key.encode()))
            except:
                return None
        else:
            return None
            #raise KeyError

    def get(self, key, default):
        return self[key]

    @property
    def size(self):
        return self.dn_cpp.GetSize()

    @property
    def fromConnectivity(self):
        return self.dn_cpp.GetFromConnectivity()

    @fromConnectivity.setter
    def fromConnectivity(self,val):
        self.dn_cpp.SetFromConnectivity(val)

    #@cython.boundscheck(False)  # Deactivate bounds checking
    #@cython.wraparound(False)   # Deactivate negative indexing.
    def ComputeNumberingFromConnectivity(self,mesh,notused):
        cdef cNUM.UnstructuredMesh obj = NUM.UnstructuredMesh()
        obj.SetDataFromPython(mesh)
        #cdef cNUM.NativeUnstructuredMesh* cnum_object=  (cNUM.NativeUnstructuredMesh)()
        self.dn_cpp.ComputeNumberingFromConnectivity( obj.GetCppObject())

        """
        self._doftopointLeft = range(self.size)
        self._doftopointRight = range(self.size)
        self._doftocellLeft =  []
        self._doftocellRight = []

        almanac = {}
        for i in range(self.size):
            almanac[('P', i, None)] = i
        self.almanac = almanac
        """
        return self

    def ComputeNumberingGeneral(self,mesh,space,elementFilter=None,discontinuous=False):
        self.fromConnectivity = False
        cdef cNUM.UnstructuredMesh obj = NUM.UnstructuredMesh()
        obj.SetDataFromPython(mesh)


        cdef cNS.WrapedSpace s = cNS.WrapedSpace()
        s.SetDataFromPython(space)


        cdef cNF.WrapElementFilterEvaluated ef = cNF.WrapElementFilterEvaluated()
        ef.SetIdsToTreat(mesh,elementFilter)
        self.dn_cpp.ComputeNumberingGeneral(obj.GetCppObject(), s.GetCppObject(), ef.GetCppObject())

        """
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
        """
        return self
    """
    def ApplyNumberingToMesh(self,mesh,ef,space,numbering):
        almanac = self.almanac
        cctt = 0
        for name,data,elids in ef:
            cctt += len(elids)
            res = self.GetHashFor(data,space[name],elids,discontinuous=False)
            if name in numbering:
                dofs = numbering[name]
            else:
                dofs = numbering.setdefault(name,np.zeros((data.GetNumberOfElements(),space[name].GetNumberOfShapeFunctions()), dtype=np.int_) -1)
            for i in range(len(res)):
                lres = res[i]
                ldofs = dofs[:,i]
                for j,elid in enumerate(elids):
                    if ldofs[elid] >= 0:
                        continue
                    ldofs[elid] = almanac.get(lres[j],-1)
            numbering[name] = dofs
        numbering.size = self.size

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


    def __init__(self):
        super(DofNumberingDict,self).__init__()
        self.numbering = dict()
        self._doftopointLeft = None
        self._doftopointRight= None
        self._doftocellLeft =  None
        self._doftocellRight = None
        self.almanac = dict()
        self.newAlmanac = dict()



    def __contains__(self, k):
        return k in self.numbering




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
                #""G is for global ""
                key = (on,0,idxII)
                res.append( [ key for x in elids ]  )
            elif on == "IP":
                res.append( [ (on,tuple(lcoon),idxI) for lcoon,i in zip(elidsConnectivity,elids) ]  )
            else:
                print(on)
                raise
        return res
"""
    def __str__(self):
        res = self.dn_cpp.ToStr().decode('UTF-8')
        return res

def CheckIntegrity(GUI=False):
    import BasicTools.FE.DofNumbering  as DN
    return DN.CheckIntegrityUsingAlgo("DictBase",GUI)

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
