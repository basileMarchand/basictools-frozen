# distutils: language = c++
#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, annotate=True
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np
cimport numpy as cnp
cnp.import_array()

from eigency.core cimport *


from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.SymWeakForm import testcharacter
import BasicTools.Containers.ElementNames as EN
cimport BasicTools.FE.WeakForms.NativeNumericalWeakForm as NNWF
cimport BasicTools.FE.Integrators.NativeIntegration as NI
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject


from BasicTools.CythonDefs cimport CBasicIndexType, CBasicFloatType
from BasicTools.NumpyDefs import PBasicIndexType

cimport numpy as cnp
cimport cython

from libcpp.vector cimport vector
from libcpp.string cimport string
cnp.import_array()

cdef class PyMonoElementsIntegralCpp():
    cdef NI.MonoElementsIntegralCpp NativeIntegrator
    #this are internal references to the python objects
    cdef list __ufs__
    cdef list __tfs__
    cdef list __efs__
    cdef list __ipefs__
    cdef dict __cfs__
    cdef list constantsNames

    cdef dict integrationRule

    cdef CBasicIndexType numberOfVIJ
    cdef CBasicIndexType maxNumberOfElementVIJ

    cdef cnp.ndarray nodes
    cdef dict geoSpace
    cdef int geoSpaceNumber

    cdef list __usedSpaces__
    cdef list __usedNumbering__
    cdef list __usedValues__
    cdef list __spaces_ipvalues__

    cdef cnp.ndarray testDofsNumbering
    cdef cnp.ndarray unkownDofsNumbering

    cdef int numberOfTerms

    cdef object BOO

    # we keep reference to the data send to C++

    cdef cnp.ndarray __vK
    cdef cnp.ndarray __iK
    cdef cnp.ndarray __jK
    cdef cnp.ndarray __F
    cdef list __values
    cdef list __ipvalues
    cdef cnp.ndarray __nodes
    cdef cnp.ndarray __localUnkownDofsOffset
    cdef cnp.ndarray __localtTestDofsOffset
    cdef cnp.ndarray __conn
    cdef list __valN
    cdef list __dphidxi
    cdef list __numbering

    def __init__(self):
        self.NativeIntegrator = NI.MonoElementsIntegralCpp()
        self.geoSpace = LagrangeSpaceGeo
        self.__ufs__ = list()
        self.__tfs__ = list()
        self.__efs__ = list()
        self.__ipefs__ = list()
        self.maxNumberOfElementVIJ = 0
        self.BOO = BaseOutputObject()

    def Reset(self):
        self.maxNumberOfElementVIJ = 0
        with nogil:
            self.NativeIntegrator.Reset()

    def IsMultiThread(self):
        return True

    def SetUnkownFields(self,list ufs): # OK

        self.__ufs__ = ufs
        if ufs is None:
            with nogil:
                self.NativeIntegrator.SetNumberOfUnkownFields(0)
                self.NativeIntegrator.SetTotalUnkownDofs(0)

        cdef int cpt = len(ufs)
        with nogil:
            self.NativeIntegrator.SetNumberOfUnkownFields(cpt)

        cdef int totalUnkownDofs = 0
        cpt = 0
        for uf in ufs:
            with nogil:
               self.NativeIntegrator.SetUnkownOffset(cpt,totalUnkownDofs)
            totalUnkownDofs += uf.numbering["size"]
            cpt += 1
        with nogil:
            self.NativeIntegrator.SetTotalUnkownDofs(totalUnkownDofs)

    def SetTestFields(self,list tfs = None): #ok
        if tfs is None:
            self.__tfs__ = [ uf.GetTestField() for uf in self.__ufs__ ]
        else:
            self.__tfs__ = tfs
        cdef int nbtfs
        nbtfs = len(self.__tfs__)
        with nogil:
            self.NativeIntegrator.SetNumberOfTestFields(nbtfs)
        cdef int totalTestDofs = 0
        cdef int cpt = 0
        for tf in self.__tfs__:
            with nogil:
                self.NativeIntegrator.SetTestOffset(cpt,totalTestDofs)
            totalTestDofs += tf.numbering["size"]
            cpt += 1
        with nogil:
            self.NativeIntegrator.SetTotalTestDofs(totalTestDofs)

    def SetExtraFields(self,efs): #ok
        self.__efs__ = []
        self.__ipefs__ = []
        for ef in efs:
            if isinstance(ef,FEField):
                self.__efs__.append(ef)
            else:
                self.__ipefs__.append(ef)

    def SetConstants(self,cfs):
        self.__cfs__ = cfs

        cdef int cpt =len(cfs)

        with nogil:
            self.NativeIntegrator.SetNumberOfConstants(cpt)

        cpt =0
        self.constantsNames = []
        cdef double value =0
        for name,value in cfs.items():
            self.constantsNames.append(name)
            with nogil:
                self.NativeIntegrator.SetConstants(cpt,value)
            cpt += 1

    def ComputeNumberOfVIJ(self,mesh,elementFilter)-> CBasicIndexType: # OK
       """
        Compute and return the number triplets to be calculated during integration
       """
       self.numberOfVIJ = 0
       self.maxNumberOfElementVIJ = 0
       cdef  CBasicIndexType numberOfUsedElements;
       cdef  CBasicIndexType us;
       cdef  CBasicIndexType ts;
       for name,data,ids in elementFilter:
            if len(ids) == 0:
                continue
            numberOfUsedElements = len(ids)

            us = np.sum([f.space[name].GetNumberOfShapeFunctions() for f in self.__ufs__] )
            ts = np.sum([f.space[name].GetNumberOfShapeFunctions() for f in self.__tfs__ ] )

            self.maxNumberOfElementVIJ = max(self.maxNumberOfElementVIJ,numberOfUsedElements*(us*ts))
            self.numberOfVIJ += numberOfUsedElements*(us*ts)

       return self.numberOfVIJ

    def SetIntegrationRule(self,itegrationRuleOrName):  # OK

        if itegrationRuleOrName is None :
            from BasicTools.FE.IntegrationsRules import LagrangeP1
            self.integrationRule = LagrangeP1
        elif  type(itegrationRuleOrName) == dict :
            self.integrationRule = itegrationRuleOrName
        elif  type(itegrationRuleOrName) == str :
            from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac
            self.integrationRule = IntegrationRulesAlmanac[itegrationRuleOrName]
        else:
            raise(Exception("Error seting the integration rule.."))





    def PrepareFastIntegration(self,
                               mesh,
                               wform,
                               cnp.ndarray[CBasicFloatType, ndim=1, mode="c"] vK not None ,
                               cnp.ndarray[CBasicIndexType, ndim=1, mode="c"] iK not None ,
                               cnp.ndarray[CBasicIndexType, ndim=1, mode="c"] jK not None ,
                               int cpt,
                               cnp.ndarray[CBasicFloatType, ndim=1, mode="c"] F not None):

        self.__vK = vK
        self.__iK = iK
        self.__jK = jK
        self.__F = F

        self.NativeIntegrator.vK = &vK[0]
        self.NativeIntegrator.iK = &iK[0]
        self.NativeIntegrator.jK = &jK[0]
        self.NativeIntegrator.F = &F[0]

        self.numberOfTerms = wform.GetNumberOfTerms()

        cdef bint hasnormal = False
        for monom in wform:
            for term in monom:
                if "Normal" in term.fieldName:
                    hasnormal = True
                    break
            if hasnormal == True:
                break

        with nogil:
            self.NativeIntegrator.SetComputeNormal(hasnormal)

        ## spaces treatement
        spacesId = {}
        spacesNames = {}
        spacesId[id(self.geoSpace)] = self.geoSpace
        spacesNames["Geometry_Space_internal"] = id(self.geoSpace)

        for uf in self.__ufs__:
            spacesId[id(uf.space)] = uf.space
            spacesNames[uf.name] = id(uf.space)
        for tf in self.__tfs__:
            spacesId[id(tf.space)] = tf.space
            spacesNames[tf.name] = id(tf.space)
        for ef in self.__efs__:
            spacesId[id(ef.space)] = ef.space
            spacesNames[ef.name] = id(ef.space)

        sId = list(spacesId.keys())
        self.__usedSpaces__ =  [ spacesId[k] for k in sId]

        self.geoSpaceNumber = sId.index(spacesNames["Geometry_Space_internal"])
        self.NativeIntegrator.geoSpaceNumber = self.geoSpaceNumber

        spacesNames = { sn:sId.index(spacesNames[sn]) for sn in spacesNames }

        # Numbering treatement
        numberingId = {}
        numberingNames = {}
        for uf in self.__ufs__:
            numberingId[id(uf.numbering)] = uf.numbering
            numberingNames[uf.name] = id(uf.numbering)
        for tf in self.__tfs__:
            numberingId[id(tf.numbering)] = tf.numbering
            numberingNames[tf.name] = id(tf.numbering)
        for ef in self.__efs__:
            numberingId[id(ef.numbering)] = ef.numbering
            numberingNames[ef.name] = id(ef.numbering)

        nId = list(numberingId.keys())
        self.__usedNumbering__ =  [ numberingId[k] for k in nId]

        numberingNames = { sn:nId.index(numberingNames[sn]) for sn in numberingNames}

        self.testDofsNumbering = np.array([ numberingNames[f.name] for f in self.__tfs__], dtype=PBasicIndexType)
        self.unkownDofsNumbering = np.array([ numberingNames[f.name] for f in self.__ufs__] ,dtype=PBasicIndexType)

        # Values treatement
        valuesId = {}
        valuesNames = {}
        for ef in self.__efs__:
            valuesId[id(ef.data)] = ef.data
            valuesNames[ef.name] = id(ef.data)

        vId = list(valuesId.keys())
        self.__usedValues__ =  [ valuesId[k] for k in vId]

        valuesNames = { sn:vId.index(valuesNames[sn]) for sn in valuesNames}


        for monom in wform:
            for term in monom:
                if term.fieldName in [f.name for f in self.__ufs__] :
                    term.spaceIndex_= spacesNames[term.fieldName]
                    term.numberingIndex_= numberingNames[term.fieldName]
                    #used for the offset
                    term.valuesIndex_= [uf.name for uf in  self.__ufs__ ].index(term.fieldName)
                    term.internalType = term.EnumUnknownField

                elif term.fieldName in [f.name for f in self.__tfs__]:
                    term.spaceIndex_= spacesNames[term.fieldName]
                    term.numberingIndex_= numberingNames[term.fieldName]
                    term.valuesIndex_= [uf.name for uf in  self.__tfs__ ].index(term.fieldName)
                    term.internalType = term.EnumTestField

                elif term.fieldName in [f.name for f in self.__efs__]:
                    term.spaceIndex_= spacesNames[term.fieldName]
                    term.numberingIndex_= numberingNames[term.fieldName]
                    term.valuesIndex_= valuesNames[term.fieldName]
                    term.internalType = term.EnumExtraField

                elif term.fieldName in [f.name for f in self.__ipefs__]:
                    term.valuesIndex_= [ipef.name for ipef in  self.__ipefs__ ].index(term.fieldName)
                    term.internalType = term.EnumExtraIPField

                elif "Normal" in term.fieldName :
                    term.internalType = term.EnumNormal

                elif term.constant:
                    if term.fieldName in self.constantsNames:
                        term.valuesIndex_ = self.constantsNames.index(term.fieldName)
                        term.internalType = term.EnumConstant
                    elif term.fieldName in [f.name for f in self.__efs__]:
                        term.spaceIndex_= spacesNames[term.fieldName]
                        term.numberingIndex_= numberingNames[term.fieldName]
                        term.valuesIndex_= valuesNames[term.fieldName]
                        term.internalType = term.EnumExtraField
                    elif term.fieldName in [f.name for f in self.__ipefs__]:
                        term.valuesIndex_= [ipef.name for ipef in  self.__ipefs__ ].index(term.fieldName)
                        term.internalType = term.EnumExtraIPField
                    else:
                        raise(Exception("Field : '{}' not found in Constants nor FEField nor IPFIelds".format(term.fieldName)))
                else :
                    term.internalType = term.EnumError
                    raise(Exception("Term " +str(term.fieldName) + " not found in the database " ))

        #self.SetPoints(np.ascontiguousarray(mesh.nodes))
        self.SetPoints(np.require(mesh.nodes,requirements=["C","A"]))
        #self.SetPoints(mesh.nodes)

        ########### sending values  ################################
        cdef int i = len(self.__usedValues__)
        with nogil:
            self.NativeIntegrator.SetNumberOfValues(i)
        self.__values = [None]*len(self.__usedValues__)
        self.__ipvalues = [None]*len(self.__ipefs__)
        cdef cnp.ndarray[CBasicFloatType, ndim=1, mode = 'c' ] vdata
        for i in range(len(self.__usedValues__)):
            try:
                self.SetValues(i,self.__usedValues__[i])
            except:
                print("error setting field " +str(vId[i]))
                for ef in self.__efs__:
                    print(ef.name,id(ef.data), ef.data)

                raise

    def SetValues(self,int i,cnp.ndarray[CBasicFloatType, ndim=1,mode="c"] vdata not None):
        self.__values[i] = vdata
        self.NativeIntegrator.SetValueI(i,vdata.shape[0], 1, &vdata[0])

    def SetIPValues(self,int i,cnp.ndarray[CBasicFloatType, ndim=2,mode="c"] vdata not None):
        self.__ipvalues[i] = vdata
        self.NativeIntegrator.SetIPValueI(i,vdata.shape[0],vdata.shape[1], &vdata[0,0])

    #@cython.boundscheck(False)
    #@cython.wraparound(False)
    def SetPoints(self,cnp.ndarray[CBasicFloatType, ndim=2,mode="c"] nodes not None):
        ## we make a copy in case is not continuous
        self.__nodes = nodes

        cdef CBasicIndexType m, n
        m, n = self.__nodes.shape[0], self.__nodes.shape[1]
        with nogil:
            self.NativeIntegrator.SetPoints(&nodes[0,0],m,n)

        return None

    cpdef SetOnlyEvaluation(self, onlyEvaluation: bool) :
        self.NativeIntegrator.onlyEvaluation = onlyEvaluation

    def ActivateElementType(self, domain):

        elementType = domain.elementType

        if not elementType in EN.geoSupport :
            print("Dont know this element : ", elementType)

        if not elementType in self.integrationRule :
            print("Integration rule incomplete for this type of geo element : "+  str(elementType) + "  integration rule : "+ str(list(self.integrationRule.keys())))

        p, w = self.integrationRule[elementType]
        cdef CBasicIndexType numberOfIntegrationPoints = len(w)


        ufssize = [f.space[elementType].GetNumberOfShapeFunctions() for f in self.__ufs__]
        tfssize = [f.space[elementType].GetNumberOfShapeFunctions() for f in self.__tfs__]

        self.__localUnkownDofsOffset = np.zeros(len(ufssize),dtype=PBasicIndexType)
        self.__localtTestDofsOffset = np.zeros(len(tfssize),dtype=PBasicIndexType)
        su = 0
        for i in range(len(ufssize)):
            self.__localUnkownDofsOffset[i] = su
            su += ufssize[i]

        st = 0.
        for i in range(len(tfssize)):
            self.__localtTestDofsOffset[i] = st
            st += tfssize[i]

        self.NativeIntegrator.SetLocalOffsets(su,self.__localUnkownDofsOffset,self.unkownDofsNumbering,st,self.__localtTestDofsOffset,self.testDofsNumbering)

        with nogil:
            self.NativeIntegrator.SetNumberOfIntegrationPoints(numberOfIntegrationPoints)

        for i in range(numberOfIntegrationPoints):
            if p.shape[1] == 1:
                self.NativeIntegrator.SetIntegrationPointI(i, w[i],p[i,0],0.,0.)
            elif p.shape[1] == 2:
                self.NativeIntegrator.SetIntegrationPointI(i, w[i],p[i,0],p[i,1],0.)
            elif p.shape[1] == 3:
                self.NativeIntegrator.SetIntegrationPointI(i, w[i],p[i,0],p[i,1],p[i,2])
            else:
                raise(Exception("error"))

        # sending the connectivity matrix for this elements
        cdef int m, n
        m = domain.connectivity.shape[0]
        n = domain.connectivity.shape[1]

        cdef cnp.ndarray[CBasicIndexType, ndim=2, mode="c"] conn = domain.connectivity
        self.__conn = conn

        with nogil:
            self.NativeIntegrator.SetConnectivity(&conn[0,0],m,n)

        self.NativeIntegrator.SetNumberOfSpaces(len(self.__usedSpaces__))
        cdef cnp.ndarray[CBasicFloatType , ndim=1, mode = 'c' ] valN
        cdef cnp.ndarray[CBasicFloatType , ndim=2, mode = 'c' ] dphidxi

        #self.BOO.PrintDebug("InitSpaceS ")

        self.__valN = [None]*len(self.__usedSpaces__)
        self.__dphidxi = [None]*len(self.__usedSpaces__)

        self.__spaces_ipvalues__ = []
        for s in range(len(self.__usedSpaces__)):
            current_space = self.__usedSpaces__[s][elementType]
            space_ipvalues = current_space.SetIntegrationRule(p,w)
            self.__spaces_ipvalues__.append(space_ipvalues)
            dimensionality = current_space.GetDimensionality()
            NOSF = current_space.GetNumberOfShapeFunctions()

            self.NativeIntegrator.InitSpaceS(s,dimensionality,NOSF,numberOfIntegrationPoints)

            self.__valN[s] = [None]*numberOfIntegrationPoints
            self.__dphidxi[s] = [None]*numberOfIntegrationPoints

            for i in range(numberOfIntegrationPoints):
                valN = space_ipvalues.valN[i]
                self.__valN[s][i] = valN
                self.NativeIntegrator.SetSpaceSvalNI(s,i,&valN[0])
                dphidxi = space_ipvalues.valdphidxi[i]
                self.__dphidxi[s][i] = dphidxi
                if dphidxi.size != 0:
                    self.NativeIntegrator.SetSpaceSvaldphidxiI(s,i,&dphidxi[0,0])
        ############ sending lnumbering ###############################
        self.NativeIntegrator.SetNumberOfNumberings(len(self.__usedNumbering__))
        cdef cnp.ndarray[CBasicIndexType , ndim=2, mode = 'c' ] ndata
        #print(self.__usedNumbering__[0][elementType].flags)

        self.__numbering = [None] * len(self.__usedNumbering__)
        for i in range(len(self.__usedNumbering__)):
            ndata =  self.__usedNumbering__[i].get(elementType,None)
            if ndata is None:
                continue
            self.__numbering[i] = ndata
            self.NativeIntegrator.SetNumberingI(i,ndata.shape[0],ndata.shape[1], &ndata[0,0])

        self.NativeIntegrator.SetNumberOfIPValues(len(self.__ipefs__))
        for i in range(len(self.__ipefs__)):
            self.SetIPValues(i,self.__ipefs__[i].data[domain.elementType] )

    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    def Integrate(self,NNWF.PyWeakForm wform,
                  cnp.ndarray[CBasicIndexType, ndim=1, mode="c"] idstotreat not None ):

        cdef NNWF.WeakForm* wfobject =  wform.GetCppPointer()
        cdef CBasicIndexType idsize = idstotreat.shape[0]
        cdef CBasicIndexType* idpointer = &idstotreat[0]
        with nogil:
            self.NativeIntegrator.Integrate(wfobject, idsize, idpointer)

    def GetTotalTestDofs(self):
        return self.NativeIntegrator.totalTestDofs

    def GetTotalUnkownDofs(self):
        return self.NativeIntegrator.totalUnkownDofs

    def GetNumberOfUsedIvij(self):
        return self.NativeIntegrator.totalvijcpt

    def AddToNumbefOfUsedIvij(self,data):
        self.NativeIntegrator.totalvijcpt += data

def CheckIntegrity(GUI=False):
    obj = PyMonoElementsIntegralCpp()
    print(obj)
    return "OK"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity())
