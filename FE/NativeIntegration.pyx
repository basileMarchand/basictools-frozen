# distutils: language = c++

##https://stackoverflow.com/questions/21851985/difference-between-np-int-np-int-int-and-np-int-t-in-cython
import numpy as np

from BasicTools.FE.Fields.NodalField import NodalField
from BasicTools.FE.WeakForm import testcharacter
import BasicTools.Containers.ElementNames as EN
#from BasicTools.FE.WeakFormNumerical cimport PyWeakForm
cimport BasicTools.FE.WeakFormNumerical as WFN
cimport BasicTools.FE.WeakFormNumericalWrapper as WFNW
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject


cimport numpy as cnp
cimport cython

int_DTYPE   = np.int
float_DTYPE = np.float

ctypedef cnp.int_t     int_DTYPE_t
ctypedef cnp.float64_t float_DTYPE_t

from libcpp.vector cimport vector
from libcpp.string cimport string
np.import_array()

cdef class PyMonoElementsIntegralCpp():
    cdef MonoElementsIntegralCpp NativeIntegrator
    #this are internal references to the python objects
    cdef list __ufs__
    cdef list __tfs__
    cdef list __efs__
    cdef dict __cfs__
    cdef list constantsNames

    cdef dict integrationRule

    cdef int numberOfVIJ
    cdef int maxNumberOfElementVIJ

#    cdef object mesh
    cdef np.ndarray nodes
    cdef dict geoSpace
    cdef int geoSpaceNumber

    cdef list __usedSpaces__
    cdef list __usedNumbering__
    cdef list __usedValues__

    cdef np.ndarray testDofsNumbering;
    cdef np.ndarray unkownDofsNumbering;

    cdef int numberOfTerms

    cdef object BOO

    # we keep reference to the data send to C++

    cdef np.ndarray __vK
    cdef np.ndarray __iK
    cdef np.ndarray __jK
    cdef np.ndarray __F
    cdef list __values
    cdef np.ndarray __nodes
    cdef np.ndarray __localUnkownDofsOffset
    cdef np.ndarray __localtTestDofsOffset
    cdef np.ndarray __conn
    cdef list __valN
    cdef list __dphidxi
    cdef list __numbering

    def __init__(self):
        self.NativeIntegrator = MonoElementsIntegralCpp()
        self.geoSpace = LagrangeSpaceGeo
        self.__ufs__ = list()
        self.__tfs__ = list()
        self.__efs__ = list()
        self.maxNumberOfElementVIJ = 0
        self.BOO = BaseOutputObject()




    def SetUnkownFields(self,list ufs): # OK
      self.__ufs__ = ufs

      self.NativeIntegrator.SetNumberOfUnkownFields(len(ufs));
      cdef int totalUnkownDofs = 0
      cdef int cpt = 0
      for uf in ufs:
        self.NativeIntegrator.SetUnkownOffset(cpt,totalUnkownDofs);
        totalUnkownDofs += uf.numbering["size"]
        cpt += 1
      self.NativeIntegrator.SetTotalUnkownDofs(totalUnkownDofs);

    def SetTestFields(self,list tfs = None): #ok
      if tfs is None:
         self.__tfs__  = []
         for f in self.__ufs__:
            self.__tfs__.append(NodalField(name=f.name+testcharacter,mesh=f.mesh,space=f.space,numbering=f.numbering,data=f.data) )
      else:
          self.__tfs__ = tfs

      self.NativeIntegrator.SetNumberOfTestFields(len(self.__tfs__));
      cdef int totalTestDofs = 0
      cdef int cpt = 0
      for tf in self.__tfs__:
        self.NativeIntegrator.SetTestOffset(cpt,totalTestDofs);
        totalTestDofs += tf.numbering["size"]
        cpt += 1
      self.NativeIntegrator.SetTotalTestDofs(totalTestDofs);

    def SetExtraFields(self,efs): #ok
        self.__efs__ = efs

    def SetConstants(self,cfs):
        self.__cfs__ = cfs


        self.NativeIntegrator.SetNumberOfConstants(len(cfs))

        self.constantsNames = []
        cdef int cpt =0
        for name,value in cfs.items():
            self.constantsNames.append(name)
            self.NativeIntegrator.SetConstants(cpt,value)

    def ComputeNumberOfVIJ(self,mesh,tag): # OK
       #Total number of ikv calculated
       self.numberOfVIJ = 0
       self.maxNumberOfElementVIJ = 0
       for name,data in mesh.elements.items():
          if tag == "ALL":
             numberOfUsedElements = data.GetNumberOfElements()
          elif tag in data.tags:
             numberOfUsedElements = len(data.tags[tag])
          else:
             continue

          us = np.sum([f.space[name].GetNumberOfShapeFunctions() for f in self.__ufs__] )
          ts = np.sum([f.space[name].GetNumberOfShapeFunctions() for f in self.__tfs__ ] )

          self.maxNumberOfElementVIJ = max(self.maxNumberOfElementVIJ,numberOfUsedElements*(us*ts))
          self.numberOfVIJ += numberOfUsedElements*(us*ts)

       return self.numberOfVIJ

    def SetIntegrationRule(self,itegrationRuleOrName=None):  # OK

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
                               np.ndarray[float_DTYPE_t, ndim=1, mode="c"] vK not None ,
                               np.ndarray[int_DTYPE_t, ndim=1, mode="c"] iK not None ,
                               np.ndarray[int_DTYPE_t, ndim=1, mode="c"] jK not None ,
                               int cpt,
                               np.ndarray[float_DTYPE_t, ndim=1, mode="c"] F not None):

      self.__vK = vK
      self.__iK = iK
      self.__jK = jK
      self.__F = F

      self.NativeIntegrator.vK = &vK[0]
      self.NativeIntegrator.iK = &iK[0]
      self.NativeIntegrator.jK = &jK[0]
      self.NativeIntegrator.F = &F[0]

      self.numberOfTerms = wform.GetNumberOfTerms();

      cdef bint hasnormal = False
      for monom in wform:
         for term in monom:
             if "Normal" in term.fieldName:
                 hasnormal = True
                 break
         if hasnormal == True:
             break

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

      self.testDofsNumbering = np.array([ numberingNames[f.name] for f in self.__tfs__], dtype=int)
      self.unkownDofsNumbering = np.array([ numberingNames[f.name] for f in self.__ufs__] ,dtype=int)

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
            if "Normal" in term.fieldName :
                term.internalType = 0
            elif term.constant:
                try:
                    term.valuesIndex_ = self.constantNames.index(term.fieldName)
                    term.internalType = 1
                except:
                    #print("Warning: constant '" +str(term.fieldName) + "'  not found in constants ")
                    #print("searching extra fields")
                    term.spaceIndex_= spacesNames[term.fieldName]
                    term.numberingIndex_= numberingNames[term.fieldName]
                    term.valuesIndex_= valuesNames[term.fieldName]
                    term.internalType = 4

            elif term.fieldName in [f.name for f in self.__ufs__] :
                term.spaceIndex_= spacesNames[term.fieldName]
                term.numberingIndex_= numberingNames[term.fieldName]
                #used for the offset
                term.valuesIndex_= [uf.name for uf in  self.__ufs__ ].index(term.fieldName)

                term.internalType = 2
            elif term.fieldName in [f.name for f in self.__tfs__]:
                term.spaceIndex_= spacesNames[term.fieldName]
                term.numberingIndex_= numberingNames[term.fieldName]
                term.valuesIndex_= [uf.name for uf in  self.__tfs__ ].index(term.fieldName)

                term.internalType = 3
            elif term.fieldName in [f.name for f in self.__efs__]:
                term.spaceIndex_= spacesNames[term.fieldName]
                term.numberingIndex_= numberingNames[term.fieldName]
                term.valuesIndex_= valuesNames[term.fieldName]
                term.internalType = 4
            else :
                term.internalType = -1
                raise(Exception("Term " +str(term.fieldName) + " not found in the database " ))

      self.SetPoints(mesh.nodes)

      ########### sending values  ################################
      self.NativeIntegrator.SetNumberOfValues(len(self.__usedValues__))
      self.__values = [None]*len(self.__usedValues__)
      cdef np.ndarray[float_DTYPE_t, ndim=1, mode = 'c' ] vdata
      for i in range(len(self.__usedValues__)):
          self.SetValues(i,self.__usedValues__[i])




    def SetValues(self,int i,np.ndarray[float_DTYPE_t, ndim=1,mode="c"] vdata not None):
          self.__values[i] = vdata
          self.NativeIntegrator.SetValueI(i,vdata.shape[0],vdata.shape[1], &vdata[0])

    #@cython.boundscheck(False)
    #@cython.wraparound(False)
    def SetPoints(self,np.ndarray[float_DTYPE_t, ndim=2,mode="c"] nodes not None):
        ## we make a copy in case is not continuous
        self.__nodes = nodes

        cdef int m, n
        m, n = self.__nodes.shape[0], self.__nodes.shape[1]
        self.NativeIntegrator.SetPoints(&nodes[0,0],m,n)

        return None

    def SetOnlyEvaluation(self, onlyEvaluation):
        self.NativeIntegrator.onlyEvaluation = onlyEvaluation

    def ActivateElementType(self, domain):

      elementType = domain.elementType
      p, w = self.integrationRule[EN.geoSupport[elementType]];
      cdef int numberOfIntegrationPoints = len(w)


      ufssize = [f.space[elementType].GetNumberOfShapeFunctions() for f in self.__ufs__]
      tfssize = [f.space[elementType].GetNumberOfShapeFunctions() for f in self.__tfs__]

      self.__localUnkownDofsOffset = np.zeros(len(ufssize))
      self.__localtTestDofsOffset = np.zeros(len(tfssize))
      su = 0.
      for i in range(len(ufssize)):
          self.__localUnkownDofsOffset[i] = su
          su += ufssize[i]

      st = 0.
      for i in range(len(tfssize)):
          self.__localtTestDofsOffset[i] = st
          st += tfssize[i]

      self.NativeIntegrator.SetLocalOffsets(su,self.__localUnkownDofsOffset,self.unkownDofsNumbering,st,self.__localtTestDofsOffset,self.testDofsNumbering)


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

      cdef cnp.ndarray[int_DTYPE_t, ndim=2, mode="c"] conn = domain.connectivity
      self.__conn = conn

      self.NativeIntegrator.SetConnectivity(&conn[0,0],m,n)

      self.NativeIntegrator.SetNumberOfSpaces(len(self.__usedSpaces__))
      cdef cnp.ndarray[float_DTYPE_t , ndim=1, mode = 'c' ] valN
      cdef cnp.ndarray[float_DTYPE_t , ndim=2, mode = 'c' ] dphidxi

      #self.BOO.PrintDebug("InitSpaceS ")

      self.__valN = [None]*len(self.__usedSpaces__)
      self.__dphidxi = [None]*len(self.__usedSpaces__)

      for s in range(len(self.__usedSpaces__)):
          current_space = self.__usedSpaces__[s][elementType]
          current_space.SetIntegrationRule(p,w)

          dimensionality = current_space.GetDimensionality()
          NOSF = current_space.GetNumberOfShapeFunctions()

          self.NativeIntegrator.InitSpaceS(s,dimensionality,NOSF,numberOfIntegrationPoints)

          self.__valN[s] = [None]*numberOfIntegrationPoints
          self.__dphidxi[s] = [None]*numberOfIntegrationPoints

          for i in range(numberOfIntegrationPoints):
              valN = current_space.valN[i]
              self.__valN[s][i] = valN
              self.NativeIntegrator.SetSpaceSvalNI(s,i,&valN[0])
              dphidxi = current_space.valdphidxi[i]
              self.__dphidxi[s][i] = dphidxi

              self.NativeIntegrator.SetSpaceSvaldphidxiI(s,i,&dphidxi[0,0])
      ############ sending lnumbering ###############################
      self.NativeIntegrator.SetNumberOfNumberings(len(self.__usedNumbering__))
      cdef cnp.ndarray[int_DTYPE_t , ndim=2, mode = 'c' ] ndata
      #print(self.__usedNumbering__[0][elementType].flags)

      self.__numbering = [None] * len(self.__usedNumbering__)
      for i in range(len(self.__usedNumbering__)):
          ndata =  self.__usedNumbering__[i].get(elementType,None)
          if ndata is None:
              continue
          self.__numbering[i] = ndata
          self.NativeIntegrator.SetNumberingI(i,ndata.shape[0],ndata.shape[1], &ndata[0,0])


    def Integrate(self,WFN.PyWeakForm wform,
                  np.ndarray[int_DTYPE_t, ndim=1, mode="c"] idstotreat not None ):

        self.NativeIntegrator.Integrate(<WFNW.WeakForm*> wform.GetCppObject(), idstotreat)

    def GetTotalTestDofs(self):
      return self.NativeIntegrator.totalTestDofs

    def GetTotalUnkownDofs(self):
      return self.NativeIntegrator.totalUnkownDofs

    def GetNumberOfUsedIvij(self):
        return self.NativeIntegrator.totalvijcpt

def CheckIntegrity(GUI=False):
    return "OK"
