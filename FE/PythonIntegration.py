# -*- coding: utf-8 -*-

# distutils: language = c++

import numpy as np
from scipy.sparse import coo_matrix

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO

import BasicTools.Containers.ElementNames as EN

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.WeakForm import testcharacter
from BasicTools.FE.Fields.NodalField import NodalField


class MonoElementsIntegral(BOO):

  def __init__(self):
      super(MonoElementsIntegral,self).__init__()
      #self.unkownDofsOffset = np.zeros(1,dtype=DTYPEint_t)
      self.unkownDofsOffset = None
      self.testDofsOffset = None
      self.totalTestDofs = 0
      self.totalUnkownDofs= 0
      self.geoSpace = None
      self.__cfs__ = None

      #self.constantsNumerical = None

      ##self.unknowNames = None
      #self.testNames = None
      #self.extraFieldsNames = None
      self.integrationRule = None
      self.numberOfVIJ = 0
      self.totalvijcpt = 0
      self.maxNumberOfTerms = 0;


#  cdef public ar unkownDofsOffset
#  cdef public np.ndarray testDofsOffset
#  cdef int totalTestDofs
#  cdef int totalUnkownDofs
#  cdef object geoSpace
#  cdef int geoSpaceNumber
#  cdef int totalvijcpt
#  #cdef ar constantsNumerical
#  cdef dict integrationRule
#  cdef list __ufs__
#  cdef list __tfs__
#  cdef list __efs__
#  cdef dict __cfs__
#  cdef int numberOfVIJ
#  cdef bint hasnormal
#  cdef bint onlyEvaluation
#
#  cdef list __usedSpaces__
#  cdef list __usedNumbering__
#  cdef list __usedValues__
#  cdef np.ndarray nodes
#
#  #cdef object mesh
#
#  cdef np.ndarray connectivity
#
#  cdef np.ndarray vK#[DTYPEfloat_t, ndim=1]
#  cdef np.ndarray iK#[DTYPEint_t, ndim=1]
#  cdef np.ndarray jK #[DTYPEint_t, ndim=1]
#  cdef np.ndarray F
#
#  cdef np.ndarray w
#  cdef np.ndarray p
#  cdef list localNumbering
#  cdef list localSpaces
#  cdef np.ndarray NumberOfShapeFunctionForEachSpace
#


  def SetUnkownFields(self,ufs):
      self.__ufs__ = ufs

      self.unkownDofsOffset = np.zeros(len(ufs),dtype=int)
      self.totalUnkownDofs = 0
      cpt = 0
      for uf in ufs:
        self.unkownDofsOffset[cpt] = self.totalUnkownDofs
        self.totalUnkownDofs += uf.numbering["size"]
        cpt += 1

  def SetTestFields(self,tfs=None):
      if tfs is None:
         tfs = []
         for f in self.__ufs__:
            tfs.append(NodalField(name=f.name+testcharacter,mesh=f.mesh,space=f.space,numbering=f.numbering,data=f.data) )

      self.__tfs__ = tfs

      self.testDofsOffset = np.zeros(len(tfs),dtype=int)
      self.totalTestDofs = 0
      cpt =0
      for tf in tfs:
        self.testDofsOffset[cpt] = self.totalTestDofs
        self.totalTestDofs += tf.numbering["size"]
        cpt += 1

  def SetExtraFields(self,efs):
        self.__efs__ = efs

  def SetConstants(self,cfs):
      self.__cfs__ = cfs

  def ComputeNumberOfVIJ(self,mesh,tag):
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

        self.maxNumberOfElementVIJ = max(self.maxNumberOfElementVIJ,(us*ts))
        self.numberOfVIJ += numberOfUsedElements*(us*ts)
    return self.numberOfVIJ

  def SetIntegrationRule(self,itegrationRuleOrName=None):

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



  def PrepareFastIntegration(self,mesh,wform,vK,iK,jK,cpt,F):

      self.vK = vK
      self.iK = iK
      self.jK = jK
      self.F = F

      ##we modified the internal structure (varialbe ending with _) for fast access

      self.hasnormal = False
      for monom in wform:
         for term in monom:
             if "Normal" in term.fieldName:
                 self.hasnormal = True
                 break
         if self.hasnormal == True:
             break

      constantNames = []
      for x in self.__cfs__:
          constantNames.append(x)

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

      # Values treatement
      valuesId = {}
      valuesNames = {}
      for ef in self.__efs__:
          valuesId[id(ef.data)] = ef.data
          valuesNames[ef.name] = id(ef.data)

      vId = list(valuesId.keys())
      self.__usedValues__ =  [ valuesId[k] for k in vId]

      valuesNames = { sn:vId.index(valuesNames[sn]) for sn in valuesNames}

      self.maxNumberOfTerms = 0;
      for monom in wform:
        self.maxNumberOfTerms = max(self.maxNumberOfTerms, monom.GetNumberOfProds())
        for term in monom:
            if "Normal" in term.fieldName :
                term.internalType = 0
            elif term.constant:
                try:
                    term.valuesIndex_ = constantNames.index(term.fieldName)
                    term.internalType = 1
                except:
                    print("Warning: constant '" +str(term.fieldName) + "'  not found in constants ")
                    print("searching extra fields")
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
                #term.valuesIndex_= valuesNames[term.fieldName]
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

  def SetPoints(self,nodes):
        """
        ## from https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC

        multiply (arr, value)

        Takes a numpy arry as input, and multiplies each elemetn by value, in place

        param: array -- a 2-d numpy array of np.float64

        """
        self.nodes = nodes

        return None

  def SetOnlyEvaluation(self,onlyEvaluation):
      self.onlyEvaluation = onlyEvaluation

  def ActivateElementType(self,domain):

    self.localNumbering = []
    for numbering in self.__usedNumbering__:
        self.localNumbering.append( numbering.get(domain.elementType,None) )


    self.p, self.w = self.integrationRule[EN.geoSupport[domain.elementType]];

    self.geoSpace = LagrangeSpaceGeo[domain.elementType]
    self.geoSpace.SetIntegrationRule(self.p,self.w)

    self.NumberOfShapeFunctionForEachSpace = np.zeros(len(self.__usedSpaces__), dtype=int)

    cpt = 0
    self.localSpaces = list()
    for space in self.__usedSpaces__:
        if space is None :
            self.NumberOfShapeFunctionForEachSpace[cpt] = 0
            self.localSpaces.append(None)
        else:
            self.localSpaces.append(space[domain.elementType])
            space[domain.elementType].SetIntegrationRule(self.p,self.w)
            self.NumberOfShapeFunctionForEachSpace[cpt] = space[domain.elementType].GetNumberOfShapeFunctions()
        cpt += 1

    self.connectivity = domain.connectivity

  def Integrate(self,wform,idstotreat):

    constantsNumerical = np.empty(len(self.__cfs__))
    cpt =0;
    for x in self.__cfs__:
        constantsNumerical[cpt] = self.__cfs__[x]
        cpt += 1

    NumberOfIntegrationPoints = len(self.w)


    ev = np.empty(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms()*NumberOfIntegrationPoints,dtype=np.float)
    ei = np.empty(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms()*NumberOfIntegrationPoints,dtype=np.int)
    ej = np.empty(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms()*NumberOfIntegrationPoints,dtype=np.int)

    numberOfFields = len(self.__usedSpaces__)

    BxByBzI = [None] *numberOfFields

    NxNyNzI = [None] *numberOfFields



    for n in idstotreat:
        fillcpt =0

        xcoor = self.nodes[self.connectivity[n],:]


        for ip in range(NumberOfIntegrationPoints):
            """ we recover the jacobian matrix """
            Jack, Jdet, Jinv = self.geoSpace.GetJackAndDetI(ip,xcoor)
            #print("Jack")
            #print(Jack)
            #print("Jdet")
            #print(Jdet)
            #print("Jinv")
            #print(Jinv(np.identity(Jack.shape[0])))

            for i in range(numberOfFields):
                 if self.localSpaces[i] is not None:
                     NxNyNzI[i] = self.localSpaces[i].valN[ip]

                     BxByBzI[i] = Jinv(self.localSpaces[i].valdphidxi[ip])

            if self.hasnormal:
                normal = self.geoSpace.GetNormal(Jack)

            for monom in wform:

                factor = monom.prefactor
                if self.onlyEvaluation :
                    # For the evaluation we only add the constribution without doing the integration
                    # the user is responsible of dividing by the mass matrix to get the correct values
                    # also the user can use a discontinues field to generate element surface stress (for example)
                    pass
                else:
                    # for the integration we multiply by the deteminant of the jac
                    factor *= Jdet

                hasright = False

                for term in monom:

                    if term.internalType == 0 :
                        factor *= normal[term.derDegree]
                        continue
                    elif  term.internalType == 1 :
                        factor *= constantsNumerical[term.valuesIndex_]
                        continue
                    elif  term.internalType == 2 :
                        if term.derDegree == 1:
                            right = BxByBzI[term.spaceIndex_][[term.derCoordIndex_],:]
                        else:
                            right = np.array([NxNyNzI[term.spaceIndex_],])
                        rightNumbering = self.localNumbering[term.numberingIndex_][n,:] + self.unkownDofsOffset[term.valuesIndex_]
                        hasright = True
                        l2 = self.NumberOfShapeFunctionForEachSpace[term.spaceIndex_]
                        continue
                    elif  term.internalType == 3 :
                        if term.derDegree == 1:
                            left = BxByBzI[term.spaceIndex_][term.derCoordIndex_]
                        else:
                            left = NxNyNzI[term.spaceIndex_]
                        leftNumbering = self.localNumbering[term.numberingIndex_][n,:] + self.testDofsOffset[term.valuesIndex_]
                        l1 = self.NumberOfShapeFunctionForEachSpace[term.spaceIndex_]
                        continue
                    elif  term.internalType == 4 :

                        if term.derDegree == 1:
                            func = BxByBzI[term.spaceIndex_][term.derCoordIndex_]
                        else:
                            func = NxNyNzI[term.spaceIndex_]
                        centerNumbering = self.localNumbering[term.numberingIndex_][n,:]
                        vals = self.__usedValues__[term.valuesIndex_][centerNumbering]
                        factor *= np.dot(func,vals)

                        continue
                    else :
                        raise(Exception("Cant treat term " + str(term.fieldName)))

                if factor == 0:
                    continue


                factor *= self.w[ip]


                if hasright:
                    l = l1*l2

                    l2cpt = fillcpt
                    for i in range(l1):
                        for j in range(l2) :
                            ev[l2cpt] =  left[i]*right[0,j]*factor
                            l2cpt +=1

                    l2cpt = fillcpt
                    for i in range(l1):
                        for j in range(l2) :
                          ej[l2cpt] = rightNumbering[j]
                          l2cpt += 1

                    l2cpt = fillcpt
                    for j in range(l2) :
                         for i in range(l1):
                              ei[l2cpt] = leftNumbering[j]
                              l2cpt += 1
                    fillcpt += l
                else:
                    for i in range(l1):
                      self.F[leftNumbering[i]] += left[i]*factor


        if fillcpt:
            data = coo_matrix((ev[:fillcpt], (ei[:fillcpt],ej[:fillcpt])), shape=( self.totalTestDofs,self.totalUnkownDofs))
            data.sum_duplicates()
            start = self.totalvijcpt
            stop = start+len(data.data)

            self.vK[start:stop] = data.data
            self.iK[start:stop] = data.row
            self.jK[start:stop] = data.col
            self.totalvijcpt += len(data.data)


  def GetNumberOfUsedIvij(self):
    return self.totalvijcpt

  def GetTotalTestDofs(self):
     return self.totalTestDofs

  def GetTotalUnkownDofs(self):
     return self.totalUnkownDofs



def CheckIntegrity():

    return 'ok'

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity())
