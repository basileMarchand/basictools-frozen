# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


# distutils: language = c++

import numpy as np
from scipy.sparse import coo_matrix

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO,froze_it

import BasicTools.Containers.ElementNames as EN

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.SymWeakForm import testcharacter
from BasicTools.FE.Fields.FEField import FEField

@froze_it
class MonoElementsIntegral(BOO):
    """
    Class to assembly a formulation (weak form) into a matrix or a vector

    * unkownDofsOffset : offset for the unkowns dofs
    * __ufs__          : unkown fields
    * testDofsOffset   : offset for the test dofs
    * __tfs__          : unkown dofs

    * totalTestDofs    : Total number fo test dofs    (computed Automaticaly)
    * totalUnkownDofs  : Total number fo unkown dofs  (computed Automaticaly)
    * geoSpace         : Geometry aproximation space  (computed Automaticaly)

    * __efs__          : Extra Fields
    * __cfs__          : (dic(str:float) ) Constants
    * integrationRule  : integration rule for the integratoin
    * onlyEvaluation   : To force the integrator to not multiply by the detJac

    * numberOfVIJ = 0
    * F                : rhs vector
    * vK, iK, jK       : Vectors containing the values and indices for the entries
                         of the operator
    * totalvijcpt      : Number of non zero entries in the self.vK iK ant jK vector
    * maxNumberOfTerms : maximal number of terms in a monom (computed Automaticaly)
    * maxNumberOfElementVIJ : (computed Automaticaly)
    * hasnormal        : (computed Automaticaly)
    * __usedSpaces__
    * geoSpaceNumber
    """
    def __init__(self):
        super(MonoElementsIntegral,self).__init__()
        self.unkownDofsOffset = None
        self.__ufs__ = None

        self.testDofsOffset = None
        self.__tfs__ = None

        self.totalTestDofs = 0
        self.totalUnkownDofs= 0
        self.geoSpace = None
        self.__efs__ = None
        self.__ipefs__ = None
        self.__cfs__ = {}
        self.integrationRule = None
        self.onlyEvaluation = False
        """
        For the evaluation we only add the constribution without doing the integration
        the user is responsible of dividing by the mass matrix to get the correct values
        also the user can use a discontinues field to generate element surface stress (for example)
        """
        self.numberOfVIJ = 0

        self.F  = None
        self.vK = None
        self.iK = None
        self.jK = None
        self.totalvijcpt = 0
        self.maxNumberOfTerms = 0
        self.maxNumberOfElementVIJ = 0
        self.hasnormal = False
        self.__usedSpaces__ = None
        self.__usedNumbering__ = None
        self.__usedValues__ = None
        self.__usedValuesAtIP__ = None
        self.geoSpaceNumber = 0

        # internal variables dependent on the courrent element type been treatd
        # (internal use only )
        self.localSpaces = None
        self.localNumbering = None
        self.NumberOfShapeFunctionForEachSpace = None
        self.p = None
        self.w = None
        self.nodes = None
        self.connectivity = None

    def SetUnkownFields(self,ufs):
        """
        Set the fields used for the unkown space

        ufs : list(FEField) list of fields
        """
        self.__ufs__ = ufs

        self.unkownDofsOffset = np.zeros(len(ufs),dtype=int)
        self.totalUnkownDofs = 0
        cpt = 0
        for uf in ufs:
          self.unkownDofsOffset[cpt] = self.totalUnkownDofs
          self.totalUnkownDofs += uf.numbering["size"]
          cpt += 1

    def SetTestFields(self,tfs=None):
      """
      Set the fields used for the test space

      tfs : list(FEField) list of fields
      if tfs is none then the unkown fields are used (Galerkin projection)
      """
      if tfs is None:
         tfs = []
         for f in self.__ufs__:
            tfs.append(FEField(name=f.name+testcharacter,mesh=f.mesh,space=f.space,numbering=f.numbering,data=f.data) )

      self.__tfs__ = tfs

      self.testDofsOffset = np.zeros(len(tfs),dtype=int)
      self.totalTestDofs = 0
      cpt =0
      for tf in tfs:
          self.testDofsOffset[cpt] = self.totalTestDofs
          self.totalTestDofs += tf.numbering["size"]
          cpt += 1

    def SetExtraFields(self,efs):
        """
        Set the extra fields used in the weak formulation

        efs : list(FEField or IPField) list of fields

        """
        self.__efs__ = []
        self.__ipefs__ = []
        for ef in efs:
            if isinstance(ef,FEField):
                self.__efs__.append(ef)
            else:
                self.__ipefs__.append(ef)

    def SetConstants(self,cfs):
        """
        Set The constants used in the weak formulation

        cfs : dic(str:float) constants dictionary
        """
        self.__cfs__ = cfs

    def ComputeNumberOfVIJ(self,mesh,elementFilter):
        """
        Compute and return the number triplets to be calculated during integration
        """
        self.numberOfVIJ = 0
        self.maxNumberOfElementVIJ = 0
        for name,data,ids in elementFilter:
            if len(ids) == 0:
                continue
            numberOfUsedElements = len(ids)

            us = np.sum([f.space[name].GetNumberOfShapeFunctions() for f in self.__ufs__] )
            ts = np.sum([f.space[name].GetNumberOfShapeFunctions() for f in self.__tfs__ ] )

            self.maxNumberOfElementVIJ = max(self.maxNumberOfElementVIJ,numberOfUsedElements*(us*ts))
            self.numberOfVIJ += numberOfUsedElements*(us*ts)
        return self.numberOfVIJ

    def SetIntegrationRule(self,itegrationRuleOrName=None):
        """
        Function to set the integration Rule
        """
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
      """
      Function to prepare the integration procedure, this function checks:
          - if the weak form needs the normal at each integration point
          - prepare the fields to be used
          - fills each term in the weak formulation with the data about the
            fields for fast acces

      mesh : a mesh
      wform: the weak form to be integrated
      vK,iK,jK = the vectors to store the calculated values for the K op
      cpt
      """

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
                term.internalType = term.EnumNormal
            elif term.constant:
                try:
                    term.valuesIndex_ = constantNames.index(term.fieldName)
                    term.internalType = term.EnumConstant
                except:
                    print("Warning: constant '" +str(term.fieldName) + "'  not found in constants ")
                    print("searching extra fields")
                    term.spaceIndex_= spacesNames[term.fieldName]
                    term.numberingIndex_= numberingNames[term.fieldName]
                    term.valuesIndex_= valuesNames[term.fieldName]
                    term.internalType = term.EnumExtraField

            elif term.fieldName in [f.name for f in self.__ufs__] :
                term.spaceIndex_= spacesNames[term.fieldName]
                term.numberingIndex_= numberingNames[term.fieldName]
                #used for the offset
                term.valuesIndex_= [uf.name for uf in  self.__ufs__ ].index(term.fieldName)

                term.internalType = term.EnumUnknownField
            elif term.fieldName in [f.name for f in self.__tfs__]:
                term.spaceIndex_= spacesNames[term.fieldName]
                term.numberingIndex_= numberingNames[term.fieldName]
                #term.valuesIndex_= valuesNames[term.fieldName]
                term.valuesIndex_= [uf.name for uf in  self.__tfs__ ].index(term.fieldName)

                term.internalType = term.EnumTestField
            elif term.fieldName in [f.name for f in self.__efs__]:
                term.spaceIndex_= spacesNames[term.fieldName]
                term.numberingIndex_= numberingNames[term.fieldName]
                term.valuesIndex_= valuesNames[term.fieldName]
                term.internalType = term.EnumExtraField
            elif term.fieldName in [f.name for f in self.__ipefs__]:
                term.valuesIndex_= [ef.name for ef in  self.__ipefs__ ].index(term.fieldName)
                term.internalType = term.EnumExtraIPField
            else :
                term.internalType = term.EnumError
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

    def SetOnlyEvaluation(self,onlyEvaluation= True):
      """
      To activate the Only Evaluation functionality
          For the evaluation we only add the constribution without doing the integration (multiplication by the detjac )
          the user is responsible of dividing by the mass matrix to get the correct values
          . Ffr example  the user can use a discontinues field to generate element surface stress
      """
      self.onlyEvaluation = onlyEvaluation

    def ActivateElementType(self,domain):
        """
        Function to prepared the integration for a type of element
        domain : (ElementsContainer)

        """

        self.localNumbering = []
        for numbering in self.__usedNumbering__:
            self.localNumbering.append( numbering.get(domain.elementType,None) )


        self.p, self.w = self.integrationRule[domain.elementType];

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

        self.__usedValuesAtIP__ = [ipef.data[domain.elementType] for ipef in self.__ipefs__ ]

    def Integrate(self,wform,idstotreat):
        """
        Main function to execute the integration
        wform: (PyWeakForm) Python Or C++ version of the weak form to be integrated
        idstotreat:  list like (int) ids of the element to treat
        """
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
                        factor *= self.w[ip]

                    hasright = False

                    for term in monom:

                        if term.internalType == term.EnumNormal :
                            factor *= normal[term.derDegree]
                            continue
                        elif  term.internalType == term.EnumConstant :
                            factor *= constantsNumerical[term.valuesIndex_]
                            continue
                        elif  term.internalType == term.EnumUnknownField :
                            if term.derDegree == 1:
                                right = BxByBzI[term.spaceIndex_][[term.derCoordIndex_],:]
                            else:
                                right = np.array([NxNyNzI[term.spaceIndex_],])
                            rightNumbering = self.localNumbering[term.numberingIndex_][n,:] + self.unkownDofsOffset[term.valuesIndex_]
                            hasright = True
                            l2 = self.NumberOfShapeFunctionForEachSpace[term.spaceIndex_]
                            continue
                        elif  term.internalType == term.EnumTestField :
                            if term.derDegree == 1:
                                left = BxByBzI[term.spaceIndex_][term.derCoordIndex_]
                            else:
                                left = NxNyNzI[term.spaceIndex_]
                            leftNumbering = self.localNumbering[term.numberingIndex_][n,:] + self.testDofsOffset[term.valuesIndex_]
                            l1 = self.NumberOfShapeFunctionForEachSpace[term.spaceIndex_]
                            continue
                        elif term.internalType == term.EnumExtraField :

                            if term.derDegree == 1:
                                func = BxByBzI[term.spaceIndex_][term.derCoordIndex_]
                            else:
                                func = NxNyNzI[term.spaceIndex_]
                            centerNumbering = self.localNumbering[term.numberingIndex_][n,:]
                            vals = self.__usedValues__[term.valuesIndex_][centerNumbering]
                            factor *= np.dot(func,vals)

                            continue
                        elif term.internalType == term.EnumExtraIPField :
                            if term.derDegree == 1:
                                raise(Exception("Integration point field cant be derivated"))
                            val = self.__usedValuesAtIP__[term.valuesIndex_][n,ip]
                            factor *= val
                        else :
                            raise(Exception("Cant treat term " + str(term.fieldName)))

                    if factor == 0:
                        continue



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
        """
        Return number of non zero values in the vectors vK,iK,jK
        """
        return self.totalvijcpt

    def GetTotalTestDofs(self):
        """
        Return the number of dofs of the test space (number of rows of the K matrix)
        """
        return self.totalTestDofs

    def GetTotalUnkownDofs(self):
        """
        Return the number of dofs of the Unkown space (number of cols of the K matrix)
        """
        return self.totalUnkownDofs



def CheckIntegrity():
    import BasicTools.FE.Integration as Integration
    backup  = Integration.UseCpp

    Integration.UseCpp = False
    res = "ok"
    try:
        from BasicTools.FE.UnstructuredFeaSym import CheckIntegrity as CI
        res = CI()
    except:
        Integration.UseCpp = backup
        raise
    Integration.UseCpp = backup
    return res

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity())
