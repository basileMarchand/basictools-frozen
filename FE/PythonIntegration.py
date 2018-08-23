# -*- coding: utf-8 -*-

# distutils: language = c++



import numpy as np


from scipy.sparse import coo_matrix
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
import BasicTools.FE.ElementNames as EN

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


#      self.mesh = mesh
      self.SetPoints(mesh.nodes)

      #self.internalfBiform = np.array((len(self.__ufs__),2,len(self.__tfs__),2),dtype=WeakForm)
      #self.internalfLiform = np.array((len(self.__tfs__),2),dtype=WeakForm)

#      cdef WeakMonom tt = WeakMonom()
#      cdef int ud
#      cdef int td
#      for monom in wform.form:
#        u = None
#        ud = 0
#        t = None
#        td = 0
#
#        tt.prefactor = monom.prefactor
#        for term in monom.prod:
#            if term.internalType == 0:
#                tt.prod.append(term)
#            elif term.internalType == 1:
#                tt.prod.append(term)
#            elif term.internalType == 2:
#                u = term.valuesIndex_
#                ud = term.derDegree
#            elif term.internalType == 3:
#                t = term.valuesIndex_
#                td = term.derDegree
#            elif term.internalType == 4:
#                tt.prod.append(term)
#            else:
#                raise
#
#        if u is None:
#            print(self.internalfLiform.shape)
#            self.internalfLiform[t,td].form.append(tt)
#        else:
#            self.internalfBiform[u,ud,t,td].form.append(tt)


  def SetPoints(self,nodes):
        """
        ## from https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC

        multiply (arr, value)

        Takes a numpy arry as input, and multiplies each elemetn by value, in place

        param: array -- a 2-d numpy array of np.float64

        """
        self.nodes = nodes

        return None


  #@initializedcheck(False)
  #@overflowcheck (False)

  def SetOnlyEvaluation(self,onlyEvaluation):
      self.onlyEvaluation = onlyEvaluation

  def ActivateElementType(self,domain):

    self.localNumbering = []
    for numbering in self.__usedNumbering__:
        self.localNumbering.append( numbering.get(domain.elementType,None) )

    #print(domain.elementType)
    #print(EN.geoSupport[domain.elementType])
    #print([str(x) for x in self.integrationRule.keys()])
    #print(self.integrationRule[EN.geoSupport[domain.elementType]])
    self.p, self.w = self.integrationRule[EN.geoSupport[domain.elementType]];
    #NumberOfIntegrationPoints = len(w_np)

    self.geoSpace = LagrangeSpaceGeo[domain.elementType]


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

    #cdef np.ndarray[DTYPEfloat_t, ndim=1]  F = np.zeros(self.totalTestDofs,dtype=DTYPEfloat)
    #cdef DTYPEfloat_t [::1] F = F_np;

#    #w_wp = np.zeros(0,dtype=np.float)
#    p, w_np = self.integrationRule[EN.geoSupport[domain.elementType]];
#    NumberOfIntegrationPoints = len(w_np)
#    #print("w_np")
#    #print(w_np)


#    cdef DTYPEfloat_t [:] w = w_np;


#    spaces = []
#    cdef np.ndarray[DTYPEint_t, ndim=1] NumberOfShapeFunctionForEachSpace = np.zeros(len(self.__usedSpaces__),dtype=DTYPEint)
#    cpt = 0
#    for space in self.__usedSpaces__:
#        spaces.append(space[domain.elementType])
#        space[domain.elementType].SetIntegrationRule(p,w)
#        NumberOfShapeFunctionForEachSpace[cpt] = space[domain.elementType].GetNumberOfShapeFunctions()
#        cpt += 1

    NumberOfIntegrationPoints = len(self.w)


    ev = np.empty(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms()*NumberOfIntegrationPoints,dtype=np.float)
    ei = np.empty(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms()*NumberOfIntegrationPoints,dtype=np.int)
    ej = np.empty(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms()*NumberOfIntegrationPoints,dtype=np.int)

    #print("----")
    #print(self.maxNumberOfElementVIJ*wform.GetNumberOfTerms())
    numberOfFields = len(self.__usedSpaces__)

    BxByBzI = [None] *numberOfFields

    NxNyNzI = [None] *numberOfFields



    for n in idstotreat:
        fillcpt =0
        #BaseOutputObject().PrintDebug("treating : " + str(n)  )

        # the coordinates of the nodes
        xcoor = self.nodes[self.connectivity[n],:]


        for ip in range(NumberOfIntegrationPoints):
            """ we recover the jacobian matrix """
            Jack, Jdet, Jinv = self.geoSpace.GetJackAndDetI(ip,xcoor)
            #Jack, Jdet, Jinv = GetJackAndDet(self.geoSpace.valdphidxi[ip],xcoor,self.geoSpace.GetDimensionality())
#            print("geoSpace.valdphidxi["+str(ip)+"]")
#            print(self.geoSpace.valdphidxi[ip]);
#
#            print("xcoor")
#            print(xcoor)


            for i in range(numberOfFields):
                 if self.localSpaces[i] is not None:
                     NxNyNzI[i] = self.localSpaces[i].valN[ip]
#                     print("Jack")
#                     print(Jack)
#                     print("self.localSpaces[i].valdphidxi[ip]")
#                     print(self.localSpaces[i].valdphidxi[ip])
                     BxByBzI[i] = Jinv(self.localSpaces[i].valdphidxi[ip])
#                     print("BxByBzI["+str(i)+"]")
#                     print(BxByBzI[i])
#                     print("error")
#                     print(Jack.dot(BxByBzI[i]) )
#                     print(Jack.dot(BxByBzI[i]) - self.localSpaces[i].valdphidxi[ip])
#

#                     print("--------------")
#                     print("solution II")
#                     md = Jack
#                     yd = self.localSpaces[i].valdphidxi[ip]
#                     import scipy.linalg as la
#                     s = la.solve(md.T.dot(md),md.T.dot(yd) )
#                     print(s)
            if self.hasnormal:
                normal = self.geoSpace.GetNormal(Jack)

            for monom in wform:
#            for cptmono in range(NumberOfTerms):
#                monom = wform.form[cptmono]

                factor = monom.prefactor
#                print(" factor(A) "+str(factor))
                if self.onlyEvaluation :
                    # For the evaluation we only add the constribution without doing the integration
                    # the user is responsible of dividing by the mass matrix to get the correct values
                    # also the user can use a discontinues field to generate element surface stress (for example)
                    pass
                else:
                    # for the integration we multiply by the deteminant of the jac
                    factor *= Jdet
#                    print("Jdet " + str(Jdet))
#                    print("Jack " + str(Jack))

                hasright = False

#                print(" factor(B) "+str(factor))

                for term in monom:
#                for xterm in range( monom.prod.size()):
 #                   term = monom.prod[xterm]
#                    print(" factor("+str(term.internalType)+") " +str( factor));

                    if term.internalType == 0 :
#                        print("normal" + str(normal[term.derDegree]))
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
                            #left = BxByBzI[term.spaceIndex_][[term.derCoordIndex_]].T
                            left = BxByBzI[term.spaceIndex_][term.derCoordIndex_]
                        else:
                            left = NxNyNzI[term.spaceIndex_]
                        leftNumbering = self.localNumbering[term.numberingIndex_][n,:] + self.testDofsOffset[term.valuesIndex_]
                        l1 = self.NumberOfShapeFunctionForEachSpace[term.spaceIndex_]
                        #cdef int [:] leftNumbering = leftNumbering_np;
                        continue
                    elif  term.internalType == 4 :

                        if term.derDegree == 1:
                            func = BxByBzI[term.spaceIndex_][term.derCoordIndex_]
                        else:
                            func = NxNyNzI[term.spaceIndex_]
                        centerNumbering = self.localNumbering[term.numberingIndex_][n,:]
                        vals = self.__usedValues__[term.valuesIndex_][centerNumbering]
                        #print("centerNumbering");
                        #print(centerNumbering);
                        #print("vals");
                        #print(vals);
                        #print("func");
                        #print(func);
                        #s = 0.
                        #for i in range(NumberOfShapeFunctionForEachSpace[term.spaceIndex_]):
                        #    s += func[i]*vals[i]
                        #factor *= s
                        factor *= np.dot(func,vals)
                        #print("factor");
                        #print(factor);
                        #if factor:
                        #   raise(Exception("tata"))

                        continue
                    else :
                        #print(term.internalType)
                        #print(self.unkwowNames)
                        #print(self.testNames)
                        #print(self.extraFieldsNames)
                        #print([str(f) for f in  monom.prod])
                        raise(Exception("Cant treat term " + str(term.fieldName)))
                    #raise (Exception( "Error field '" + term.fieldName + "' not found in the available fields" ) )

                if factor == 0:
                    continue

                #print(factor)

                factor *= self.w[ip]

                #print(factor)
                #l1 = leftNumbering.size
                #print(monom)

                if hasright:
#                    print("right")
#                    print(right)
#                    print("rightNumbering")
#                    print(rightNumbering)
#                    print("leftNumbering")
#                    print(leftNumbering)
#                    print("left")
#                    print(left)
#                    #l2 = len(rightNumbering)
                    l = l1*l2

                    #print("fillcpt")
                    #print(fillcpt)
                    l2cpt = fillcpt
                    #print(l1)
                    #print(l2)
                    #print(l2cpt)
                    for i in range(l1):
                        for j in range(l2) :
                            ev[l2cpt] =  left[i]*right[0,j]*factor
                            l2cpt +=1


                    #print(left)
                    #print(right)
                    #vals = np.dot(left,right).ravel()
                    #vals *= factor
                    #raise


                    #print(vals.ravel())
                    #print(ev[fillcpt:fillcpt+l])
                    #for i in range(l):
                    #    ev[fillcpt+i] =  vals[i]
                    #ev[fillcpt:fillcpt+l] =  vals.ravel()


                    #ej[fillcpt:fillcpt+l] = np.tile(rightNumbering,len(leftNumbering))

                    l2cpt = fillcpt
                    for i in range(l1):
                        for j in range(l2) :
                          ej[l2cpt] = rightNumbering[j]
                          l2cpt += 1
                    #for i in leftNumbering:
                    #    ej[l2cpt:l2cpt+l2] = rightNumbering
                    #    l2cpt += l2
                    l2cpt = fillcpt
                    for j in range(l2) :
                         for i in range(l1):
                              ei[l2cpt] = leftNumbering[j]
                              l2cpt += 1
                    #ei[fillcpt:fillcpt+l] = np.repeat(leftNumbering,l2).ravel()
#                    for i in range(l):
#                        print("m({0},{1}) = {2}".format(ei[i+fillcpt],ej[i+fillcpt],ev[i+fillcpt] ))
                    fillcpt += l
                else:
                     #print("***")
                    #print(w[ip])
                    #print(factor)
                    #print(left)
                    #print((w[ip]*factor)*left)
                    #print(F[leftNumbering])
                    #vals = (factor*left).ravel()
                    for i in range(l1):
                      self.F[leftNumbering[i]] += left[i]*factor

#                    print(leftNumbering)
                    #print(left)
                    #print(factor)
                    #raise(Exception("tata"))
                    #print(left*factor)
                    #F[leftNumbering] += (factor*left).ravel()


        #if len (ev):
        if fillcpt:
            data = coo_matrix((ev[:fillcpt], (ei[:fillcpt],ej[:fillcpt])), shape=( self.totalTestDofs,self.totalUnkownDofs))
            #.tocsr().tocoo()
            data.sum_duplicates()
            #data =
            #data = coo_matrix((ev, (ei,ej)), shape=( totaldofs,totaldofs)).tocsr().tocoo()
#            print(data)
            start = self.totalvijcpt
            stop = start+len(data.data)

            self.vK[start:stop] = data.data
            self.iK[start:stop] = data.row
            self.jK[start:stop] = data.col
            #vij_v[start:stop] = data.data
            #vij_i[start:stop] = data.row
            #vij_j[start:stop] = data.col
            self.totalvijcpt += len(data.data)
#        raise(Exception("toto"))


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
