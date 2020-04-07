# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np
from sympy import pprint
from sympy.core.containers import Tuple
from sympy import Symbol,Function

from BasicTools.FE.SymWeakForm import GetNormal


UseCpp = False
try:
    from BasicTools.FE.WeakForms.NativeNumericalWeakForm import PyWeakTerm as PyWeakTerm
    from BasicTools.FE.WeakForms.NativeNumericalWeakForm import PyWeakMonom as PyWeakMonom
    from BasicTools.FE.WeakForms.NativeNumericalWeakForm import PyWeakForm as PyWeakForm
    UseCpp = True
except:
    UseCpp = False
    print("Waring NativeNumericalWeakForm (cpp) not avilable, using python variant")

    class PyWeakForm(object):
        def __init__(self):
            super(PyWeakForm,self).__init__()
            self.form = []
        def AddTerm(self,monom):
            self.form.append(monom)
        def GetNumberOfTerms(self):
            return len(self.form)
        def GetMonom(self,i):
            return self.form[i]
        def GetRightPart(self,unknownvars):
            res = PyWeakForm()
            for p in self:
                for uv in unknownvars:
                    if p.hasVariable(uv):
                        break
                else:
                    res.AddTerm(p)
            return res

        def  GetLeftPart(self,unknownvars):
            res = PyWeakForm()
            for p in self:
                tocopy = False
                for uv in unknownvars:
                   if p.hasVariable(uv):
                        tocopy =True
                        break
                if tocopy:
                    res.AddTerm(p)
            return res

        def __str__(self):
            res = ""
            for i in range(self.GetNumberOfTerms()):
                res +=str(self.GetMonom(i))   + "\n"
            return res

        def __iter__(self):
            return iter(self.form)

    class PyWeakMonom(object):
        def __init__(self):
            super(PyWeakMonom,self).__init__()
            self.prefactor = 1
            self.prod = []
        def AddProd(self,term):
            self.prod.append(term)

        def GetNumberOfProds(self):
            return len(self.prod)

        def GetProd(self, n):
            return self.prod[n]

        def  hasVariable(self,var):
            for m in self :
                if m.fieldName == str(var):
                    return True
            return False


        def __str__(self):
            res = str(self.prefactor)
            for i in range(self.GetNumberOfProds()):
                res += "*"
                res += str(self.GetProd(i))
            return res

        def __iter__(self):
            return iter(self.prod)

        def copy(self):
            import copy
            return copy.deepcopy(self)

    class PyWeakTerm(object):
        EnumError = -1
        EnumNormal = 0
        EnumConstant = 1
        EnumUnknownField = 2
        EnumTestField = 3
        EnumExtraField = 4
        EnumExtraIPField = 5

        def __init__(self):
            super(PyWeakTerm,self).__init__()
            self.fieldName = ""
            self.derCoordName = ""
            self.derDegree = -1
            self.constant = False
            self.normal = False

            # internal data repr for integration
            self.internalType = PyWeakTerm.EnumError
            self.spaceIndex_ = None
            self.numberingIndex_ = None
            self.valuesIndex_ = None
            self.modeIndex_ = None
            self.derCoordIndex_ = None

        def __str__(self):
            res = ""
            if self.derDegree > 0 and self.normal == 0 :
    #            #res += "d" + self.fieldName + "/"  + "d"  + str(self.derCoordName)
                res += "Derivative("+str(self.fieldName)+", "+str(self.derCoordName)+","+str(self.derDegree)+")"
            else:
                res += self.fieldName
            return res

        def copy(self):
            import copy
            return copy.deepcopy(self)

def SymWeakMonomToNumWeakMono(exp):
    from  sympy.core.mul import Mul
    from  sympy.core.power import Pow

    if exp.func == Mul:
        res = PyWeakMonom()
        for arg in exp.args:
            if arg.is_Number:
                res.prefactor = float(arg)
                continue

            if isinstance(arg,Pow):
                term = ConverTermToProd(arg.args[0])
                if term is None:
                    print(type(arg.args[0]))
                    print(arg.args[0])
                    raise( Exception("Unable to treat term " + str(arg.args[0]) ))
                for i in range(arg.args[1]):
                    res.AddProd(term)
                continue

            term = ConverTermToProd(arg)
            if term is not None:
                #res.prod.append(term)
                res.AddProd(term)
                continue

            print(type(arg))
            print(arg)

            raise
        return res
    else:
        #only one term no product
        res = PyWeakMonom()
        term = ConverTermToProd(exp)
        if term is not None:
            res.AddProd(term)
        return res

def ConverTermToProd(arg):
    if isinstance(arg,Symbol):
        t = PyWeakTerm()
        t.constant = True
        t.derDegree = 0
        t.fieldName = str(arg)
        return t

    from sympy.core.function import Derivative
    if type(arg) == Derivative:
        t = PyWeakTerm()

        t.fieldName = str(arg.args[0].func)
        #Python 3
        t.derDegree = 1
        if Tuple == type(arg.args[1]):
            #sympy 1.2
            t.derCoordName = str(arg.args[1][0])
        else:
            #sympy 1.1.1
            t.derCoordName = str(arg.args[1])

        sn = []
        for i in range(0,len(arg.args[0].args)):
            sn.append(str(arg.args[0].args[i] ) )

        t.derCoordIndex_ =  sn.index(t.derCoordName)
        return t

    if isinstance(arg,Function):
        t = PyWeakTerm()
        N = GetNormal(3)
        #print(str(arg.func)+"* * * * * * ")
        #print(arg.func)
        #print (N)
        #print ([arg == nc for nc in N])
        if np.any([arg == nc for nc in N]):
            t.normal = True
            t.derDegree = int(str(arg.func).split("_")[1])
        else:
            t.derDegree = 0

        t.fieldName = str(arg.func)
        return t

    raise

def SymWeakToNumWeak(exp):
    from  sympy.core.add import Add
    from  sympy.core.mul import Mul

    exp = exp.expand()
    res = PyWeakForm()
    try:
        if exp.shape[0] == 1 and exp.shape[1] == 1:
            exp = exp[0,0]
    except:
        pass

    if exp.func == Mul:
        #res.AddTermform.append(SymWeakMonomToNumWeakMono(exp))
        res.AddTerm(SymWeakMonomToNumWeakMono(exp))

    elif exp.func == Add:
        for monom in exp.args:
            #res.form.append(SymWeakMonomToNumWeakMono(monom))
            res.AddTerm(SymWeakMonomToNumWeakMono(monom))

    else:
        mono = PyWeakMonom()
        mono.AddProd(ConverTermToProd(exp))
        res.AddTerm(mono)
        #raise(Exception("error treating formulation term"))

    return res



def CheckIntegrity(GUI=False):

    from BasicTools.FE.MaterialHelp import HookeIso
    from BasicTools.FE.SymWeakForm import Symbol,GetField,GetTestField,Strain,ToVoigtEpsilon

    u = GetField("u",3)
    u0 = GetField("u0",3)
    ut = GetTestField("u",3)
    f = GetField("f",3)
    alpha = Symbol("alpha")

    K = HookeIso(1,0.3)
    pprint(K,use_unicode=GUI)

    ener = ToVoigtEpsilon(Strain(u+u0)).T*K*ToVoigtEpsilon(Strain(ut))+ f.T*ut*alpha
    pprint(ener,use_unicode=GUI)

    wf = SymWeakToNumWeak(ener)

    print([str(wf.GetMonom(i)) for i in range(wf.GetNumberOfTerms())])

    unknames = ["u_0", "u_1", "u_2"]

    rwf = wf.GetRightPart(unknames )
    print(rwf)

    lwf = wf.GetLeftPart(unknames)
    print(lwf)


    pointdataT = GetTestField("pointdata",1)
    J_prim = (1*pointdataT)[0]

    print(J_prim)
    numwform = SymWeakToNumWeak(J_prim)
    print(numwform)

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))# pragma: no cover
