# -*- coding: utf-8 -*-

import numpy as np
from sympy.matrices import Matrix
from sympy import Symbol,Function
from sympy import pprint

space = Matrix([Symbol('x'),Symbol('y'), Symbol("z")])
testcharacter = "'"

UseCpp = False
try:
    from BasicTools.FE.WeakFormNumerical import PyWeakTerm as PyWeakTerm
    from BasicTools.FE.WeakFormNumerical import PyWeakMonom as PyWeakMonom
    from BasicTools.FE.WeakFormNumerical import PyWeakForm as PyWeakForm
    UseCpp = True
except:
    UseCpp = False
    print("Waring WeakFormNumerical (cpp) not avilable, using python variant")

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

    class PyWeakTerm(object):
        def __init__(self):
            super(PyWeakTerm,self).__init__()
            self.fieldName = ""
            self.derCoordName = ""
            self.derCoordIndex = 0
            self.derDegree = -1
            self.constant = False
            self.normal = False
        def __str__(self):
            res = ""
            if self.derDegree > 0 and self.normal == 0 :
    #            #res += "d" + self.fieldName + "/"  + "d"  + str(self.derCoordName)
                res += "Derivative("+str(self.fieldName)+"(x, y, z), "+str(self.derCoordName)+")"
            else:
                res += self.fieldName
            return res


def GetNormal(size):
    return GetField("Normal",size)

def GetConstant(name):
    return Symbol(name)

def GetTestField(name,size,sdim=3,extraCoordinates=[]):
    return GetField(name,size,star=True,sdim=sdim,extraCoordinates=extraCoordinates)

def GetField(name,size,star=False,sdim=3,extraCoordinates=[]):
    res = []
    suffix = ""
    if star:
        suffix = testcharacter
    s = space[0:sdim]
    s.extend(extraCoordinates)

    if size == 1:
        if sdim == 0:
            res.append(Function(name+suffix))
        else:
            res.append(Function(name+suffix)(*s))
    else:
        for i in range(size):
            res.append(Function(name+"_"+str(i)+suffix)(*s))
    return (Matrix([res])).T

def Divergence(arg,dim=3):
    res = 0
    for i in range(dim):
        res += arg.diff(space[i])

def Gradient(arg,sdim=3):
    shape = arg.shape[0]
    res = [[0]*shape for i in range(sdim)]

    for s in range(shape):
        for d in range(sdim):
            res[s][d] = arg[s].diff(space[d])
    return Matrix(res)

def Strain(arg ,sdim=3):
    G = Gradient(arg,sdim)
    return (G+G.T)/2

def ToVoigtEpsilon(arg):
    """ we use yamma for shear

    """
    if arg.shape[0] ==3:
        return Matrix([arg[0,0],arg[1,1],arg[2,2],2*arg[1,2],2*arg[0,2],2*arg[0,1], ])
    if arg.shape[0] ==2:
        return Matrix([arg[0,0],arg[1,1],2*arg[0,1]])
    if arg.shape[0] ==1:
        return Matrix([arg[0,0]])
    raise()

def ToVoigtSigma(arg):
    if arg.shape[0] ==3:
        return Matrix([arg[0,0],arg[1,1],arg[2,2],arg[1,2],arg[0,2],arg[0,1], ])
    if arg.shape[0] ==2:
        return Matrix([arg[0,0],arg[1,1],arg[0,1]])
    if arg.shape[0] ==1:
        return Matrix([arg[0,0]])
    raise()

def GetMecaElasticProblem(name="u",dim=3,K=None,planeStress=True):
    u = GetField("u",dim)
    ut = GetTestField("u",dim)
    if K is None:
        from MaterialHelp import HookeIso
        K = HookeIso(1,0.3,dim, planeStress)
    ener = ToVoigtEpsilon(Strain(u)).T*K*ToVoigtEpsilon(Strain(ut))
    return ener

def GetMecaNormalPressure(flux="p",name="u", dim=3):
    ut = GetTestField("u",dim)
    if isinstance(flux,str):
        p = GetConstant(flux)
    else:
        p = float(flux)

    from BasicTools.FE.WeakForm import GetNormal
    Normal = GetNormal(dim)

    wflux = p*Normal.T*ut
    return wflux




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

        pprint(exp)
        raise ()

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
        t.derDegree = 1
        t.fieldName = str(arg.args[0].func)
        t.derCoordName = str(arg.args[1])
        sn = [ str(c) for c in space]
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
    from sympy import pprint
    #init_session()

    print(space)

    u = GetField("u",3)
    u0 = GetField("u0",3)

    ut = GetTestField("u",3)
    f = GetField("f",3)
    alpha = Symbol("alpha")

    globalconstant = GetField("g",1,sdim=0)
    print(globalconstant )


    print(u)
    print(u.shape)
    print(u.diff(Symbol("x")))
    print(ut.diff(Symbol("x")))
    print("-----------------")
    pprint(u)
    pprint(Gradient(u))

    pprint(Strain(u))
    pprint(u[0].diff(space[1]))


    from BasicTools.FE.MaterialHelp import HookeIso
    K = HookeIso(1,0.3)
    pprint(K)

    ener = ToVoigtEpsilon(Strain(u+u0)).T*K*ToVoigtEpsilon(Strain(ut))+ f.T*ut*alpha
    pprint(ener)

    wf = SymWeakToNumWeak(ener)
    print("coucou")
    print([str(wf.GetMonom(i)) for i in range(wf.GetNumberOfTerms())])
    print("coucouII")

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
    print(CheckIntegrity())# pragma: no cover
