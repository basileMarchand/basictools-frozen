# -*- coding: utf-8 -*-

import numpy as np
from sympy.matrices import Matrix
from sympy import Symbol,Function
from sympy import pprint

space = Matrix([Symbol('x'),Symbol('y'), Symbol("z")])
testcharacter = "'"

def GetNormal(size):
    return GetField("Normal",size)

def GetConstant(name):
    return Symbol(name)

def GetTestField(name,size,extraCoordinates=[]):
    return GetField(name,size,star=True,extraCoordinates=extraCoordinates)

def GetField(name,size,star=False,sdim=3,extraCoordinates=[]):
    res = []
    suffix = ""
    if star:
        suffix = testcharacter
    s = space[0:sdim]
    s.extend(extraCoordinates)
    if size == 1:
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



class Weakterm(object):
    def __init__(self):
        self.fieldName = ""
        self.derCoordName = ""
        self.derCoordIndex = 0
        self.derDegree = -1
        self.constant = False
        self.normal = False

    def __str__(self):
        res = ""
        if self.derDegree > 0 :
            #res += "d" + self.fieldName + "/"  + "d"  + str(self.derCoordName)
            res += "Derivative("+self.fieldName+"(x, y, z), "+self.derCoordName+")"
        else:
            res += self.fieldName
        return res

class WeakMonom(object):
    def __init__(self):
        self.prefactor = 1
        self.prod = []

    def hasTestFunc(self):
        for p in self.prod :
            if p.fieldName[-1] == testcharacter:
                return True

    def hasVariable(self,var):
        for p in self.prod :
            if p.fieldName == str(var):
                return True

    def __str__(self):
        res = str(self.prefactor)
        for p in self.prod:
            res += "*"
            res += str(p)

        return res

class WeakForm(object):
    def __init__(self):
        self.OnTag = ""
        self.form = []
    def GetLeftPart(self,unknownvars):
        import copy
        res = WeakForm()
        res.OnTag = self.OnTag
        for p in self.form:
            tocopy = False
            for uv in unknownvars:
                if p.hasVariable(uv):
                    tocopy =True
                    break
            if tocopy:
                res.form.append(copy.deepcopy(p))
        return res
    def GetRightPart(self,unknownvars):
        import copy
        res = WeakForm()
        res.OnTag = self.OnTag
        for p in self.form:
            for uv in unknownvars:
                if p.hasVariable(uv):
                    break
            else:
                res.form.append(copy.deepcopy(p))
        return res
    def __str__(self):
        res = "Weak Form:\n"
        res += "  Nb Terms : " + str(len(self.form)) +"\n"
        fields = set()
        const = set()
        coords = set()
        for p in self.form:
            for f in p.prod:
                if f.constant  :
                    const.add(f.fieldName)
                    continue
                fields.add(f.fieldName)
                if f.derDegree > 0:
                    coords.add(f.derCoordName)
        res += "  Constants : " + " ".join(const) + "\n"
        res += "  Fields : " + " ".join(fields) + "\n"
        res += "  Coordinates used : " + " ".join(coords) + "\n"
        return res

def SymWeakMonomToNumWeakMono(exp):
    from  sympy.core.mul import Mul
    from  sympy.core.power import Pow

    if exp.func == Mul:
        res = WeakMonom()
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
                    res.prod.append(term)
                continue

            term = ConverTermToProd(arg)
            if not term is None:
                res.prod.append(term)
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
        t = Weakterm()
        t.constant = True
        t.derDegree = 0
        t.fieldName = str(arg)
        return t

    from sympy.core.function import Derivative
    if type(arg) == Derivative:
        t = Weakterm()
        t.derDegree = 1
        t.fieldName = str(arg.args[0].func)
        t.derCoordName = str(arg.args[1])
        sn = [ str(c) for c in space]
        t.derCoordIndex =  sn.index(t.derCoordName)
        return t

    if isinstance(arg,Function):
        t = Weakterm()
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
    res = WeakForm()
    try:
        if exp.shape[0] == 1 and exp.shape[1] == 1:
            exp = exp[0,0]
    except:
        pass

    if exp.func == Mul:
        res.form.append(SymWeakMonomToNumWeakMono(exp))

    elif exp.func == Add:
        for monom in exp.args:
            res.form.append(SymWeakMonomToNumWeakMono(monom))
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
    ener = ToVoigtEpsilon(Strain(u+u0)).T*K*ToVoigtEpsilon(Strain(ut))+ f.T*ut*alpha
    pprint(ener)

    wf = SymWeakToNumWeak(ener)
    print([str(r) for r in wf.form])

    rwf = wf.GetRightPart(["u0", "u1", "u2"])
    print([str(r) for r in rwf.form])

    lwf = wf.GetLeftPart(["u0", "u1", "u2"])
    print([str(r) for r in lwf.form])

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
