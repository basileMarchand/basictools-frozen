# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import sympy
from sympy import Symbol,Function
from sympy.matrices import Matrix

import numpy as np

testcharacter = "'"
space = Matrix([Symbol('x'),Symbol('y'), Symbol('z')])

def GetTestSufixChar():
    return testcharacter

def GetNormal(size):
    return GetField("Normal",size)

def GetConstant(name,size=1):
    if size == 1:
      return Symbol(name)
    else:
        res = []
        for i in range(size):
            res.append(Symbol(name+"_"+str(i)))
        return (Matrix([res])).T

def GetScalarField(alpha):

    if alpha is None:
        a = 1.
    elif isinstance(alpha,str):
        a = GetConstant(alpha)
    elif isinstance(alpha,(float,int)):
        a = float(alpha)
    else:
        a = alpha
    return a

def GetTestField(name,size,sdim=3,extraCoordinates=[]):
    return GetField(name,size,star=True,sdim=sdim,extraCoordinates=extraCoordinates)

def GetField(name,size,star=False,sdim=3,extraCoordinates=[]):
    res = []
    suffix = ""
    if star:
        suffix = testcharacter
    s = space[0:sdim]
    s.extend(extraCoordinates)

    if type(size) == int:
        if size == 1:
            if len(s) == 0:
                res.append(Function(name+suffix))
            else:
                res.append(Function(name+suffix)(*s))
        else:
            for i in range(size):
                res.append(Function(name+"_"+str(i)+suffix)(*s))
    else:
        res = [[None]*size[1] for x in range(size[0])]
        for i in range(size[0]):
            for j in range(size[1]):
                t = Function(name+"_"+str(i)+"_"+str(j)+suffix)(*s)
                res[i][j] = t
        return Matrix(res).T


    return (Matrix([res])).T

########################## Mathematical operators ##############################
def Diag(arg):
    shape = len( arg)
    res = [[0]*shape for i in range(shape)]
    for i,v in enumerate(arg):
        res[i][i] = v

    if type(arg).__module__ == np.__name__ or type(arg)==list:
        return np.array(res)
    else :
        return Matrix(res)

def Inner(a,b):
    return a@b

def Trace(arg):
    if type(arg).__module__ == np.__name__:
        return np.trace(arg)
    else:
        return sympy.trace(arg)

def Divergence(arg,sdim=3):
    return Trace(Gradient(arg,sdim=sdim) )

def Gradient(arg,sdim=3):
    shape = arg.shape[0]
    res = [[0]*shape for i in range(sdim)]
    for s in range(shape):
        for d in range(sdim):
            res[d][s] = arg[s].diff(space[d])

    if type(arg).__module__ == np.__name__:
       return np.array(res)
    else :
        return Matrix(res)

def Cross(a,b):
    shape = a.shape[0]
    res = [[0]*shape]

    if shape == 3:
        res[0][0] =   a[1]*b[2]-a[2]*b[1]
        res[0][1] = -(a[0]*b[2]-a[2]*b[0])
        res[0][2] =   a[0]*b[2]-a[1]*b[0]
    else:
        raise

    if type(a).__module__ == np.__name__:
       return np.array(res)
    else :
        return Matrix(res)

def Strain(arg ,sdim=3):
    G = Gradient(arg,sdim)
    return (G+G.T)/2

def StrainAxyCol(arg,radius):
    # strain definition for axisymmetric mechanical problem
    u = arg[0]
    v = arg[1]
    r = space[0]
    h = space[1]
    dudr = u.diff(r)
    dudh = u.diff(h)
    dvdr = v.diff(r)
    dvdh = v.diff(h)
    # E_r, E_z, E_theta, E_rz
    res = [ dudr, dvdh, u/radius, dudh+dvdr  ]
    if type(arg).__module__ == np.__name__:
        return np.array(res)
    else:
        return Matrix(res)


def ToVoigtEpsilon(arg):
    """ https://en.wikipedia.org/wiki/Voigt_notation
    """
    if arg.shape[0] ==3:
        res = [arg[0,0],arg[1,1],arg[2,2],2*arg[1,2],2*arg[0,2],2*arg[0,1], ]
    elif arg.shape[0] ==2:
        res = [arg[0,0],arg[1,1],2*arg[0,1]]
    elif arg.shape[0] ==1:
        res = [arg[0,0]]
    else:
        raise()

    if type(arg).__module__ == np.__name__:
        return np.array(res)
    else:
        return Matrix(res)

def ToVoigtSigma(arg):
    """ https://en.wikipedia.org/wiki/Voigt_notation
    """
    if arg.shape[0] ==3:
        res = [arg[0,0],arg[1,1],arg[2,2],arg[1,2],arg[0,2],arg[0,1], ]
    elif arg.shape[0] ==2:
        res = [arg[0,0],arg[1,1],arg[0,1]]
    elif arg.shape[0] ==1:
        res = [arg[0,0]]
    else:
        raise()

    if type(arg).__module__ == np.__name__:
        return np.array(res)
    else:
        return Matrix(res)

def FromVoigtSigma(arg):
    """ https://en.wikipedia.org/wiki/Voigt_notation
    """
    if arg.shape[0] == 6:
        res = [[arg[0], arg[5], arg[4] ],
               [arg[5], arg[1], arg[3] ],
               [arg[4], arg[3], arg[2] ],]
    elif arg.shape[0] == 3:
        res =  [[arg[0] ,arg[2]  ],
                [arg[2], arg[1] ]]
    elif arg.shape[0] ==1:
        res = [arg[0]]
    else:
        raise()

    if type(arg).__module__ == np.__name__:
        return np.array(res)
    else:
        return Matrix(res)


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
    pprint(u,use_unicode=GUI)
    pprint(Gradient(u),use_unicode=GUI)

    pprint(Strain(u),use_unicode=GUI)
    pprint(u[0].diff(space[1]),use_unicode=GUI)


    from BasicTools.FE.MaterialHelp import HookeIso
    K = HookeIso(1,0.3)
    pprint(K,use_unicode=GUI)

    ener = ToVoigtEpsilon(Strain(u+u0)).T*K*ToVoigtEpsilon(Strain(ut))+ f.T*ut*alpha
    pprint(ener,use_unicode=GUI)


    pprint(Gradient(u))
    pprint(Trace(Gradient(u)))
    print(type(Gradient(u)))
    print(type(u))
    print(type(u).__module__)
    print(type(sympy))
    print(sympy.__name__)

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))# pragma: no cover
