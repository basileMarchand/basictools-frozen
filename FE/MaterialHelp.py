# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

def HookeIso(E,nu, dim = 3, planeStress = True):
    hl= HookeLaw({"E":E,"nu":nu})
    return hl.HookeIso(dim=dim ,planeStress= planeStress )



def LaplaceOrtho(k1,k2,k3=1, dim = 3):
    if dim == 2:
        return np.array([[k1 , 0 ],
                         [0  , k2]]);

    res= np.array([[k1 , 0 ,  0  ],
                   [0  , k2,  0  ],
                   [0  ,0  ,  k3 ]]);
    return res


class HookeLaw() :
    def __init__(self,opt=None):

        if not opt is None:
            self.Read(opt)

    def HookeIso(self, dim = 3, planeStress = True):
        E = self.Get("E")
        nu = self.Get("nu")
        if dim == 1:
            return np.array([[E],])
        if dim == 2:
            if planeStress:
              return ((float(E)/(1.-nu**2))*
              np.array([[1 , nu, 0     ],
                        [nu, 1 , 0     ],
                        [0 , 0 , (1-nu)/2.]]));
            else:
              return  ((float(E)/((1+nu)*(1-2*nu)))*
              np.array([[1-nu, nu  , 0 ],
                        [nu  , 1-nu, 0 ],
                        [0   , 0   , 0.5-nu]]));

        res= ((float(E)/((1+nu)*(1-2*nu)))*
            np.array([[1-nu, nu  ,  nu  , 0      ,0      ,0 ],
                      [nu  , 1-nu,  nu  , 0      ,0      ,0 ],
                      [nu  , nu  ,  1-nu, 0      ,0      ,0 ],
                      [0   , 0   ,  0   , 0.5-nu ,0      ,0 ],
                      [0   , 0   ,  0   , 0      ,0.5-nu ,0 ],
                      [0   , 0   ,  0   , 0      ,0      ,0.5-nu]]));
        return res

    def Read(self,opt):

        if len(opt) != 2:
            raise(Exception("I need only 2 parameter to define a elastic material"))

        self.data = {}
        parser = {  "lambda":"l",# Premier coefficient de Lamé (?)
                  "G":"G",       # Module de cisaillement (G)
                  "E":"E",       # Module de Young (E)
                  "K":"K",       # Module d'élasticité isostatique (K)
                  "mu":"K",       # Module d'élasticité isostatique (K) (mu) deuxiemme coeff de Lamé
                  "nu":"nu",     # Coefficient de Poisson (?)
                  "M":"M"}       # Module d'onde de compression (M, P-wave modulus)

        for key,data in opt.items():
            if key in parser:
                self.data[parser[key]] = float(data)
            else:
                raise(Exception("dont know how to treat the key " + str(key) ) )

        # from  https://fr.wikipedia.org/wiki/Coefficient_de_Lam%C3%A9
        #                              lambda & K          E & G                K & lambda        K & G             lambda & nu              G & nu                   E & nu                   K & nu         K & E                    M & G
        self.formulas ={"K":[     "K","l+2.*g/3."        ,"G*(3*l+2*G)/(l+G)",                                    "l*(1+nu)/(3*nu)"       ,"2*G*(1+nu)/(3*(1-2*nu))","E/(3*(1-2*nu))"                                              ,"M - 4*G/3"],
                        "E":[     "E","G(3*l+2*G)/(l+G)" ,                    "9*K*(K-l)/(3*K-l)","9*K*G/(3*K-G)","l*(1+nu)*(1-2*nu)/(nu)","2*G*(1+nu)"                                      ,"3*K*(1-2*v)"                        ,"G*(3*M-4*G)/(M-G)"]    ,
                        "lambda":["l",                    "G*(E-2*G)/(3*G-E)",                    "K-2*G/3"                               ,"2*G*nu/(2*nu)"          ,"E*nu/((1+nu)*(1-2*nu))","3*K*nu/(1+nu)","3*K*(3*k-E)/(9*K-E)","M-2*G"],
                        "G":["G"],# to be completed
                        "nu":["nu"],# to be completed
                        "M":["M"]# to be completed
                         }
    def Get(self,name):
        formulas = self.formulas[name]
        for f in formulas :
            try:
                return eval(f,self.data)
            except:
                pass
        raise (Exception ("Unable to calculate " + str(name)))


def CheckIntegrity():
    HI3D = HookeIso(1,0.3)
    HI2D = HookeIso(1,0.3,dim=2)
    HI2DP = HookeIso(1,0.3,dim=2,planeStress= False)

    hl= HookeLaw({"E":1,"nu":0.3})
    print(hl.Get("K"))
    print(hl.Get("E"))
    print(hl.Get("nu"))
    print(hl.Get("lambda"))
    LO3DII = hl.HookeIso()
    diff = HI3D-LO3DII
    print(diff)

    HI2DII = hl.HookeIso(dim=2)
    print(HI2D-HI2DII)

    HI2DPII = hl.HookeIso(dim=2,planeStress= False )
    print(HI2DP-HI2DPII)

    LO3D = LaplaceOrtho(1,2,3)
    LO2D= LaplaceOrtho(1,2,dim=2)


    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover