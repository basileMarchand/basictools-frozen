# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       
import numpy as np
from sympy.matrices import Matrix
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase

class Quad_P0_Global(SymSpaceBase):
    def __init__(self):
        super(Quad_P0_Global,self).__init__()
        self.geoSupport = EN.GeoQuad


        self.symN = Matrix([1])
        self.posN = np.array([ [ None,None] ])
        self.dofAttachments = [("G",None,None)]
        self.Create()

class Quad_P0_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Quad_P0_Lagrange,self).__init__()
        self.geoSupport = EN.GeoQuad


        self.symN = Matrix([1])
        self.posN = np.array([ [ 0.5, 0.5] ])
        self.dofAttachments = [("C",0,None)]
        self.Create()



class Quad_P1_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Quad_P1_Lagrange,self).__init__()
        self.geoSupport = EN.GeoQuad


        xi = self.xi
        eta = self.eta

        self.symN = Matrix([(1-xi)*(1-eta),
                            ( +xi)*(1-eta),
                            ( +xi)*( +eta),
                            (1-xi)*( +eta)])
        self.posN = np.array([[ 0, 0],
                              [ 1, 0],
                              [ 1, 1],
                              [ 0, 1]])
        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ]
        self.Create()

class Quad8_P2_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Quad8_P2_Lagrange,self).__init__()
        self.geoSupport = EN.GeoQuad

        xi = self.xi
        x1 = 1-xi
        x2 = xi

        eta = self.eta
        e1 = 1-eta
        e2 = eta

        XI = 2*xi-1
        ETA = 2*eta-1

        self.symN = Matrix([x1*e1*(+1-2*xi-2*eta), # linear part # (1-XI)*(1-ETA)*(-1-XI-ETA)/4
                            x2*e1*(-1+2*xi-2*eta),               # (1+XI)*(1-ETA)*(-1+XI-ETA)/4
                            x2*e2*(-3+2*xi+2*eta),               # ...
                            x1*e2*(-1-2*xi+2*eta),

                            (1-ETA   )*(1-XI**2)/2, #edges of base
                            (1-ETA**2)*(1+XI)/2,
                            (1+ETA   )*(1-XI**2)/2,
                            (1-ETA**2)*(1-XI)/2,

                            ])
        self.posN = np.array([[ 0, 0],
                              [ 1, 0],
                              [ 1, 1],
                              [ 0, 1],
                              [ 0.5, 0 ],
                              [ 1.0, 0.5],
                              [ 0.5, 1.],
                              [ 0, 0.5],
                              ])
        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ("F",0,None),
                               ("F",1,None),
                               ("F",2,None),
                               ("F",3,None),
                               ]
        self.Create()


class Quad_P2_Lagrange(SymSpaceBase):
    def __init__(self):
        super(Quad_P2_Lagrange,self).__init__()
        self.geoSupport = EN.GeoQuad



        xi = self.xi
        x1 = 1-xi
        x2 = xi

        eta = self.eta
        e1 = 1-eta
        e2 = eta

        xA = x1*(2*x1-1)
        xB = 4*x1*x2
        xC = x2*(2*x2-1)

        eA = e1*(2*e1-1)
        eB = 4*e1*e2
        eC = e2*(2*e2-1)



        self.symN = Matrix([xA*eA, # linear part
                            xC*eA,
                            xC*eC,
                            xA*eC,

                            xB*eA, #edges of base
                            xC*eB,
                            xB*eC,
                            xA*eB,

                            xB*eB, # center
                            ])
        self.posN = np.array([[ 0, 0],
                              [ 1, 0],
                              [ 1, 1],
                              [ 0, 1],
                              [ 0.5, 0 ],
                              [ 1.0, 0.5],
                              [ 0.5, 1.],
                              [ 0, 0.5],
                              [ 0.5, 0.5],
                              ])
        self.dofAttachments = [("P",0,None),
                               ("P",1,None),
                               ("P",2,None),
                               ("P",3,None),
                               ("F",0,None),
                               ("F",1,None),
                               ("F",2,None),
                               ("F",3,None),
                               ("C",0,None),

                               ]
        self.Create()

def plot2DSquare(Space):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm


    ep = np.array([[ 0, 0],[ 0.5, 0],[ 1, 1.2],[ 0, 0.5]])
    ep = Space.posN
    #ep = Space.posN
    X = np.empty((11,11),dtype=np.float)
    Y = np.empty((11,11),dtype=np.float)
    Z = np.empty((11,11),dtype=np.float)

    for xi in range(11):
        xp = xi/10.
        for yi in range(11):
            yp = yi/10.
            l1 = xp*ep[1,:]+(1-xp)*ep[0,:]
            l2 = xp*ep[2,:]+(1-xp)*ep[3,:]
            p  = (yp)*l2 +(1-yp)*l1

            X[xi,yi] =  p[0]
            Y[xi,yi] =  p[1]

    from mpl_toolkits.mplot3d import Axes3D
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = None
    for sf in range(len(Space.posN)):
        #cpt =0
        #for cpt in range(len(X)):
        for i in range(11):
            for j in range(11):
                p = [X[i,j],Y[i,j]]
                xi = i/10.
                eta = j/10.
                Space.SetIntegrationRule([[xi, eta]],[1])
                Jack, Jdet, Jinv = Space.GetJackAndDetI(0,ep)


                #Z[i,j] = Space.valN[0][sf]
                Z[i,j] = Space.GetShapeFunc([xi, eta])[sf]
                #Z[i,j] = Space.GetShapeFuncDer([xi, eta])[0,sf]
                #Z[i,j] = Space.GetShapeFuncDer([xi, eta])[1,sf]
                #Z[i,j] = Space.Eval_FieldI(0,ep[:,0],Jack,Jinv,-1)
                #Z[i,j] = Space.Eval_FieldI(0,ep[:,1],Jack,Jinv,-1)
                #Z[i,j] = Space.Eval_FieldI(0,ep[:,0],Jack,Jinv,0)




        #surf = ax.plot_surface(X, Y, z, cmap=cm.coolwarm,
        #               linewidth=1, antialiased=False)

        if surf is not None:
            surf.remove()
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,cmap=cm.jet)

        #fig.colorbar(surf, shrink=0.5, aspect=5)


        plt.title('Quad grid')


        plt.draw()
        plt.pause(1)


def CheckIntegrity(GUI=False):
    if GUI:
       plot2DSquare(Quad_P1_Lagrange())
       plot2DSquare(Quad_P2_Lagrange())
       plot2DSquare(Quad8_P2_Lagrange())
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
