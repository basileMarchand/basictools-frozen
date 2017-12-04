# -*- coding: utf-8 -*-
from BasicTools.FE.Spaces.SymSpace import SymSpaceBase
import BasicTools.FE.ElementNames as EN
from sympy.matrices import Matrix
import numpy as np

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


def plot2DSquare(Space):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm

    # Create triangulation.
    #x = np.arange(0, 1.1, 0.1)
    #y = np.arange(0, 1.1, 0.1)
    #X, Y = np.meshgrid(x,y)
    #
    #print(Y)
    #z = X*.0;
    #np.empty((len(x),len(y)),dtype=np.float)

    ep = np.array([[ 0, 0],[ 0.5, 0],[ 1, 1.2],[ 0, 0.5]])
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
                #Z[i,j] = Space.GetShapeFunc([xi, eta])[sf]
                #Z[i,j] = Space.GetShapeFuncDer([xi, eta])[0,sf]
                #Z[i,j] = Space.GetShapeFuncDer([xi, eta])[1,sf]
                #Z[i,j] = Space.Eval_FieldI(0,ep[:,0],Jack,Jinv,-1)
                #Z[i,j] = Space.Eval_FieldI(0,ep[:,1],Jack,Jinv,-1)
                Z[i,j] = Space.Eval_FieldI(0,ep[:,0],Jack,Jinv,0)


        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        #surf = ax.plot_surface(X, Y, z, cmap=cm.coolwarm,
        #               linewidth=1, antialiased=False)
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,cmap=cm.jet)


        plt.title('Quad grid')

        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
        plt.pause(1)

def CheckIntegrity(GUI=False):
    return "ok"
