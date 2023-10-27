# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Optional, Union
import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.NumpyDefs import ArrayLike,PBasicFloatType
from BasicTools.ImplicitGeometry.ImplicitGeometryFactory import RegisterClass

def dsin(x):
    r = np.mod(x+np.pi/2,2*np.pi)-np.pi
    return -(abs(r)-np.pi/2)

def dcos(x):
    return dsin(x+np.pi/2)

def smoothmin(a,b,k=1):
    h = np.maximum( k-np.abs(a-b), 0.0 )/k;
    return np.minimum( a, b ) - h*h*k*(1.0/4.0);

def smoothmax(a,b,k=1):
    return -smoothmin(-a,-b,k=k)


class ImplicitGeometryBase(BaseOutputObject):
    def __init__(self):
        super(ImplicitGeometryBase,self).__init__()
        self.insideOut = False

    def ApplyVector(self, support,cellCenter=False,dim=None):
       if type(support).__module__ == np.__name__:
           if  len(support.shape) == 2 and support.shape[1] ==2:
                  support = np.hstack((support,np.zeros((support.shape[0],1) ,dtype=support.dtype ) ) )
           d  = self.GetDistanceToPoint(support)
       else:
           if cellCenter:
              from BasicTools.Containers.MeshTools import GetElementsCenters
              pos = GetElementsCenters(support,dim=dim)
           else:
              pos = support.GetPosOfNodes()
              if pos.shape[1] ==2:
                  pos = np.hstack((pos,np.zeros((pos.shape[0],1) ,dtype=pos.dtype ) ) )
           d  = self.GetDistanceToPoint(pos)
       return d

    def GetGradientDistanceToPoint(self, pos:np.ndarray, dx:Optional[Union[PBasicFloatType,ArrayLike] ] = None) -> np.ndarray:
        """Compute the numerical gradient of the implicit geometry
        the child classes can overwrite it for a more efficient or exact result.

        Parameters
        ----------
        pos : np.ndarray
            the position to evaluate the gradient
        dx : Optional[PBasicFloatType,Arraylike], optional
            the step (per coordinate) to be used to compute the numerical gradient, by default dx is 1e-6

        Returns
        -------
        np.ndarray
            _description_
        """

        if dx is None:
            dx = [1.e-6]*3
        elif not hasattr(dx, "__len__"):
            dx = [dx]*3

        res = np.zeros_like(pos)
        for i in range(res.shape[1]):
            if dx[i] != 0:
                pos[:,i] -= dx[i]
                res[:,i] = -self.GetDistanceToPoint(pos)
                pos[:,i] += 2*dx[i]
                res[:,i] += self.GetDistanceToPoint(pos)
                pos[:,i] -= dx[i]
                res[:,i] /= dx[i]
        return res

    def ApplyInsideOut(self,res):
        if self.insideOut :
           return -res
        return res

    def __call__(self,support,cellCenter=False):
        return self.ApplyVector(support,cellCenter=cellCenter)

    def GetDistanceToPoint(self,pos):
        raise(Exception("Please overload this function"))
    #def __str__(self):
    #    return "+-+-"

####################### cache object ################################
class ImplicitGeometryCachedData(ImplicitGeometryBase):
    def __init__(self,internalImplicitGeometry):
        super(ImplicitGeometryCachedData,self).__init__()
        self.internalImplicitGeometry = internalImplicitGeometry
        self.cacheData = None
        self.cachedInputGetDistanceToPoint = None
        self.cachedInputApplyVector = None
        self.cachedInputcellCenter = False
        self.cachedInputcelldim = 0

    def GetDistanceToPoint(self,pos):
        if self.cacheData is None:
            self.PrintDebug("building Cache")
            self.cacheData = self.internalImplicitGeometry.GetDistanceToPoint(pos)
        elif pos is not self.cachedInputGetDistanceToPoint:
            self.PrintDebug("rebuilding Cache")
            self.cacheData = self.internalImplicitGeometry.GetDistanceToPoint(pos)
        else :
            self.PrintDebug("Using Cache")

        self.cachedInputGetDistanceToPoint = pos
        self.cachedInputApplyVector = None

        return self.cacheData


    def ApplyVector(self, support,cellCenter=False,dim=None):
       if self.cacheData is not None:
           if support is self.cachedInputApplyVector and cellCenter == self.cachedInputcellCenter and dim == self.cachedInputcelldim:
               return self.cacheData

       res = super(ImplicitGeometryCachedData,self).ApplyVector(support,cellCenter=cellCenter,dim=dim)
       self.cachedInputApplyVector = support
       self.cachedInputcellCenter  = cellCenter
       self.cachedInputcelldim = dim

       return  res
####################### cache object ################################
class ImplicitGeometryDelayedInit (ImplicitGeometryBase):
    def __init__(self,name,ops={}):
        super(ImplicitGeometryDelayedInit,self).__init__()
        self.internalImplicitGeometry = None
        self.__name = name
        self.__ops = ops

    def Init(self):
        if self.internalImplicitGeometry is None:
          from BasicTools.ImplicitGeometry.ImplicitGeometryFactory import Create
          self.internalImplicitGeometry = Create(self.__name,self.__ops)
          ## release memory
          self.__name = None
          self.__ops = None
          if self.internalImplicitGeometry is None:
             raise (Exception('Error in the initialisation of  : ' + str(self.__name)) ) # pragma: no cover

    def ApplyVector(self, support,cellCenter=False,dim=None):
        self.Init()
        return self.internalImplicitGeometry.ApplyVector(support,cellCenter=cellCenter)

    def GetDistanceToPoint(self,pos):
        self.Init()
        return self.internalImplicitGeometry.GetDistanceToPoint(pos)

    def __str__(self):
        res = "ImplicitGeometryDelayedInit :\n"
        if self.__name is None:
            res += self.internalImplicitGeometry.__str__()
        else:
            res = str(self.__name) + "\n"
            res += str(self.__ops) + "\n"
        return res


def TryToCreate(name,ops=None):
    from BasicTools.ImplicitGeometry.ImplicitGeometryFactory import Create
    res = Create(name,ops)
    if res is None:
        raise Exception("Unable to Create a " + str(name))# pragma: no cover
    return res
#-----------------------------.

def CheckIntegrity(GUI=False):

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh3D = ConstantRectilinearMesh()
    myMesh3D.SetDimensions([30,30,30]);
    myMesh3D.SetSpacing(1./(myMesh3D.GetDimensions()-1)*2);
    myMesh3D.SetOrigin([-1.,-1.,-1.]);

    myMesh2D = ConstantRectilinearMesh(dim=2)
    myMesh2D.SetDimensions([30,30]);
    myMesh2D.SetSpacing(1./(myMesh2D.GetDimensions()-1)*2);
    myMesh2D.SetOrigin([-1.,-1.]);


    class DummyImplicitGeometry(ImplicitGeometryBase):
        def __init__(self):
            super(DummyImplicitGeometry,self).__init__()

        def GetDistanceToPoint(self, pos):
            return self.ApplyInsideOut(np.zeros(pos.shape[0]))

    RegisterClass("Dummy",DummyImplicitGeometry,withError=False)

    TryToCreate("Dummy")(myMesh3D)
    TryToCreate("Dummy")(myMesh2D)

    TryToCreate("Dummy").ApplyVector(myMesh3D,cellCenter=True)

    TwoPoints3D = np.array([[0.,0.,0.],[1.,2.,3.]], dtype=float)
    TryToCreate("Dummy")(TwoPoints3D)

    try:
        obj = ImplicitGeometryBase()
        obj.GetDistanceToPoint(TwoPoints3D)
        raise # pragma: no cover
    except:
        pass

    TwoPoints2D = np.array([[0.,0.],[1.,2.]], dtype=float)

    obj = TryToCreate("Dummy")
    obj.insideOut = True
    obj(TwoPoints2D)

    ####################################################################

    myObj6 = ImplicitGeometryCachedData(obj)
    myObj6(myMesh3D)
    myObj6(myMesh3D)
    myObj6.GetDistanceToPoint(TwoPoints3D)
    myObj6.GetDistanceToPoint(TwoPoints3D)

    myObj6.GetGradientDistanceToPoint(TwoPoints3D)
    myObj6.GetGradientDistanceToPoint(TwoPoints3D,1)
    #######################################################################

    res = ImplicitGeometryDelayedInit("Dummy")
    print(res)
    res(myMesh3D)
    res.GetDistanceToPoint(TwoPoints3D)
    print(res)

    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=False))

