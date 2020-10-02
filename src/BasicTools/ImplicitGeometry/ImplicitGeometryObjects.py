# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
import math

import BasicTools.Helpers.ParserHelper as PH

from BasicTools.ImplicitGeometry.ImplicitGeometryFactory import RegisterClass
from BasicTools.ImplicitGeometry.ImplicitGeometryBase import ImplicitGeometryBase,dsin,dcos
from BasicTools.ImplicitGeometry.ImplicitGeometryOperators import ImplicitGeometryIntersection,ImplicitGeometryUnion



class ImplicitGeometryWrapped(ImplicitGeometryBase):
    """Wrapper to add the ImplicitGeometry API arround a numpy array
        field : numpy array
    """

    def __init__(self, field=None ):
        super(ImplicitGeometryWrapped,self).__init__()
        self.field = field

    def ApplyVector(self, support,cellCenter=False):
        if type(support).__module__ == np.__name__:
            if support.shape[0] != len(self.field):
                raise(Exception("suport and field not compatible"))
        else:
            if cellCenter :
                if support.GetNumberOfElements() != len(self.field):
                    raise(Exception("suport and field not compatible"))
            elif support.GetNumberOfNodes() != len(self.field):
                raise(Exception("suport and field not compatible"))

        return self.ApplyInsideOut(self.field)

    def GetDistanceToPoint(self,_pos):
        raise Exception("Cant use this function") # pragma: no cover

RegisterClass("Wrapped",ImplicitGeometryWrapped)

def CreateImplicitGeometryExternalSurface(ops):
     """
     ImplicitGeometry of the external surface of a mesh (support)
     support : Grid to compute the surface

     <All id="x" support="" >
     """
     res = ImplicitGeometryExternalSurface(ops["support"])
     if "support" in ops:
         res.SetSupport(ops["support"])
     return res

class ImplicitGeometryExternalSurface(ImplicitGeometryBase):
    """ImplicitGeometry of the external surface of a mesh (support)
        support : mesh to compute the surface
    """
    def __init__(self,support = None):
        super(ImplicitGeometryExternalSurface,self).__init__()
        self.internalImplicitGeometry = None
        self.offset = 1E-3
        if support is not None:
            self.SetSupport(support)

    def SetSupport(self, support):
      self.internalImplicitGeometry = ImplicitGeometryStl()
      self.internalImplicitGeometry.SetMesh(support)
      self.internalImplicitGeometry.ismanifold = True
      support.ComputeBoundingBox()
      self.offset = 1.*np.linalg.norm(support.boundingMax-support.boundingMin)/1000


    def GetDistanceToPoint(self,pos):
        return self.internalImplicitGeometry.GetDistanceToPoint(pos) + self.offset

RegisterClass("All",ImplicitGeometryExternalSurface,CreateImplicitGeometryExternalSurface)
####################### objects ################################
def CreateImplicitGeometryByETag(ops):
     res = ImplicitGeometryByETag()
     if "eTags" not in ops:
         raise(Exception('Need a eTags'))
     if ("support" in ops or "ls" in ops ) :
         if "ls" in ops:
             sup = ops["ls"].support
         else:
             sup = ops["support"]
         res.SetSupportAndZones(sup,PH.ReadStrings(ops["eTags"] ))
     else:
         raise(Exception('Need a (ls or support) '))

     res.offset = (PH.ReadFloat(ops.get("offset",res.offset)))
     return res

class ImplicitGeometryByETag(ImplicitGeometryBase):
    """ImplicitGeometry based in a eTags from a mesh (support) or levelset
        support : mesh from where to extract the eTags
        ls : levelset from where to extract the support to extract the eTags
        etags : list of etags
        offset : offset from the eTag (negative will grow the zone )

        Note: For the moment works only for volume and surface eTag (no lines)
        the implementation compute the distance on the the skin of the etag

    """
    def __init__(self ):
        super(ImplicitGeometryByETag,self).__init__()
        self.op = ImplicitGeometryStl()
        # we add an epsiloin to treat cases where all 4 nodes on one tetra are
        # in the iso (case of a tetra in a corrner/edge of a zone)
        self.eps = -1.e-4
        self.offset = 0.

    def SetSupportAndZones(self,support, etags):
        from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementByTags

        sup = ExtractElementByTags(support, etags)
        self.op.SetSurface(sup)

    def __str__(self):
        res = "ImplicitGeometryByETag\n"
        #res += "  support: " + str(self.support)+"\n"
        return res

    def ApplyVector(self, support,cellCenter=False):
        vals = self.op.ApplyVector(support,cellCenter=cellCenter)
        return self.ApplyInsideOut(vals+(self.eps+self.offset))

    def GetDistanceToPoint(self,_pos):
        vals = self.op.GetDistanceToPoint(_pos)
        return self.ApplyInsideOut(vals+(self.eps+self.offset))

RegisterClass("ETag",ImplicitGeometryByETag,CreateImplicitGeometryByETag)

def CreateImplicitGeometryByETagII(ops):
     res = ImplicitGeometryByETagII()
     if "eTags" not in ops:
         raise(Exception('Need a eTags'))
     if ("support" in ops or "ls" in ops ) :
         if "ls" in ops:
             sup = ops["ls"].support
         else:
             sup = ops["support"]
         res.SetSupportAndZones(sup,PH.ReadStrings(ops["eTags"] ))
     else:
         raise(Exception('Need a (ls or support) '))

     res.offset = (PH.ReadFloat(ops.get("offset",res.offset)))
     return res

class ImplicitGeometryByETagII(ImplicitGeometryBase):
    """ImplicitGeometry based in a eTags from a mesh (support) or levelset
        support : mesh from where to extract the eTags
        ls : levelset from where to extract the support to extract the eTags
        etags : list of etags
        offset : offset from the eTag (negative will grow the zone )

        Note: for surface and line the resulting fields is possitive (no inside)
        the implementation puts zeros in the interface beetwen the elements in the
        etag and the other elements (1 outside, 0 in the interface, -1 inside)
        Some special treatement is done for elements with all the nodes with zero values
        (presence of values -0.5 and 0.5 in the resulting field)

    """
    def __init__(self ):
        super(ImplicitGeometryByETagII,self).__init__()
        self.offset = 0.
        from BasicTools.Containers.Filters import ElementFilter
        self.elementFilter  = ElementFilter()
        self.offset = 0.
        self.numberOfElementSeudoDistance = 10

    def SetSupportAndZones(self,support, etags):
        self.elementFilter.mesh = support
        self.elementFilter.tags = etags

    def __str__(self):
        res = "ImplicitGeometryByETagII\n"
        res += "   " + str(self.elementFilter)
        return res

    def ApplyVector(self, support, cellCenter=False):
        if cellCenter:
            raise (Exception("Not implemented"))

        from BasicTools.Containers.Filters import ElementFilterToImplicitField
        self.elementFilter.mesh = support
        vals = ElementFilterToImplicitField(self.elementFilter,self.numberOfElementSeudoDistance)
        res = self.ApplyInsideOut(vals+(self.offset))
        return res

    def GetDistanceToPoint(self,_pos):
        if _pos is self.elementFilter.mesh.nodes:
            from BasicTools.Containers.Filters import ElementFilterToImplicitField
            vals = ElementFilterToImplicitField(self.elementFilter,self.numberOfElementSeudoDistance)
            return self.ApplyInsideOut(vals+(self.offset))
        raise


RegisterClass("ETagII",ImplicitGeometryByETagII,CreateImplicitGeometryByETagII)

class ImplicitGeometryAxisAlignBox(ImplicitGeometryBase):
    """ImplicitGeometry based ona Axis Alinged Box
        origin : origin in the mox (Xmin Ymin Zmin)
        dimensions : size of the box (Xl, Yl Zl)
    """
    def __init__(self, origin=None, size=None ):
        super(ImplicitGeometryAxisAlignBox,self).__init__()
        if origin is None:
            self.origin = np.array([0.,0.,0.], dtype= float)
        else:
            self.origin = np.array(origin, dtype= float)

        if size is None:
            self.size = np.array([1.,1.,1.], dtype= float)
        else:
            if type(size) is float:
                self.size = np.array([1.,1.,1.], dtype= float)*size
            else:
                self.size = size

    def __str__(self):
        res = "ImplicitGeometryAxisAlignBox \n"
        res += "  origin: " + str(self.origin)+"\n"
        res += "  dimensions: " + str(self.size)+"\n"
        return res

    def GetBoundingMin(self):
        return self.origin

    def GetBoundingMax(self):
        return self.origin + self.size

    def SetSupport(self,support):
        support.ComputeBoundingBox()
        self.origin = support.boundingMin
        self.size = support.boundingMax - support.boundingMin

    def GetDistanceToPoint(self,_pos):

        walls = []
        data = [[-1, 0, 0],
                [ 0,-1, 0],
                [ 0, 0,-1]]

        for normal in data:
            Obj = ImplicitGeometryPlane()
            Obj.point=self.GetBoundingMin();
            Obj.normal=np.array(normal,dtype=float);
            walls.append(Obj)

        data = [[ 1, 0, 0],
                [ 0, 1, 0],
                [ 0, 0, 1]]

        for normal in data:
            Obj = ImplicitGeometryPlane()
            Obj.point=self.GetBoundingMax();
            Obj.normal=np.array(normal,dtype=float);
            walls.append(Obj)

        return self.ApplyInsideOut(ImplicitGeometryIntersection(walls).ApplyVector(_pos))

RegisterClass("AABox",ImplicitGeometryAxisAlignBox)


def CreateImplicitGeometrySphereFromNTag(ops):
    res = ImplicitGeometrySphereFromNTag()
    res.SetRadius(PH.ReadFloat(ops.get("radius",1.)))
    res.nTag = PH.ReadString(ops["nTag"])
    res.SetSupport(ops["ls"].support)
    return res

class ImplicitGeometrySphereFromNTag(ImplicitGeometryBase):
    """ImplicitGeometry based on a ntags
        ls : levelset to extrat ntag
        nTag : name of the nTag to retrieve the centers
        radius : radius of the sphere (r)
    """
    def __init__(self,radius= 1. ):
        super(ImplicitGeometrySphereFromNTag,self).__init__()
        self.SetRadius(radius)
        self.centers = None
        self.nTag = None

    def SetLs(self,ls):
        self.SetSupport(ls.support)

    def SetSupport(self,support):
        self.centers = support.nodes[support.nodesTags[self.nTag].GetIds(),:]

    def SetRadius(self,val):
        self.radius = val

    def __str__(self):
        res = "ImplicitGeometrySphereFromNTag \n"
        res += "  radius: " + str(self.radius) + "\n"
        res += "  nTags:" + str(self.nTag) + "\n"
        return res

    def GetDistanceToPoint(self,_pos):

        if len(_pos.shape) == 1:
            pos = np.array([_pos])
        else:
            pos = _pos

        res = np.sqrt(np.sum((pos-self.centers[0,:])**2,axis=1))-self.radius
        for i in range(1,self.centers.shape[0]):
            res2 = np.sqrt(np.sum((pos-self.centers[i,:])**2,axis=1))-self.radius
            res = np.minimum(res,res2)

        return self.ApplyInsideOut(res)

RegisterClass("SphereFromNTag",ImplicitGeometrySphereFromNTag,CreateImplicitGeometrySphereFromNTag )

class ImplicitGeometrySphere(ImplicitGeometryBase):
    """ImplicitGeometry based on a sphere
        center : center of the sphere (X Y Z)
        radius : radius of the sphere (r)
    """
    def __init__(self,radius= 1.,center= None ):
        super(ImplicitGeometrySphere,self).__init__()
        if center is not None:
            self.SetCenter(PH.ReadFloats(center))
        else:
            self.SetCenter(PH.ReadFloats([0., 0., 0.]) )
        self.SetRadius(PH.ReadFloat(radius))

    def __str__(self):
        res = "ImplicitGeometrySphere \n"
        res += "  center: " + str(self.center)+"\n"
        res += "  radius: "+str(self.radius)
        return res

    def SetCenter(self,center):
        self.center = PH.ReadFloats(center)

    def SetRadius(self,val):
        self.radius = PH.ReadFloat(val)

    def GetDistanceToPoint(self,_pos):

        if len(_pos.shape) == 1:
            pos = np.array([_pos])
        else:
            pos = _pos

        res = np.sqrt(np.sum((pos-self.center)**2,axis=1))-self.radius

        return self.ApplyInsideOut(res)

    def GetBoundingMin(self):
        res = self.center -  self.radius
        return res

    def GetBoundingMax(self):
        res = self.center + self.radius
        return res

RegisterClass("Sphere",ImplicitGeometrySphere )


class ImplicitGeometryCylinder(ImplicitGeometryBase):
    """ImplicitGeometry based on a cylinder
        center1 : first point on the axis (X Y Z)
        center2 : second point on the axis (X Y Z)
        radius : radius of the cylinder (r)
        wcups : if we compute the distance to the cups of the cylinder
    """
    def __init__(self,center1=None,center2=None,radius=1.,wcups=True):
        super(ImplicitGeometryCylinder,self).__init__()
        if center1 is None:
           self.center1 = np.zeros((3,))
        else:
           self.center1 = np.array(center1,copy=True,dtype=np.float)
        if center2 is None:
           self.center2 = np.array([1.,0.,0.])
        else:
           self.center2 = np.array(center2,copy=True,dtype=np.float)

        self.radius = radius
        self.wcups = wcups
        #self.infinit = infinit

    def __str__(self):
        res = "ImplicitGeometryCylinder \n"
        res += "  center1: " + str(self.center1)+"\n"
        res += "  center2: " + str(self.center2)+"\n"
        res += "  radius: "+str(self.radius)
        res += "  width cups: "+str(self.wcups)
        #res += "  infinit: "+str(self.infinit)
        return res

    def GetDistanceToPoint(self,_pos):
        if len(_pos.shape) == 1:
            pos = np.array([_pos])
        else:
            pos = _pos

        u = (self.center2 -self.center1).astype(np.float) # director vector
        nu = np.linalg.norm(u)
        u /= nu

        # distance to the center of the cylinder
        d = np.linalg.norm(np.cross(pos-self.center1,u),axis=1) - self.radius

        if self.wcups:
            # distance to the planes (the cups)
            pA1 = ImplicitGeometryPlane(point=self.center1, normal=-u)
            pA2 = ImplicitGeometryPlane(point=self.center2, normal=u)
            pA =  ImplicitGeometryIntersection([pA1,pA2]).GetDistanceToPoint(pos)
            res = np.maximum(d,pA)
        else:
            res = d

        return self.ApplyInsideOut(res)

    def GetBoundingMin(self):
        res = np.minimum(self.center1, self.center2)
        diff = self.center2 -self.center1
        l = math.sqrt(np.linalg.norm(diff))
        if diff[0] > 0:
            res[0] = res[0]-l*self.radius/abs(diff[0])
        if diff[1] > 0:
            res[1] = res[1]-l*self.radius/abs(diff[1])
        if diff[2] > 0:
            res[2] = res[2]-l*self.radius/abs(diff[2])
        return res

    def GetBoundingMax(self):
        res = np.maximum(self.center1, self.center2)
        diff = self.center2 - self.center1
        l = math.sqrt(np.linalg.norm(diff))
        if diff[0] > 0:
            res[0] = res[0]+l*self.radius/abs(diff[0])
        if diff[1] > 0:
            res[1] = res[1]+l*self.radius/abs(diff[1])
        if diff[2] > 0:
            res[2] = res[2]+l*self.radius/abs(diff[2])
        return res

RegisterClass("Cylinder",ImplicitGeometryCylinder )

def CreateImplicitGeometryPlane(ops):
    obj = ImplicitGeometryPlane()
    PH.ReadProperties(ops,["point","offset"], obj)
    obj.SetNormal(PH.ReadFloats(ops["normal"]))
#    _setProps = ["normal":np.array([1.,0,0])]
    return obj


class ImplicitGeometryPlane(ImplicitGeometryBase):
    """ImplicitGeometry based on a plane
        point : point on the plane (X Y Z)
        normal : normale of the plane (X Y Z)
        offset : offset to the plane  (r)
    """
    def __init__(self, point=None, normal=None, offset = 0.0):
        super(ImplicitGeometryPlane,self).__init__()

        if point is None:
            self.point = np.array([0.,0.,0.],dtype=float,copy=True)
        else:
            self.point = np.array(point,copy=True)

        if normal is None:
            self.SetNormal(np.array([1,0,0], dtype =float ) )
        else :
            self.SetNormal(normal)

        self.offset = float(offset)

    def SetNormal(self,invec):
        self.normal = np.array(invec,copy=True)
        self.normal = self.normal/np.linalg.norm(self.normal)

    #normal = property(fset=SetNormal)


    def __str__(self):
        res = "ImplicitGeometryPlane\n"
        res += "  point: " + str(self.point) + "\n"
        res += "  normal: " + str(self.normal) + "\n"
        res += "  offset: " + str(self.offset) + "\n"
        return res

    def GetDistanceToPoint(self,pos):
        d = self.normal.dot(self.point)
        res = np.sum(pos*self.normal,axis=1) - d - self.offset
        #res = np.sum(self.normal*(pos-self.point),axis=1)
        return self.ApplyInsideOut(res)

RegisterClass("Plane",ImplicitGeometryPlane,CreateImplicitGeometryPlane)

class ImplicitGeometryGyroid(ImplicitGeometryBase):
    """ImplicitGeometry based on a gyroid
        scale : scale of the gyroid (1. defaul)
        offset : offset of the gyroid (X Y Z)
        type : [0-1] version of gyroid (0 default)
        wall : [True-False] wall or body (false default)
        wallThickness : parameter to control the wall thickness
    """

    def __init__(self, scale=1., offset=[0.,0.,0.] ):
        super(ImplicitGeometryGyroid,self).__init__()
        self.scale = scale
        self.offset = np.array(offset,dtype=float)
        self.type = 0
        self.wall = False;
        self.wallThickness = 0.5

    def __str__(self):
        res = "ImplicitGeometryGyroid \n"
        res += "  scale: " + str(self.scale)+"\n"
        res += "  offset: "+str(self.offset)+"\n"
        return res

    def GetDistanceToPoint(self,pos):

        npos = pos/self.scale +self.offset

        x = npos[:,0]*np.pi
        y = npos[:,1]*np.pi
        z = npos[:,2]*np.pi

        if self.type == 0:
            res = (dsin(x)*dcos(y) + dsin(y)*dcos(z) + dsin(z)*dcos(x))
            res *= self.scale/(1.55*np.pi)
        elif self.type == 1:
            res = np.sin(x)*np.cos(y) + np.sin(y)*np.cos(z) + np.sin(z)*np.cos(x)
            #sin = np.sin
            #cos = np.cos
            #ng = 0.01+np.sqrt(abs(sin(x)*sin(y) - cos(y)*cos(z))**2 + abs(sin(x)*sin(z) - cos(x)*cos(y))**2 + abs(sin(y)*sin(z) - cos(x)*cos(z))**2)
            #res /= ng
            # aproximation of the signed distance  dist = f/norm(grad(f))

        else:
            res = 0.5*(+np.sin(2.*x)*np.cos(y)*np.sin(z)
                       +np.sin(2.*y)*np.cos(z)*np.sin(x)
                       +np.sin(2.*z)*np.cos(x)*np.sin(y)
                      )-0.5*(
                       +np.cos(2.*x)*np.cos(2.*y)
                       +np.cos(2.*y)*np.cos(2.*z)
                       +np.cos(2.*z)*np.cos(2.*x)
                      )+0.15;

        if self.wall:
            res = np.abs(res)-self.wallThickness/2.

        return self.ApplyInsideOut(res)

RegisterClass("Gyroid",ImplicitGeometryGyroid)

class ImplicitGeometry60D(ImplicitGeometryBase):
    """ImplicitGeometry using walls at 60 degrees in the (x,y) plane
        lx : scale
        w : wall thickness
    """
    def __init__(self, lx=1., w=0.1 ):
        super(ImplicitGeometry60D,self).__init__()
        self.lx = lx
        self.w = w

    def __str__(self):
        res = "ImplicitGeometry60D \n"
        res += "  lx: " + str(self.lx)+"\n"
        res += "  width: "+str(self.w)+"\n"
        return res

    def GetDistanceToPoint(self,pos):

        v = math.sqrt(3.)

        llx = self.lx/2.
        lw =self.w/2.

        aw = lw/llx;

        mpos = abs((np.mod(pos/llx,[2.,2*v,1.])-[1,v,0]) )

        pA1 = ImplicitGeometryPlane(point=np.array([0.,0,0.]), normal=np.array([0,1.,0,]), offset = aw )

        normal = np.array([-v,1.,0])
        pB1 =ImplicitGeometryPlane(point=[0.,0,0], normal=-normal,  offset= aw)
        pB2 =ImplicitGeometryPlane(point=[0.,0,0], normal= normal,  offset= aw)
        pB =  ImplicitGeometryIntersection([pB1,pB2])

        pC1 =ImplicitGeometryPlane(point=[0,v,0], normal=np.array([0,-1,0,]), offset = aw)

        res = ImplicitGeometryUnion([pA1,pB,pC1]).GetDistanceToPoint(mpos)

        return self.ApplyInsideOut(res)

RegisterClass("60D",ImplicitGeometry60D)


class ImplicitGeometryHoneycomb(ImplicitGeometryBase):
    """ImplicitGeometry Honycomb  (in the xy plane)
        lx : scale (1)
        w : wall thickness (0.1)
        translate : offset of (0 0 0)
    """

    def __init__(self, lx=1., w=0.1 ):
        super(ImplicitGeometryHoneycomb,self).__init__()
        self.lx = float(lx)
        self.w = float(w)
        self.translate = np.array([0.,0.,0.], dtype= float)

    def __str__(self):
        res = "ImplicitGeometryHoneycomb \n"
        res += "  lx: " + str(self.lx)+"\n"
        res += "  width: "+str(self.w)+"\n"
        res += "  translate: "+str(self.translate)+"\n"
        return res

    def GetDistanceToPoint(self,pos):

        v = math.sqrt(3.)
        v2 = math.sqrt(1./3.)

        llx = self.lx/2.
        lw =self.w/2.

        mpos = abs((np.mod((pos+self.translate)/llx,[2.,2*v,1.])-[1.,v,0.]) )

        pA1 = ImplicitGeometryPlane(point=np.array([(lw/llx),0.,0.]), normal=np.array([1.,0.,0.,]))
        pA2 = ImplicitGeometryPlane(point=[0.,v2,0.], normal=[0.,1.,0.])
        pA =  ImplicitGeometryIntersection([pA1,pA2])


        normal = np.array([-v2,1.,0.])
        pB1 =ImplicitGeometryPlane(point=[0.,v2,0.], normal=-normal, offset= lw/llx)
        pB2 =ImplicitGeometryPlane(point=[0.,v2,0.], normal=normal, offset= lw/llx)
        pB =  ImplicitGeometryIntersection([pB1,pB2])

        pC1 =ImplicitGeometryPlane(point=[1-(lw/llx),0.,0.], normal=np.array([-1.,0.,0.,]))
        pC2 =ImplicitGeometryPlane(point=[1.,v-v2,0], normal=[0.,-1.,0.,])
        pC =  ImplicitGeometryIntersection([pC1,pC2])

        #res = ImplicitGeometryUnion([pA,pC]).GetDistanceToPoint(mpos)
        res = ImplicitGeometryUnion([pA,pB,pC]).GetDistanceToPoint(mpos)
        #res = ImplicitGeometryUnion([pA,pB,pC]).GetDistanceToPoint(mpos)

        #res = ImplicitGeometryUnion([pA,pB,pC]).GetDistanceToPoint(mpos)

        return self.ApplyInsideOut(res)

RegisterClass("Honeycomb",ImplicitGeometryHoneycomb)


import os
## hack to make work loading stl using  point as decimal separator (,/.)
os.environ["LANG"] = "en_UK"

class ImplicitGeometryStl(ImplicitGeometryBase):
    """ImplicitGeometry based on a external stlfile
        filename : stl filename to be loaded
    """

    def __init__(self):
        super(ImplicitGeometryStl,self).__init__()
        self.implicitFunction = None

        self.filename = ''
        self.boundingMin = [0,0,0];
        self.boundingMax = [0,0,0];
        self.surface = None
        self.ismanifold = True
        self.onLines = False



    def GetBoundingMin(self):
        return self.boundingMin

    def GetBoundingMax(self):
        return self.boundingMax

    def SetFileName(self,filenameSTL):
        self.LoadFromFile(filenameSTL)

    def LoadFromFile(self,filenameSTL):

        import vtk
        if filenameSTL.split(".")[-1] == "stl":
            readerSTL = vtk.vtkSTLReader()
        else:
            readerSTL = vtk.vtkXMLPolyDataReader()

        readerSTL.SetFileName(filenameSTL)

        self.filename = filenameSTL
        # 'update' the reader i.e. read the .stl file
        readerSTL.Update()
        polydata = readerSTL.GetOutput()

        if polydata.GetNumberOfPoints() == 0:# pragma: no cover
            raise ValueError( "No point data could be loaded from '" + filenameSTL)

        self.SetSurfaceUsingVtkPolyData(polydata)
        return self

    def SetMesh(self,mesh):
        # check if we have only element of dimensionality 2 or less
        from BasicTools.Containers.ElementNames import dimension
        for name,data in mesh.elements.items():
            if data.GetNumberOfElements() == 0: continue
            if dimension[name] == 3:
                break
        else:
             return self.SetSurface(mesh)


        import vtk
        from BasicTools.Containers.vtkBridge import MeshToVtk
        vtkmesh = MeshToVtk(mesh)

        filt = vtk.vtkDataSetSurfaceFilter()
        filt.SetInputData(vtkmesh)
        filt.Update()
        self.SetSurfaceUsingVtkPolyData(filt.GetOutput())


    def SetSurface(self,mesh):
        # check if we have element of dimensionality 3
        from BasicTools.Containers.ElementNames import dimension
        for name,data in mesh.elements.items():
            if data.GetNumberOfElements() == 0: continue
            if dimension[name] >2:
                return self.SetMesh(mesh)

        from BasicTools.Containers.vtkBridge import MeshToVtk
        vtkmesh = MeshToVtk(mesh)
        self.SetSurfaceUsingVtkPolyData(vtkmesh)

    def SetSurfaceUsingVtkPolyData(self,polydata):
        import vtk

        if polydata.GetNumberOfPoints() == 0:# pragma: no cover
            raise ValueError( "No points " )

        bounds = [0 for i in range(6)]
        polydata.GetBounds(bounds)
        self.boundingMin = np.array([bounds[0],  bounds[2], bounds[4]])
        self.boundingMax = np.array([bounds[1],  bounds[3], bounds[5]])

        # check if we have a closed surface or not
        checkFilter = vtk.vtkFeatureEdges()
        checkFilter.SetInputData(polydata)
        checkFilter.SetFeatureEdges(False)
        checkFilter.SetBoundaryEdges(True)
        checkFilter.SetNonManifoldEdges(True)
        checkFilter.SetManifoldEdges(False)
        checkFilter.Update()

        numberOfOpenEdges = checkFilter.GetOutput().GetNumberOfCells();
        if numberOfOpenEdges > 0  :
            self.ismanifold = False
        else:
            self.ismanifold = True

        if polydata.GetNumberOfPolys() > 0:
            self.implicitFunction = vtk.vtkImplicitPolyDataDistance()
            self.implicitFunction.SetNoValue(-100.)
            self.implicitFunction.SetInput(polydata);
            self.onLines = False
        elif polydata.GetNumberOfLines() >0:
            #we are working on lines only
            # The current installation of tk does not have thevtkSpheres() class
            self.implicitFunction = vtk.vtkSpheres()
            self.implicitFunction.SetCenters(polydata.GetPoints())
            self.onLines = True
            self.ismanifold = False
        else:
            raise(Exception("internal Error"))

    def GetDistanceToPoint(self,pos):
        if len(pos.shape) == 1:
            res = np.zeros(1,dtype=np.float)
        else:
            res = np.zeros(pos.shape[0],dtype=np.float)

        if len(pos.shape) == 1:
            res[0]  =  self.implicitFunction.EvaluateFunction(pos)
        else:
            res = np.zeros(pos.shape[0],dtype=np.float)
            cpt =0;
            self.PrintDebug("in ImplicitGeometryStl::GetDistanceToPoint")
            for point in pos:
               res[cpt]  =  self.implicitFunction.EvaluateFunction(point)
               cpt +=1;

        if not self.ismanifold:
            np.abs(res,out=res)

        return self.ApplyInsideOut(res)

    def __str__(self):
        res = "ImplicitGeometryStl \n"
        res += "  fileName : " + str(self.filename)+"\n"
        res += "  InsideOut : " + str(self.insideOut)+"\n"
        res += "  surface : " + str(self.surface)+"\n"
        return res

def CreateImplicitGeometryStl(ops):
       res =  ImplicitGeometryStl()
       if "filename" in ops:
           res.LoadFromFile(ops["filename"])
       return res

RegisterClass("StlFile",ImplicitGeometryStl,CreateImplicitGeometryStl)

class ImplicitGeometryAnalytical(ImplicitGeometryBase):
    """ImplicitGeometry based on a function f(x,y,z) or f(pos)
        expr : function expression
    """
    def __init__(self):
        super(ImplicitGeometryAnalytical,self).__init__()
        self.expression = "0"

    def SetExpression(self,string):
        self.expression = string

    def GetDistanceToPoint(self,pos):

        if len(pos.shape) == 1:
            res = np.zeros(1,dtype=np.float)
            x, y, z = pos
        else:
            res = np.zeros(pos.shape[0],dtype=np.float)
            x = pos[:,0]
            y = pos[:,1]
            z = pos[:,2]

        from sympy import sympify
        ls = {"x":x,"y":y,"z":z,"pos":pos}
        res = np.array(sympify(self.expression,ls),dtype=float)
        return self.ApplyInsideOut(res)

    def __str__(self):
        res = "ImplicitGeometryAnalytical \n"
        res += "  expression : " + str(self.expression)+"\n"
        res += "  InsideOut : " + str(self.insideOut)+"\n"
        return res

def CreateImplicitGeometryAnalytical(ops):
       res =  ImplicitGeometryAnalytical()
       if "expr" in ops:
           res.SetExpression(ops["expr"])
       return res

RegisterClass("Analytical",ImplicitGeometryAnalytical,CreateImplicitGeometryAnalytical)


def InitHoles(ls, nx, ny, nz, r):
    nNodes = ls.support.GetDimensions()
    grids = np.ix_(
            np.linspace(0.0, 2 * nx * np.pi, nNodes[0]),
            np.linspace(0.0, 2 * ny * np.pi, nNodes[1]),
            np.linspace(0.0, 2 * nz * np.pi, nNodes[2]))
    ls.phi = -np.cos(grids[0]) * np.cos(grids[1]) * np.cos(grids[2]) + r - 1.0
    #self.phi = 0.2 * np.ceil(np.maximum(self.phi, 0.0)) - 0.1
    ls.phi.shape = (np.prod(ls.phi.shape),)

class ImplicitGeometryHoles(ImplicitGeometryBase):
    """ImplicitGeometry to create holes
        n : number of hole in each direction (3 3 3)
        r : radius of holes (0.5)
        offset : offset (0 0 0 )
        support : (grid) to compute the bounding box to generate the holes
    """
    def __init__(self):
        super(ImplicitGeometryHoles,self).__init__()
        self.r = 0.5
        self.n = np.array([3.,3.,3.])
        self.offset = np.array([0.,0.,0.])
        self.type= ''
        self.boundingMin = np.array([0.,0.,0.])
        self.boundingMax = np.array([1.,1.,1.])

    def SetSupport(self,support):
        support.ComputeBoundingBox()
        self.boundingMin = support.boundingMin[:]
        self.boundingMax = support.boundingMax[:]

    def GetBoundingMin(self):
        return self.boundingMin

    def GetBoundingMax(self):
        return self.boundingMax

    def SetNumberOfHoles(self,data):
        if len(data) != 3:
            raise
        self.n = np.array(data,dtype=int)

    def GetDistanceToPoint(self,pos):
        l = self.boundingMax - self.boundingMin

        l[l ==0] = 1.
        l = list(l)
        if len(l) == 2:
            l.append(1.)

        if self.type == "Original":
            if len(pos.shape) == 1:
                res = np.zeros(1,dtype=np.float)
                res[0] = -np.cos(pos[0]) * np.cos(pos[1]) * np.cos(pos[2]) + self.r - 1.0
            else:
                res = -np.cos(self.offset[0] + (self.n[0]*np.pi/l[0])*pos[:,0]) * \
                       np.cos(self.offset[1] + (self.n[1]*np.pi/l[1])*pos[:,1]) * \
                       np.cos(self.offset[2] + (self.n[2]*np.pi/l[2])*pos[:,2]) + self.r - 1.0
        else:
            dl = l/self.n
            mpos = abs(np.mod((pos+self.offset)+dl,2*dl))-dl
            sA = ImplicitGeometrySphere(radius=self.r,center=[0,0,0])
            sB0 = ImplicitGeometrySphere(radius=self.r,center=[dl[0],dl[1],0])
            sB1 = ImplicitGeometrySphere(radius=self.r,center=[-dl[0],-dl[1],0])
            sB2 = ImplicitGeometrySphere(radius=self.r,center=[dl[0],-dl[1],0])
            sB3 = ImplicitGeometrySphere(radius=self.r,center=[-dl[0],dl[1],0])

            sC0 = ImplicitGeometrySphere(radius=self.r,center=[dl[0],0,dl[2]])
            sC1 = ImplicitGeometrySphere(radius=self.r,center=[-dl[0],0,-dl[2]])
            sC2 = ImplicitGeometrySphere(radius=self.r,center=[-dl[0],0,dl[2]])
            sC3 = ImplicitGeometrySphere(radius=self.r,center=[dl[0],0,-dl[2]])

            sD0 = ImplicitGeometrySphere(radius=self.r,center=[0,dl[1],dl[2]])
            sD1 = ImplicitGeometrySphere(radius=self.r,center=[0,-dl[1],-dl[2]])
            sD2 = ImplicitGeometrySphere(radius=self.r,center=[0,-dl[1],dl[2]])
            sD3 = ImplicitGeometrySphere(radius=self.r,center=[0,dl[1],-dl[2]])

            Ores = ImplicitGeometryUnion([sA,sB0,sB1,sB2,sB3,sC0,sC1,sC2,sC3,sD0,sD1,sD2,sD3])
            Ores.insideOut=True
            res = Ores.GetDistanceToPoint(mpos)

        return self.ApplyInsideOut(res)


    def __str__(self):
        res = "ImplicitHoles \n"
        res += "  self.r : " + str(self.r)+"\n"
        res += "  InsideOut : " + str(self.insideOut)+"\n"
        return res

RegisterClass("Holes",ImplicitGeometryHoles)



def CheckIntegrity(GUI=False):

    def MustFail(func):
        try:
            func()
            raise #pragma: no cover
        except:
            pass

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([5,6,7]);
    myMesh.SetSpacing(1./(myMesh.GetDimensions()-1)*2);
    myMesh.SetOrigin([-1.,-1.,-1.]);
    myMesh.elements["hex8"].tags.CreateTag("2elems").SetIds([0,1])
    myMesh.nodesTags.CreateTag("3points").SetIds([2,3,4])
    print(myMesh)

    OnePoint3D = np.array([1,2,3])
    TwoPoints3D = np.array([[0.,0.,0.],[1.,2.,3.]], dtype=np.float)

    import BasicTools.TestData as TestData
    from BasicTools.Helpers.Tests import TestTempDir

    testDataPath = TestData.GetTestDataPath()

    ##################### ImplicitGeometrySphere ################################
    IGSphere = ImplicitGeometrySphere(center=[0,0,0],radius=0.1)
    IGSphere = ImplicitGeometrySphere()
    IGSphere.center = np.array([.5,-0.3,0.5])
    IGSphere.radius = 0.7
    IGSphere_Data = IGSphere.ApplyVector(myMesh)
    print(IGSphere)
    print(IGSphere.GetBoundingMin()  )
    print(IGSphere.GetBoundingMax()  )
    IGSphere(np.array([1,2,3]))
    ###########################################################################
    IGWrapped = ImplicitGeometryWrapped(IGSphere_Data)
    IGWrapped(myMesh)

    ImplicitGeometryWrapped(TwoPoints3D[:,0])(TwoPoints3D)
    from functools import partial
    MustFail(partial(ImplicitGeometryWrapped(TwoPoints3D[0:1,0]),TwoPoints3D))
    MustFail(partial(ImplicitGeometryWrapped(TwoPoints3D[0:1,0]),myMesh,cellCenter=True))
    MustFail(partial(ImplicitGeometryWrapped(TwoPoints3D[0:1,0]),myMesh,cellCenter=False))

    #########################ImplicitGeometryExternalSurface#########################
    IGExternalSurface = CreateImplicitGeometryExternalSurface({"support":myMesh})
    IGExternalSurface(TwoPoints3D)
    ############################ImplicitGeometryByETag###############################
    IGByETag = CreateImplicitGeometryByETag({"support":myMesh, "eTags":["2elems"]})
    print(IGByETag)
    IGByETag(TwoPoints3D)
    IGByETag.GetDistanceToPoint(TwoPoints3D)

    class LS():
        pass
    LS.support  = myMesh
    IGByETag = CreateImplicitGeometryByETag({"ls":LS, "eTags":["2elems"]})

    MustFail(partial(CreateImplicitGeometryByETag,{"support":myMesh}  ))
    MustFail(partial(CreateImplicitGeometryByETag,{"eTags":["2elems"]}  ))
    #########################  ImplicitGeometryAxisAlignBox ###################
    IGAxisAlignBox = ImplicitGeometryAxisAlignBox()
    IGAxisAlignBox = ImplicitGeometryAxisAlignBox(origin=[0,1,2],size=[0.1,0.2,0.3])
    IGAxisAlignBox = ImplicitGeometryAxisAlignBox(origin=[0,1,2],size=3.)
    IGAxisAlignBox.GetBoundingMin()
    IGAxisAlignBox.GetBoundingMax()
    IGAxisAlignBox.SetSupport(myMesh)
    IGAxisAlignBox(myMesh)
    print(IGAxisAlignBox)
    ###################### ImplicitGeometrySphereFromNTag  ###################

    IGSphereFromNTag= CreateImplicitGeometrySphereFromNTag({"radius":0.1,"nTag":"3points","ls":LS})
    IGSphereFromNTag.SetLs(LS)
    IGSphereFromNTag(myMesh)
    IGSphereFromNTag(OnePoint3D)
    print(IGSphereFromNTag)

    ###################### ImplicitGeometryCylinder  ###################

    IGCylinder =  ImplicitGeometryCylinder(wcups=False)
    IGCylinder(OnePoint3D)
    IGCylinder =  ImplicitGeometryCylinder(center1=[0,0,0], center2=[1,2,3],radius=0.1)
    IGCylinder(TwoPoints3D)
    print(IGCylinder)
    IGCylinder.GetBoundingMin()
    IGCylinder.GetBoundingMax()

    ###################### ImplicitGeometryPlane  ###################


    IGGeometryPlane = CreateImplicitGeometryPlane({"point":[0,0,0],"normal":[1,1,1]})
    print(IGGeometryPlane)

    ###################### ImplicitGeometryGyroid  ###################
    IGGyroid = ImplicitGeometryGyroid()
    IGGyroid = ImplicitGeometryGyroid(scale=0.1,offset=[0.,0.2,0.3])
    print(IGGyroid)
    IGGyroid(myMesh)
    IGGyroid.type = 1
    IGGyroid(myMesh)
    IGGyroid.type = 2
    IGGyroid(myMesh)
    IGGyroid.wall = True
    IGGyroid(myMesh)

    ###################### ImplicitGeometry60D ###################
    IG60D = ImplicitGeometry60D()
    print(IG60D)
    IG60D(myMesh)

    ###################### ImplicitGeometryHoneycomb ###################
    IGHoneycomb = ImplicitGeometryHoneycomb()
    print(IGHoneycomb)
    IGHoneycomb(myMesh)

    ###################### ImplicitGeometryStl ###################

    IGStl = CreateImplicitGeometryStl({"filename":testDataPath+"stlsphere.stl"})
    IGStl.SetFileName(testDataPath+"stlsphere.stl")
    IGStl.GetBoundingMin()
    IGStl.GetBoundingMax()
    IGStl(myMesh)
    IGStl.SetFileName(testDataPath+"vtkPolySphere.vtp")
    print(IGStl)

    ###################### ImplicitGeometryAnalytical ###################

    IGAnalytical = CreateImplicitGeometryAnalytical({"expr":"2*y+1+z"})
    IGAnalytical(myMesh)
    IGAnalytical(OnePoint3D)
    print(IGAnalytical)


    ###################### ImplicitGeometryHoles ###################
    IGHoles = ImplicitGeometryHoles()
    IGHoles.SetSupport(myMesh)
    IGHoles.GetBoundingMin()
    IGHoles.GetBoundingMax()
    IGHoles.SetNumberOfHoles([2,3,4])
    IGHoles(myMesh)
    print(IGHoles)
    IGHoles.type = "Original"
    IGHoles(myMesh)
    IGHoles(OnePoint3D)

    MustFail(partial(IGHoles.SetNumberOfHoles,[1,2]))
    ############################ImplicitGeometryByETagII###############################
    IGByETagII = CreateImplicitGeometryByETagII({"support":myMesh, "eTags":["2elems"]})
    print(IGByETagII)
    ETagII = IGByETagII(myMesh)

    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=True))

