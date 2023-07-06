# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class Transform(BaseOutputObject):
    def __init__(self, offset=None, first=None, second=None):
       super(Transform,self).__init__()
       self.offset  = np.array([0.0, 0.0, 0.0], dtype=PBasicFloatType)
       self.RMatrix = np.array([[1.0,0,0],[0,1,0],[0,0,1]], dtype=PBasicFloatType)
       self.keepOrthogonal = True
       self.keepNormalised = True

       # offset off the new origin with respect to the old
       if offset is not None:
           self.SetOffset(offset)

       if first is not None:
           self.SetFirst(first)

       if second is not None:
           self.SetSecond(second)

    def GetDirection(self,i, pos=None, direction=None):
        return self.RMatrix[i,:]

    def SetOffset(self,data):
        # this point define the origin of the new coordinate system
        self.offset  = np.array(data, dtype=PBasicFloatType)

    def SetFirst(self,data):
        # this point define the x coordinate (direction) with respect to the new origin
        first  = np.array(data, dtype=PBasicFloatType)
        if self.keepNormalised :
            first /= np.linalg.norm(first)
        self.RMatrix[0,:] = first
        self.SetSecond()
        self.SetThird()

    def SetSecond(self,data=None):
        # this point define the y coordinate (direction) with respect to the new origin
        # the z direction is calculated with a cross product
        first = self.RMatrix[0,:]

        if data is None:
            second = self.RMatrix[1,:]-first*np.dot(first,self.RMatrix[1,:])/np.dot(first,first)
            # if norm is zero compute a vector not colinear to first and we restart
            if np.linalg.norm(second) == 0:
                second = first*1
                for i,v in enumerate(second):
                    if v == 0:
                        second[i] += 1
                        break
                else:
                    second[0] += 1
                return self.SetSecond(second)
        else:
            second = np.array(data, dtype=PBasicFloatType)

        if self.keepOrthogonal:
            second -= first*np.dot(first,second)/np.dot(first,first)

        second_norm = np.linalg.norm(second)
        if  second_norm == 0: # pragma: no cover
            raise Exception("cant set second to " + str(data) + " This vector ir colinear to first :" + str(first))

        if self.keepNormalised :
            second /= second_norm
        self.RMatrix[1,:] = second
        self.SetThird()

    def SetThird(self,data=None):
        first = self.RMatrix[0,:]
        second = self.RMatrix[1,:]

        if data is None:
            third = np.cross(first, second)
        else:
            third = np.array(data, dtype=PBasicFloatType)
            if self.keepOrthogonal:
                third -= first*np.dot(first,third)/np.dot(first,first)
                third -= second*np.dot(second,third)/np.dot(second,second)

        third_norm = np.linalg.norm(third)
        if  third_norm == 0: # pragma: no cover
            raise Exception("cant set third to " + str(data) + " This vector ir colinear to first and/or sencond :" + str(first) + " " + str(second))

        if self.keepNormalised :
            third /= third_norm
        self.RMatrix[2,:] = third

    def SetOpUsingThird(self,third):
        self.SetFirst(third)
        self.SetOperator(op=np.roll(self.RMatrix, -1, axis=0))

    def SetOperator(self, first=None, second=None, third=None, op=None)    :
        if op is None:
            self.SetFirst(first)
            self.SetSecond(second)
            self.SetThird(third)
        else:
            if (first is not None) or (second is not None) or  (third is not None): # pragma: no cover
                raise(Exception("Cant define operanter and direction verctor at the same time") )
            self.SetOperator(first=op[0,:],second=op[1,:], third=op[2,:])

    def ApplyTransform(self,point):
        #p = M*(point-self.offset)
        point = np.asarray(point)
        op = self.RMatrix
        return self.__ApplyOpForEveryLine(op,point-self.offset)

    def ApplyInvTransform(self,point):
        # we apply inverse of the transformation
        #p = M^-1*point+self.offset
        point = np.asarray(point)

        if self.keepNormalised == False or self.keepOrthogonal == False:
            op = np.linalg.inv(self.RMatrix)
        else:
            op = self.RMatrix.T

        return self.__ApplyOpForEveryLine(op,point)+self.offset

    def ApplyTransformDirection(self,point):
        point = np.asarray(point)
        op = self.RMatrix
        return self.__ApplyOpForEveryLine(op,point)

    def ApplyInvTransformDirection(self,point):
        # we apply inverse of the transformation
        #p = point+self.offset
        point = np.asarray(point)

        if self.keepNormalised == False or self.keepOrthogonal == False:
            op = np.linalg.inv(self.RMatrix)
        else:
            op = self.RMatrix.T

        return self.__ApplyOpForEveryLine(op,point)

    def ApplyTransformTensor(self,tensor):
        tensor = np.asarray(tensor)
        op = self.RMatrix
        return np.dot(op,np.dot(tensor,op.T))

    def ApplyInvTransformTensor(self,tensor):
        # we apply inverse of the transformation

        if self.keepNormalised == False or self.keepOrthogonal == False:
            op = np.linalg.inv(self.RMatrix)
        else:
            op = self.RMatrix.T
        return np.dot(op,np.dot(tensor,op.T))


    def __ApplyOpForEveryLine(self,op,points):

        if len(points.shape) == 2:
            return np.dot(op,points.T).T
        else:
            return np.dot(op,points)


    def GetOrthoNormalBase(self):
        return Transform(self.offset, self.RMatrix[0,:], self.RMatrix[1,:])

    def __str__(self):
        res = "Transform : \n"
        res += "    keepNormalised : " + str(self.keepNormalised)+ "\n"
        res += "    keepOrthogonal : " + str(self.keepOrthogonal)+ "\n"
        res += "    offset : " + str(self.offset)+ "\n"
        res += "    first  : " + str(self.RMatrix[0,:])+ "\n"
        res += "    second : " + str(self.RMatrix[1,:])+ "\n"
        res += "    third : " + str(self.RMatrix[2,:]) + "\n"
        return res

class Transform2D(BaseOutputObject):

    def __init__(self, offset=None, first=None, second=None):
        super().__init__()
        offset = self.EnsureIn3D(offset)
        first = self.EnsureIn3D(first)
        second = self.EnsureIn3D(second)
        self.transform3D = Transform(offset=offset, first=first, second=second)

    def SetOffset(self,data):
        self.transform3D.SetOffset(self.EnsureIn3D(data))

    def SetFirst(self,data):
        self.transform3D.SetFirst(self.EnsureIn3D(data))

    @property
    def keepOrthogonal(self):
        return self.transform3D.keepOrthogonal

    @property
    def keepNormalised(self):
        return self.transform3D.keepNormalised

    @keepOrthogonal.setter
    def keepOrthogonal(self, value):
        self.transform3D.keepOrthogonal = value

    @keepNormalised.setter
    def keepNormalised(self, value):
        self.transform3D.keepNormalised = value

    def GetDirection(self,i, pos=None, direction=None):
        return self.EnsureIn2D(self.transform3D.GetDirection(pos, direction))

    def ApplyTransform(self, point):
        return self.EnsureIn2D(self.transform3D.ApplyTransform(self.EnsureIn3D(point)))

    def ApplyInvTransform(self, point):
        return self.EnsureIn2D(self.transform3D.ApplyInvTransform(self.EnsureIn3D(point)))

    def ApplyTransformDirection(self, point):
        return self.EnsureIn2D(self.transform3D.ApplyTransformDirection(self.EnsureIn3D(point)))

    def ApplyInvTransformDirection(self, point):
        return self.EnsureIn2D(self.transform3D.ApplyInvTransformDirection(self.EnsureIn3D(point)))

    def GetOrthoNormalBase(self):
        return Transform2D(self.transform3D.offset, self.transform3D.RMatrix[0,:], self.transform3D.RMatrix[1,:])

    def EnsureIn2D(self, data):
        if data is None:
            return None
        data = np.asarray(data)

        if len(data.shape) == 1:
            return data[:2]
        else:
            return np.ascontiguousarray(data[:,:2])

    def EnsureIn3D(self, data):

        if data is None:
            return None
        data = np.asarray(data)

        if len(data.shape) == 1:
            if data.shape[0] == 3:
                return data
            res = np.empty((3),dtype=data.dtype)
            res[:data.shape[0]] = data
            res[data.shape[0]:] = 0
            return res
        else:
            if data.shape[1] == 3:
                return data
            res = np.empty((data.shape[0],3),dtype=data.dtype)
            res[:,:data.shape[1]] = data
            res[:,data.shape[1]:] =0
            return res


def CheckIntegrity(GUI=False):
    orient = Transform()
    orient.SetOpUsingThird([5,2,1])
    orient.SetOperator(op=orient.RMatrix)

    orient = Transform(offset=[1,1,1],
                        first=[1,0,0],
                        second =[0,1,0])

    print(orient.GetDirection(0) )
    #test with a list of points
    orient.ApplyTransform([[0,1,2],[3,4,5]])
    print(orient.GetOrthoNormalBase())
    def CheckOperations(orient,p,v,t):
        pt = orient.ApplyTransform(p)
        pr = orient.ApplyInvTransform(pt)

        if np.linalg.norm(p-pr) > 1e-14:# pragma: no cover
            print(orient)
            print("Error ", np.linalg.norm(p-pr) )
            print("diff ", (p-pr) )
            print('p ',p)
            print('pt ',pt)
            print('pr ',pr)
            raise Exception("Error in the position transform")

        vt = orient.ApplyTransformDirection(v)
        vr = orient.ApplyInvTransformDirection(vt)

        if np.linalg.norm(v-vr) > 1e-14:# pragma: no cover
            print(orient)
            print("Error ", np.linalg.norm(v-vr) )
            print("diff ", (v-vr) )
            print('v ',v)
            print('vt ',vt)
            print('vr ',vr)
            raise Exception("Error in the direction transform")

        tt = orient.ApplyTransformTensor(t)
        tr = orient.ApplyInvTransformTensor(tt)

        if np.linalg.norm(t-tr) > 1e-14:# pragma: no cover
            print(orient)
            print("Error ", np.linalg.norm(t-tr) )
            print("diff ", (t-tr) )
            print('t ',t)
            print('tt ',tt)
            print('tr ',tr)
            raise Exception("Error in the tensor transform")

    for nor in (True,False):
        for ort in (True,False):
            orient = Transform()
            orient.keepNormalised = nor
            orient.keepOrthogonal = ort
            orient.SetFirst( [1., 0., 0.])
            orient.SetSecond([1., 1., 0.])
            orient.SetThird( [1., 1., 1.])

            p = [1., 3., 3.]
            v = [1., 4., 5.]
            t = np.array([[1., 2., 3.],[4., 5., 6.],[7., 8., 9.]])
            CheckOperations(orient,p,v,t)

    orient = Transform2D(offset=[1,1],
                        first=[1,0],
                        second =[1,1])

    assert(orient.EnsureIn2D(None) == None)
    np.testing.assert_equal(orient.EnsureIn2D([0,1,2]) , [0,1])

    assert(orient.EnsureIn3D(None) == None)
    np.testing.assert_equal(orient.EnsureIn3D([[0,1,2]]) , [[0,1,2]])

    orient.keepOrthogonal = True
    assert(orient.keepOrthogonal == True)

    orient.keepNormalised = True
    assert(orient.keepNormalised == True)

    orient.SetOffset([1,1])
    orient.SetOffset([1,1,0])

    orient.SetFirst([1,0])
    orient.SetFirst([1,0,0])



    print(orient.GetDirection(0) )
    #test with a list of points
    data = np.array([[0,1],[3,4]],dtype=PBasicFloatType)
    data_new = orient.ApplyTransform(data)
    data_new_neo = orient.ApplyInvTransform(data_new)
    np.testing.assert_allclose(data, data_new_neo)

    dir_data = np.array([[0,1]],dtype=PBasicFloatType)
    dir_data_new = orient.ApplyTransformDirection(dir_data)
    dir_data_new_neo = orient.ApplyInvTransformDirection(dir_data_new)
    np.testing.assert_allclose(dir_data, dir_data_new_neo)


    print(orient.GetOrthoNormalBase())

    print(orient )


    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity(GUI=True)))# pragma: no cover
