gi# distutils: language = c++

from libcpp.vector cimport vector

cimport numpy as np
import numpy as np
from scipy.sparse import coo_matrix


cdef extern from "src_cpp/NativeEigenSolvers.h" :
    cdef cppclass NativeEigenSolvers:
        NativeEigenSolvers() except +
        void SetSolverType(int)
        void SetOp(int, vector[double],vector[int],vector[int])
        void Solve(int,double*,double*)

cdef class EigenSolvers():
     cdef NativeEigenSolvers c_solver  # Hold a C++ instance which we're wrapping
     def __cinit__(self):
        self.c_solver = NativeEigenSolvers()
     def SetSolverType(self,int stype):
         self.c_solver.SetSolverType(stype)
     def SetOp(self,matrix):
         data = coo_matrix(matrix)
         self.c_solver.SetOp(matrix.shape[0],data.data,data.row,data.col)
     def Solve(self,np.ndarray[double, ndim=1, mode="c"] rhs):
         cdef np.ndarray[double, ndim=1, mode="c"] sol = np.zeros_like(rhs)
         #print("inside pyx")
         #print(rhs)
         self.c_solver.Solve(rhs.shape[0], &rhs[0], &sol[0])
         return sol

