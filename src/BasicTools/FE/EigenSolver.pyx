# distutils: language = c++

from libcpp.vector cimport vector

cimport numpy as np
import numpy as np
from scipy.sparse import coo_matrix


cdef extern from "src_cpp/NativeEigenSolvers.h" :
    cdef cppclass NativeEigenSolvers:
        NativeEigenSolvers() except +
        void SetSolverType(int)
        void SetOp(int, int, double*,int*,int*) nogil
        void Solve(int,double*,double*) nogil

cdef class EigenSolvers():
     cdef NativeEigenSolvers c_solver  # Hold a C++ instance which we're wrapping
     def __cinit__(self):
        self.c_solver = NativeEigenSolvers()
     def SetSolverType(self,int stype):
         self.c_solver.SetSolverType(stype)
     def SetOp(self,matrix):
         data = coo_matrix(matrix)
         cdef int shape  = matrix.shape[0]
         cdef np.ndarray[double, ndim=1, mode="c"] ddata =  data.data
         cdef np.ndarray[int, ndim=1, mode="c"] drow =  data.row
         cdef np.ndarray[int, ndim=1, mode="c"] dcol =  data.col
         self.SetOp_internal(shape,data.data,data.row,data.col)

     def SetOp_internal(self,int shape,
                np.ndarray[double, ndim=1, mode="c"] data,
               np.ndarray[int, ndim=1, mode="c"] row,
               np.ndarray[int, ndim=1, mode="c"] col):

         cdef int data_shape = data.shape[0]
         cdef double* datap =  &data[0]
         cdef int* rowp =  &row[0]
         cdef int* colp =  &col[0]

         with nogil:
             self.c_solver.SetOp(shape,data_shape,datap,rowp,colp)
     def Solve(self,np.ndarray[double, ndim=1, mode="c"] rhs) :
         cdef np.ndarray[double, ndim=1, mode="c"] sol = np.zeros_like(rhs)
         cdef int s = rhs.shape[0]
         cdef double* rhsp = &rhs[0]
         cdef double* solp = &sol[0]
         cdef NativeEigenSolvers* solver = &self.c_solver
         with nogil:
             solver.Solve(s, rhsp, solp)
         return sol

