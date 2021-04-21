#distutils: language = c++
#cython: language_level = 3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

from libcpp.vector cimport vector

cimport numpy as np
import numpy as np
from scipy.sparse import coo_matrix


cdef extern from "LinAlg/EigenSolvers.h"  namespace "BasicTools"  :
    cdef cppclass EigenSolvers:
        EigenSolvers() except +
        void SetSolverType(int)
        void SetOp(int, int,int, double*,int*,int*, const double& ) nogil
        void Solve(int,double*,double*) nogil
        int GetSPQRRank()
        void GetSPQR_R(long*, long*, double*,long*,long*)
        void GetSPQR_Q(long*, long*, double*, long*, long*)
        int GetSPQR_R_nonZeros()
        int GetSPQR_Q_nonZeros()
        void GetSPQR_P(long*)

cdef class CEigenSolvers():
     cdef EigenSolvers cpp_object  # Hold a C++ instance which we're wrapping
     cdef int m
     cdef int n
     cdef double tol
     def __cinit__(self):
        self.m = 0
        self.n = 0
        self.tol = 1.e-6

     def SetTolerance(self,double tol):
         self.tol = tol

     def SetSolverType(self,str stype):
         if stype == "CG":
             self.cpp_object.SetSolverType(1)
         elif stype == "LU":
             self.cpp_object.SetSolverType(2)
         elif stype == "SPQR":
             self.cpp_object.SetSolverType(3)
         elif stype == "BiCGSTAB":
             self.cpp_object.SetSolverType(4)
         else:
             print("SolverType Not Available")
             raise

     def SetOp(self,matrix):
         data = coo_matrix(matrix)
         self.m = matrix.shape[0]
         self.n = matrix.shape[1]
         cdef np.ndarray[double, ndim=1, mode="c"] ddata =  data.data
         cdef np.ndarray[int, ndim=1, mode="c"] drow =  data.row
         cdef np.ndarray[int, ndim=1, mode="c"] dcol =  data.col
         self.SetOp_internal(self.m,self.n,data.data,data.row,data.col)


     def SetOp_internal(self,int m,int n,
               np.ndarray[double, ndim=1, mode="c"] data,
               np.ndarray[int, ndim=1, mode="c"] row,
               np.ndarray[int, ndim=1, mode="c"] col):

         cdef int data_shape = data.shape[0]
         cdef double* datap =  &data[0]
         cdef int* rowp =  &row[0]
         cdef int* colp =  &col[0]
         cdef double tol  = self.tol
         with nogil:
             self.cpp_object.SetOp(m,n,data_shape,datap,rowp,colp,tol)

     def Solve(self,np.ndarray[double, ndim=1, mode="c"] rhs,np.ndarray[double, ndim=1, mode="c"] sol ) :
         cdef int s = rhs.shape[0]
         cdef double* rhsp = &rhs[0]
         cdef double* solp = &sol[0]
         cdef EigenSolvers* solver = &self.cpp_object
         with nogil:
             solver.Solve(s, rhsp, solp)
         return sol

     def GetSPQRRank(self):
         cdef EigenSolvers* solver = &self.cpp_object
         return solver.GetSPQRRank()

     def GetSPQR_R(self):

         cdef EigenSolvers* solver = &self.cpp_object
         nzsize = solver.GetSPQR_R_nonZeros();
         cdef np.ndarray[double, ndim=1, mode="c"] vals = np.ndarray(nzsize,dtype=np.double);
         cdef np.ndarray[long, ndim=1, mode="c"] rows = np.ndarray(nzsize,dtype=int);
         cdef np.ndarray[long, ndim=1, mode="c"]cols = np.ndarray(nzsize,dtype=int);
         cdef np.ndarray[long, ndim=1, mode="c"] sizes = np.ndarray(2,dtype=int);
         from scipy.sparse import coo_matrix

         if nzsize == 0:
             return coo_matrix(([],([],[])), shape=(sizes[0],sizes[1] ))

         cdef long* si = &sizes[0]
         cdef long* sj = &sizes[1]
         cdef double* vp = &vals[0]
         cdef long* rp = &rows[0]
         cdef long* cp = &cols[0]

         solver.GetSPQR_R(si,sj, vp,rp,cp);
         data = (vals, (rows,cols))
         K = coo_matrix(data, shape=(sizes[0],sizes[1] ))
         return K

     def GetSPQR_P(self):
        cdef np.ndarray[long, ndim=1, mode="c"] p = np.ndarray(self.n,dtype=int);
        cdef long* pp = &p[0]
        cdef EigenSolvers* solver = &self.cpp_object
        solver.GetSPQR_P(pp);
        return p

     def GetSPQR_Q(self):

         cdef EigenSolvers* solver = &self.cpp_object
         nzsize = solver.GetSPQR_Q_nonZeros();
         cdef np.ndarray[double, ndim=1, mode="c"] vals = np.ndarray(nzsize,dtype=np.double);
         cdef np.ndarray[long, ndim=1, mode="c"] rows = np.ndarray(nzsize,dtype=int);
         cdef np.ndarray[long, ndim=1, mode="c"] cols = np.ndarray(nzsize,dtype=int);
         cdef np.ndarray[long, ndim=1, mode="c"] sizes = np.ndarray(2,dtype=int);

         from scipy.sparse import coo_matrix

         if nzsize == 0:
             return coo_matrix(([],([],[])), shape=(sizes[0],sizes[1] ))

         cdef long* si = &sizes[0]
         cdef long* sj = &sizes[1]
         cdef double* vp = &vals[0]
         cdef long* rp = &rows[0]
         cdef long* cp = &cols[0]

         solver.GetSPQR_Q(si,sj, vp,rp,cp);

         data = (vals, (rows,cols))
         Q = coo_matrix(data, shape=(sizes[0],sizes[1]) )
         return Q
