//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <iostream>



#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <Linalg/src_cpp/Eigentype.h>

using Eigen::Success;
using Eigen::Upper;
using Eigen::Lower;

//# solverType 0 CG

struct proxy;
struct container;
struct iterator;

struct proxy {
   proxy(const container* cont,long int c_pos ):cont(cont),c_pos (c_pos){}
   proxy* operator->() { return this; }
   const proxy* operator->() const { return this; }
   const container* cont;
   long int c_pos;
   double& value() const; // the value
   INT_TYPE row() const;   // the row index i
   INT_TYPE col() const;   // the column index j
 };

class iterator{
    public:
     iterator(const container* cont,long int c_pos ):cont(cont),c_pos (c_pos){}
     typedef iterator self_type;

        //reference operator*() { return *ptr_; }

        self_type operator++(int junk) { c_pos++; return *this; }
        self_type operator++() { self_type i = *this; c_pos++; return i; }

        bool operator!=(const self_type& rhs) const  { return c_pos != rhs.c_pos; }
        bool operator==(const self_type& rhs) const  { return c_pos == rhs.c_pos; }

        proxy operator->() { return proxy(cont,c_pos); }

    private:
        const container* cont;
        long int c_pos;
};

struct container {
        iterator begin()const {return  iterator(this,0); }
        iterator end()const {return  iterator(this,size); }
        double& value(long int cpt) const {return ev[cpt];}; // the value
        INT_TYPE row(long int cpt) const {return ei[cpt];};   // the row index i
        INT_TYPE col(long int cpt) const {return ej[cpt];};   // the column index j
        double* ev;
        int* ei;
        int* ej;
        long int size;
};
double& proxy::value() const {return cont->value(c_pos);}; // the value
INT_TYPE proxy::row() const {return cont->row(c_pos);};   // the row index i
INT_TYPE proxy::col() const {return cont->col(c_pos);};   // the column index j


struct NativeEigenSolvers {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::ConjugateGradient<SpMatD, Upper|Lower >* cgSolver;
    Eigen::SparseLU<SpMatD >*   luSolver;

    SpMatD* A;
    int solverType;

    NativeEigenSolvers(){
       this->solverType = 0;
       this->cgSolver = 0;
       this->luSolver = 0;
       this->A = 0;
    }
    ~NativeEigenSolvers(){
        this->Clean();
    }
    void Clean(){
       if(this->cgSolver) delete cgSolver;
       if(this->luSolver) delete luSolver;
       if(this->A) delete A;
    }
    void SetSolverType(int i ){
        this->Clean();
        this->solverType = i;
        if(i == 1){
            this->cgSolver = new Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Upper|Lower >();
        } else if(i == 2) {
            this->luSolver = new Eigen::SparseLU<SpMatD>();
        } else{
            std::cout << "Solver type " << i << " not avilable " << std::endl;
        }
    };
    void SetOp(const int& size,const int& ev_size ,FLOAT_TYPE* ev, int* ei, int* ej){
        this->A = new SpMatD(size,size);
        container cont;
        cont.ev = ev;
        cont.ei = ei;
        cont.ej = ej;
        cont.size = ev_size;
        this->A->setFromTriplets(cont.begin(), cont.end());
        this->A->makeCompressed ();

        if(this->solverType == 0 ) exit(0);

        if(this->solverType == 1 ){
            this->cgSolver->setTolerance(1.E-6);
            this->cgSolver->compute(*A);
            if(this->cgSolver->info()!=Success) {
                std::cout << "CG compute failed" << std::endl;
                // decomposition failed
                return;
            }
        } else if(this->solverType == 2){
            this->luSolver->analyzePattern(*A);
            // Compute the numerical factorization
            this->luSolver->factorize(*A);
            if(this->luSolver->info()!=Success) {
                std::cout << "lu decomposition failed" << std::endl;
                // decomposition failed
                return;
            }
        }
    };

    void Solve(int size, FLOAT_TYPE*  _rhs, FLOAT_TYPE*_sol) {
        MapMatrixDD1* rhs = new MapMatrixDD1(_rhs,size,1);
        MapMatrixDD1* sol = new MapMatrixDD1(_sol,size,1);

        if(this->solverType == 1 ){
            (*sol) = this->cgSolver->solve(*(rhs));
        } else if (this->solverType == 2 ){
            (*sol) = this->luSolver->solve(*(rhs));
        }

        delete rhs;
        delete sol;
    };
};
