//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <iostream>
#include <memory>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/Core>

#include <LinAlg/EigenTypes.h>


namespace BasicTools {
    
using Eigen::Success;
using Eigen::Upper;
using Eigen::Lower;

//# solverType 0 CG

struct proxy;
struct container;
struct iterator;

struct proxy {
    const container* cont;
    long int c_pos;

    proxy(const container* cont, long int c_pos) : cont(cont), c_pos(c_pos) {}

    proxy* operator->() { return this; }
    const proxy* operator->() const { return this; }

    const double& value() const; // the value
    INT_TYPE row() const;   // the row index i
    INT_TYPE col() const;   // the column index j
 };

class iterator {
public:
    typedef iterator self_type;

    iterator(const container* cont, long int c_pos) : cont(cont), c_pos(c_pos) {}

    const self_type& operator++() { ++c_pos; return *this; }
    self_type operator++(int) { const auto i = *this; ++c_pos; return i; }

    bool operator==(const self_type& rhs) const { return cont == rhs.cont && c_pos == rhs.c_pos; }
    bool operator!=(const self_type& rhs) const { return !(*this == rhs); }

    proxy operator->() const { return proxy(cont, c_pos); }

private:
    const container* cont;
    long int c_pos;
};

struct container {
    double* ev;
    int* ei;
    int* ej;
    long int size;

    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, size); }
    const double& value(long int cpt) const { return ev[cpt]; }; // the value
    double& value(long int cpt) { return ev[cpt]; }; // the value
    INT_TYPE row(long int cpt) const { return ei[cpt]; };   // the row index i
    INT_TYPE col(long int cpt) const { return ej[cpt]; };   // the column index j
};

const double& proxy::value() const { return cont->value(c_pos); } // the value
INT_TYPE proxy::row() const { return cont->row(c_pos); }   // the row index i
INT_TYPE proxy::col() const { return cont->col(c_pos); }   // the column index j

template<typename T>
void CopyMatrix(INT_TYPE* sizei, INT_TYPE* sizej, FLOAT_TYPE* ev, INT_TYPE* ei, INT_TYPE* ej, T& mat) {
    *sizei = mat.rows();
    *sizej = mat.cols();
    int cpt = 0;
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
            ev[cpt] = it.value();
            ei[cpt] = it.row();   // row index
            ej[cpt] = it.col();   // col index (here it is equal to k)
            ++cpt;
        }
    }
}

// for debugin
template<typename T>
void printMatrix(const std::string& name, T& mat) {
    Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    std::string sep = "\n----------------------------------------\n";
    std::cout << name << sep;
    std::cout << mat << sep;
}

struct NativeEigenSolvers {
    typedef Eigen::ConjugateGradient<SpMatD, Upper|Lower> EigenSpCG;
    typedef Eigen::SparseLU<SpMatD> EigenSpLU;
    typedef Eigen::SparseQR<SpMatD, Eigen::COLAMDOrdering<int>> EigenSpQR;
    typedef Eigen::BiCGSTAB<SpMatDR> EigenSpBiCGSTAB;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    std::unique_ptr<EigenSpCG> cgSolver;
    std::unique_ptr<EigenSpLU> luSolver;
    std::unique_ptr<EigenSpQR> spqrSolver;
    std::unique_ptr<EigenSpBiCGSTAB> bicgstabSolver;

    std::unique_ptr<SpMatD> A;
    std::unique_ptr<SpMatDR> A_rowmajor;
    int solverType;
    Eigen::SparseMatrix<double> Q;
    Eigen::SparseMatrix<double> R;

    NativeEigenSolvers() : solverType(0) {
        // Nothing to do
    }

    void Clean() {
        this->cgSolver = nullptr;
        this->luSolver = nullptr;
        this->spqrSolver = nullptr;
        this->bicgstabSolver = nullptr;
        this->A = nullptr;
        this->A_rowmajor = nullptr;
        this->Q.resize(0, 0);
        this->R.resize(0, 0);

        this->solverType = 0;
    }

    void SetSolverType(int i) {
        this->Clean();
        this->solverType = i;
        if (i == 1) {
            this->cgSolver.reset(new EigenSpCG);
        } else if (i == 2) {
            this->luSolver.reset(new EigenSpLU);
        } else if (i == 3) {
            this->spqrSolver.reset(new EigenSpQR);
        } else if (i == 4) {
            this->bicgstabSolver.reset(new EigenSpBiCGSTAB);
        } else {
            std::cout << "Solver type " << i << " not avilable " << std::endl;
        }
    }

    void SetOp(const int& sizem, const int& sizen, const int& ev_size, FLOAT_TYPE* ev, int* ei, int* ej, const double& tolerance) {
        container cont;
        cont.ev = ev;
        cont.ei = ei;
        cont.ej = ej;
        cont.size = ev_size;

        if (this->solverType == 4) {
            this->A_rowmajor.reset(new SpMatDR(sizem, sizen));
            this->A_rowmajor->setFromTriplets(cont.begin(), cont.end());
            this->A_rowmajor->makeCompressed();
        } else {
            this->A.reset(new SpMatD(sizem, sizen));
            this->A->setFromTriplets(cont.begin(), cont.end());
            this->A->makeCompressed();
        }

        if (this->solverType == 0) {
            std::cout << "ERROR! NativeEigenSolvers::SetOp() Solver type not set " << std::endl;
            exit(0);
        }

        if (this->solverType == 1) {
            this->cgSolver->setTolerance(tolerance);
            this->cgSolver->compute(*A);
            if (this->cgSolver->info() != Success) {
                std::cout << "CG compute failed" << std::endl;
                // decomposition failed
                return;
            }
        } else if (this->solverType == 2) {
            this->luSolver->analyzePattern(*A);
            // Compute the numerical factorization
            this->luSolver->factorize(*A);
            if (this->luSolver->info() != Success) {
                std::cout << "lu decomposition failed" << std::endl;
                // decomposition failed
                return;
            }
        } else if (this->solverType == 3) {
            this->spqrSolver->setPivotThreshold(tolerance);
            this->spqrSolver->compute(*A);
            if (this->spqrSolver->info() != Success) {
                std::cout << "spqr decomposition failed" << std::endl;
                // decomposition failed
                return;
            }
        } else if (this->solverType == 4) {
            this->bicgstabSolver->setTolerance(tolerance);
            this->bicgstabSolver->compute(*A_rowmajor);
            if (this->bicgstabSolver->info() != Success) {
                std::cout << "bicgstab compute failed" << std::endl;
                // decomposition failed
                return;
            }
        }
    }

    void Solve(int size, FLOAT_TYPE* _rhs, FLOAT_TYPE* _sol) {
        MapMatrixDD1 rhs(_rhs, size, 1);
        MapMatrixDD1 sol(_sol, size, 1);

        if (this->solverType == 1) {
            sol = this->cgSolver->solveWithGuess(rhs, sol);
        } else if (this->solverType == 2) {
            sol = this->luSolver->solve(rhs);
        } else if (this->solverType == 3) {
            sol = this->spqrSolver->solve(rhs);
        } else if (this->solverType == 4) {
            sol = this->bicgstabSolver->solveWithGuess(rhs, sol);
        }
    }

    int GetSPQRRank() {
        if (this->solverType != 3) {
            std::cout << "spqr decomposition not avilable Please Set solver to EigenSPQR " << std::endl;
            return 0;
        }
        return this->spqrSolver->rank();
    }

    int GetSPQR_R_nonZeros() {
        if (this->solverType != 3) {
            std::cout << "spqr decomposition not avilable Please Set solver to EigenSPQR" << std::endl;
            return 0;
        }
        this->R = this->spqrSolver->matrixR();
        return this->R.nonZeros();
    }

    void GetSPQR_R(INT_TYPE* sizei, INT_TYPE* sizej, FLOAT_TYPE* ev, INT_TYPE* ei, INT_TYPE* ej) {
        if (this->solverType != 3) {
            std::cout << "spqr decomposition not avilable Please Set solver to EigenSPQR" << std::endl;
        }
        CopyMatrix(sizei, sizej, ev, ei, ej, this->R);
    }

    void GetSPQR_P(INT_TYPE* p) {
        if (this->solverType != 3) {
            std::cout << "spqr decomposition not avilable Please Set solver to EigenSPQR" << std::endl;
        }
        const auto indices = this->spqrSolver->colsPermutation().indices();
        for (int i = 0 ; i < this->A->cols(); ++i) {
            p[i] = indices[i];
        }
    }

    int GetSPQR_Q_nonZeros() {
        if (this->solverType != 3) {
            std::cout << "spqr decomposition not avilable Please Set solver to EigenSPQR" << std::endl;
            return 0;
        }

        this->Q = this->spqrSolver->matrixQ();
        return this->Q.nonZeros();
    }

    void GetSPQR_Q(INT_TYPE* sizei, INT_TYPE* sizej, FLOAT_TYPE* ev, INT_TYPE* ei, INT_TYPE* ej) {
        if (this->solverType != 3) {
            std::cout << "spqr decomposition not avilable Please Set solver to EigenSPQR" << std::endl;
            return;
        }
        CopyMatrix(sizei, sizej, ev, ei, ej, this->Q);
    }

};

} // namespace BasicTools
