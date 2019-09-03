# -*- coding: utf-8 -*-

import numpy as np

# to eliminate dofs with know values

def deleterowcol(A, delrow, delcol, fixedValues ):
    # Assumes that matrix is in symmetric csc form !

    rhs = A.dot(fixedValues)
    #keep = np.delete (np.arange(0, m), delrow)
    A = A[np.logical_not(delrow) , :]
    #keep = np.delete (np.arange(0, m), delcol)
    A = A[:, np.logical_not(delcol)]

    return [A, rhs]


def CheckIntegrity():
    from scipy.sparse import csr_matrix


    fv = np.array([[1,2],]).T
    mask = np.zeros(2)
    mask[1] = True
    for sp in [True,False]:
        if sp:
            K = csr_matrix([[1, 2], [3, 4]])
        else:
            K = np.array([[1, 2], [3, 4]]);

        A,rhs = deleterowcol(K, mask, mask, fv )
        print("using Sparse : " + ("True" if sp else "False" ) )
        print("Vals")
        print(A)
        print("--")
        print(rhs)
        print("Types")
        print(type(A))
        print(type(rhs))
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
