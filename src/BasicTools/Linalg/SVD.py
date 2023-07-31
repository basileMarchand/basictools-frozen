# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np


def TruncatedSVDSymLower(matrix, epsilon = None, nbModes = None):
    """
    Computes a truncatd singular value decomposition of a symetric definite
    matrix in scipy.sparse.csr format. Only the lower triangular part needs
    to be defined

    Parameters
    ----------
    matrix : scipy.sparse.csr
        the input matrix
    epsilon : float
        the truncation tolerence, determining the number of keps eigenvalues
    nbModes : int
        the number of keps eigenvalues

    Returns
    -------
    np.ndarray
        kept eigenvalues, of size (numberOfEigenvalues)
    np.ndarray
        kept eigenvectors, of size (numberOfEigenvalues, numberOfSnapshots)
    """

    if epsilon != None and nbModes != None:# pragma: no cover
        raise("cannot specify both epsilon and nbModes")

    eigenValues, eigenVectors = np.linalg.eigh(matrix, UPLO="L")

    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:, idx]

    if nbModes == None:
        if epsilon == None:
            nbModes  = matrix.shape[0]
        else:
            nbModes = 0
            bound = (epsilon ** 2) * eigenValues[0]
            for e in eigenValues:
                if e > bound:
                    nbModes += 1
            id_max2 = 0
            bound = (1 - epsilon ** 2) * np.sum(eigenValues)
            temp = 0
            for e in eigenValues:
                temp += e
                if temp < bound:
                    id_max2 += 1  # pragma: no cover

            nbModes = max(nbModes, id_max2)

    if nbModes > matrix.shape[0]:
        print("nbModes taken to max possible value of "+str(matrix.shape[0])+" instead of provided value "+str(nbModes))
        nbModes = matrix.shape[0]

    index = np.where(eigenValues<0)
    if len(eigenValues[index])>0:
        if index[0][0]<nbModes:
            #print(nbModes, index[0][0])
            print("removing numerical noise from eigenvalues, nbModes is set to "+str(index[0][0])+" instead of "+str(nbModes))
            nbModes = index[0][0]

    return eigenValues[0:nbModes], eigenVectors[:, 0:nbModes]


def CheckIntegrity(GUI=False):

    Mat = np.arange(9).reshape(3, 3)
    Mat[np.triu_indices(3, 1)] = 0.0

    ref = (np.array([16.01240515]), np.array([[-0.38252793],[-0.53464575],[-0.7535425 ]]))


    res = TruncatedSVDSymLower(Mat)
    np.testing.assert_almost_equal(res[0], ref[0])
    np.testing.assert_almost_equal(res[1], ref[1])


    res = TruncatedSVDSymLower(Mat, 1.0e-6)
    np.testing.assert_almost_equal(res[0], ref[0])
    np.testing.assert_almost_equal(res[1], ref[1])


    res = TruncatedSVDSymLower(Mat, nbModes = 5)
    np.testing.assert_almost_equal(res[0], ref[0])
    np.testing.assert_almost_equal(res[1], ref[1])


    res = TruncatedSVDSymLower(Mat, nbModes = 11)
    np.testing.assert_almost_equal(res[0], ref[0])
    np.testing.assert_almost_equal(res[1], ref[1])

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
