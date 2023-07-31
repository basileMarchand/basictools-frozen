# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE', which is part of this source code package.
#

import numpy as np
import bisect


def BinarySearch(orderedList, item):
    """
    Inspects the sorted list "ordered_list" and returns:
        - 0 if item <= orderedList[0]
        - the rank of the largest element smaller or equal than item otherwise

    Parameters
    ----------
    ordered_list: list or one-dimensional np.ndarray
        the data sorted in increasing order from which the previous rank is\
        searched
    item : float or int
        the item for which the previous rank is searched

    Returns
    -------
    int
        0 or the rank of the largest element smaller or equal than item in\
        "orderedList"
    """
    return max(bisect.bisect_right(orderedList, item) - 1, 0)



def BinarySearchVectorized(orderedList, items):
    """
    BinarySearch for more than one call (items is now a list or one-dimensional np.ndarray)
    """
    return np.fromiter(map(lambda item: BinarySearch(orderedList, item), items), dtype = int)



def PieceWiseLinearInterpolation(item, itemIndices, vectors, tolerance = 1e-4):
    """
    Computes a item interpolation for temporal vectors defined either by
    itemIndices  and vectors at these indices

    Parameters
    ----------
    item : float
        the input item at which the interpolation is required
    itemIndices : np.ndarray
        the items where the available data is defined, of size
        (numberOfTimeIndices)
    vectors : np.ndarray or dict
        the available data, of size (numberOfVectors, numberOfDofs)
    tolerance : float
        tolerance for deciding when using the closest timestep value instead of carrying out the linear interpolation

    Returns
    -------
    np.ndarray
        interpolated vector, of size (numberOfDofs)
    """

    if item <= itemIndices[0]:
        return vectors[0]
    if item >= itemIndices[-1]:
        return vectors[-1]


    prev = BinarySearch(itemIndices, item)
    coef = (item - itemIndices[prev]) / (itemIndices[prev+1] - itemIndices[prev])

    if 0.5 - abs(coef - 0.5) < tolerance:
        coef = round(coef)

    return (
            coef * vectors[prev+1]
            + (1 - coef) * vectors[prev]
        )


def PieceWiseLinearInterpolationWithMap(item, itemIndices, vectors, vectorsMap, tolerance = 1e-4):
    """
    Computes a item interpolation for temporal vectors defined either by
    itemIndices, some tags at these item indices (vectorsMap), and vectors at those tags.

    Parameters
    ----------
    item : float
        the input item at which the interpolation is required
    itemIndices : np.ndarray
        the items where the available data is defined, of size
        (numberOfTimeIndices)
    vectors : np.ndarray or dict
        the available data, of size (numberOfVectors, numberOfDofs)
    vectorsMap : list
        list containing the mapping from the numberOfTimeIndices items indices to the numberOfVectors vectors, of size (numberOfTimeIndices,). Default is None, in which case numberOfVectors = numberOfTimeIndices.
    tolerance : float
        tolerance for deciding when using the closest timestep value instead of carrying out the linear interpolation

    Returns
    -------
    np.ndarray
        interpolated vector, of size (numberOfDofs)
    """


    if item <= itemIndices[0]:
        return vectors[vectorsMap[0]]
    if item >= itemIndices[-1]:
        return vectors[vectorsMap[-1]]


    prev = BinarySearch(itemIndices, item)
    coef = (item - itemIndices[prev]) / (itemIndices[prev+1] - itemIndices[prev])

    if 0.5 - abs(coef - 0.5) < tolerance:
        coef = round(coef)

    return (
            coef * vectors[vectorsMap[prev+1]]
            + (1 - coef) * vectors[vectorsMap[prev]]
        )



def PieceWiseLinearInterpolationVectorized(items, itemIndices, vectors):
    """
    PieceWiseLinearInterpolation for more than one call (items is now a list or one-dimensional np.ndarray)
    """
    return [PieceWiseLinearInterpolation(item, itemIndices, vectors) for item in items]
    #return np.fromiter(map(lambda item: PieceWiseLinearInterpolation(item, itemIndices, vectors), items), dtype = type(vectors[0]))



def PieceWiseLinearInterpolationVectorizedWithMap(items, itemIndices, vectors, vectorsMap):
    """
    PieceWiseLinearInterpolation for more than one call (items is now a list or one-dimensional np.ndarray)
    """
    return [PieceWiseLinearInterpolationWithMap(item, itemIndices, vectors, vectorsMap) for item in items]
    #return np.fromiter(map(lambda item: PieceWiseLinearInterpolationWithMap(item, itemIndices, vectors, vectorsMap), items), dtype = np.float)



def CheckIntegrity(GUI=False):

    timeIndices = np.array([0.0, 1.0, 2.5])
    vectors = np.array([np.ones(5), 2.0 * np.ones(5), 3.0 * np.ones(5)])
    vectorsDic = {}
    vectorsDic["vec1"] = np.ones(5)
    vectorsDic["vec2"] = 2.0 * np.ones(5)
    vectorsMap = ["vec1", "vec2", "vec1"]

    res = PieceWiseLinearInterpolation(-1.0, timeIndices, vectors)
    np.testing.assert_almost_equal(res, [1., 1., 1., 1., 1.])

    res = PieceWiseLinearInterpolationWithMap(3.0, timeIndices, vectorsDic, vectorsMap)
    np.testing.assert_almost_equal(res, [1., 1., 1., 1., 1.])

    res = PieceWiseLinearInterpolation(1.0, timeIndices, vectors)
    np.testing.assert_almost_equal(res, [2., 2., 2., 2., 2.])

    res = PieceWiseLinearInterpolationWithMap(1.0, timeIndices, vectorsDic, vectorsMap)
    np.testing.assert_almost_equal(res, [2., 2., 2., 2., 2.])

    res = PieceWiseLinearInterpolation(0.4, timeIndices, vectors)
    np.testing.assert_almost_equal(res, [1.4, 1.4, 1.4, 1.4, 1.4])

    res = PieceWiseLinearInterpolation(1.4, timeIndices, vectors)
    np.testing.assert_almost_equal(res, [6.8/3, 6.8/3, 6.8/3, 6.8/3, 6.8/3])

    res = PieceWiseLinearInterpolationWithMap(0.6, timeIndices, vectorsDic, vectorsMap)
    np.testing.assert_almost_equal(res, [1.6, 1.6, 1.6, 1.6, 1.6])

    res = PieceWiseLinearInterpolationVectorizedWithMap(np.array([-0.1, 2.0, 3.0]), timeIndices, vectorsDic, vectorsMap)
    np.testing.assert_almost_equal(res, [np.array([1., 1., 1., 1., 1.]), np.array([1.33333333, 1.33333333, 1.33333333, 1.33333333, 1.33333333]), np.array([1., 1., 1., 1., 1.])])


    timeIndices = np.array([0., 100., 200.,  300.,  400.,  500.,  600.,  700.,\
                            800.,  900., 1000., 2000.])
    coefficients = np.array([2000000., 2200000., 2400000., 2000000., 2400000.,\
                            3000000., 2500000., 2400000., 2100000., 2800000.,\
                            4000000., 3000000.])


    vals = np.array([-10., 0., 100., 150., 200., 300., 400., 500., 600., 700.,\
                    800., 900., 1000., 3000., 701.4752695491923])


    res = np.array([2000000., 2000000., 2200000., 2300000., 2400000., 2000000.,\
                    2400000., 3000000., 2500000., 2400000., 2100000., 2800000.,\
                    4000000., 3000000., 2395574.19135242])

    for i in range(vals.shape[0]):
        assert (PieceWiseLinearInterpolation(vals[i], timeIndices, coefficients) - res[i])/res[i] < 1.e-10


    res = PieceWiseLinearInterpolationVectorized(np.array(vals), timeIndices, coefficients)
    np.testing.assert_almost_equal(res, [2000000.0, 2000000.0, 2200000.0, 2300000.0, 2400000.0, 2000000.0, 2400000.0, 3000000.0, 2500000.0, 2400000.0, 2100000.0, 2800000.0, 4000000.0, 3000000.0, 2395574.1913524233])

    testlist = np.array([0.0, 1.0, 2.5, 10.])
    valList = np.array([-1., 11., 0.6, 2.0, 2.6, 9.9, 1.0])

    ref = np.array([0, 3, 0, 1, 2, 2, 1], dtype = int)
    res = BinarySearchVectorized(testlist, valList)

    for i, val in enumerate(valList):
        assert BinarySearch(testlist, val) == ref[i]
        assert res[i] == ref[i]

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
