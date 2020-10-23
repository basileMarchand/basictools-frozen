# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


"""
    Use of bisect to search in sorted lists

    from:
    https://docs.python.org/3.8/library/bisect.html
"""

import bisect


def BinarySearch(ordered_list, item):
    """
    Searches in the sorted data "ordered_list" the rank of the largest element
    smaller than item, in log(len(ordered_list)) complexity

    Parameters
    ----------
    ordered_list: list or one-dimensional np.ndarray
        the data sorted in increasing order from which the previous rank is searched
    item : float or int
        the item for which the previous rank is searched

    Returns
    -------
    int
        the rank of the largest element smaller than item in the sorted data "list"
    """

    return max(bisect.bisect_right(ordered_list, item) - 1, 0)


def CheckIntegrity(GUI=False):
    import numpy as np

    testlist = np.array([0.0, 1.0, 2.5, 10.])
    valList = np.array([-1., 0., 11., 0.6, 2.0, 2.6, 9.9, 1.0])
    expected = (0, 0, 3, 0, 1, 2, 2, 1)
    actual = tuple(BinarySearch(testlist, v) for v in valList)

    assert actual == expected, "Test failed, actual = {}, expected = {}".format(actual, expected)

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
