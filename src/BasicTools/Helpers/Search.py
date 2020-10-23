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

    first = 0
    last = len(ordered_list) - 1

    if item < ordered_list[first]:
        return first

    if item > ordered_list[last]:
        return last

    return index(ordered_list, find_gt(ordered_list, item)) - 1


def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect.bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect.bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_gt(a, x):
    'Find leftmost value greater than x'
    i = bisect.bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError


def CheckIntegrity(GUI=False):


    import numpy as np

    testlist = np.array([0.0, 1.0, 2.5, 10.])
    valList = np.array([-1., 0., 11., 0.6, 2.0, 2.6, 9.9, 1.0])


    for i, val in enumerate(valList):
        BinarySearch(testlist, val)

    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))
