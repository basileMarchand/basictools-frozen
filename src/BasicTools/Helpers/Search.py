# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
        

def BinarySearch(list, item):
    """
    Searches in the sorted data "list" the rank of the largest element smaller than item, in log(len(list)) complexity

    Parameters
    ----------
    list: list or one-dimensional np.ndarray
        the sorted data from which the previous rank is searched
    item : float or int
        the item for which the previous rank is searched

    Returns
    -------
    int
        the rank of the largest element smaller than item in the sorted data "list"
    """
    first = 0
    last = len(list)-1
    res = None

    if item < list[first]:
        return first

    if item > list[last]:
        return last

    while first<=last:
        midpoint = (first + last)//2
        if list[midpoint] == item:
            return midpoint
        else:
            if item <= list[midpoint]:
                last = midpoint-1
                res = midpoint-1
            else:
                first = midpoint+1
                res = midpoint

    return res


def CheckIntegrity(GUI=False):


    import numpy as np
    
    testlist = np.array([0.0, 1.0, 2.5, 10.])
    valList = np.array([-1., 11., 0.6, 2.0, 2.6, 9.9, 1.0])

    for i, val in enumerate(valList):
        BinarySearch(testlist, val)
    

    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))