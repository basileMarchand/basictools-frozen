# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np
from BasicTools.Containers.Filters import ElementFilter

def VerifyExclusiveFilters(listOfFilters,mesh):
    """Funtion to check if a list of ElementFilter is exclusive """

    for elemtype,data in mesh.elements.items():
        mask = np.zeros(data.GetNumberOfElements(),dtype=int)-1
        for fnb,f in enumerate(listOfFilters):
            ids = f.GetIdsToTreat(data)
            if np.any(mask[ids] > -1):
                idd = np.where(mask[ids] > -1)[0][0]
                raise(Exception(f"Filter {fnb} incompatiblewith filter {idd} "))
            mask[ids] = fnb

def ListOfElementFiltersFromETagList(taglist, checkExclusive=False, mesh=None):
    """Function to construct a list of filter from a list of tags (or a list of tags)

    taglist: is a list containing a string or a list of strings
    checkExclusive (default =False) to chech if the intersection of the filters
    is empty, in this case a mesh must be supplied

    return : It returns a list of ElementFilter with the associated tags as filters

    """

    listOfFilters = []
    for matname in taglist:
        if type(matname) in [ list, tuple ] :
            listOfFilters.append(ElementFilter(mesh=mesh,tags=matname))
        else:
            listOfFilters.append(ElementFilter(mesh=mesh,tag=matname))
    if checkExclusive:
        if mesh is None:
            raise Exception("Need a Mesh")
        else:
            VerifyExclusiveFilters(listOfFilters,mesh)

    return listOfFilters

#--------------  CheckIntegrity ---------------
def CheckIntegrity_ListOfElementFiltersFromETagList(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    mesh = CreateSquare()
    tags = ("2D",("X0","X1"),("Y0","Y1"))

    listOfFilters = ListOfElementFiltersFromETagList(tags)

    return "OK"


def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_ListOfElementFiltersFromETagList,
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
