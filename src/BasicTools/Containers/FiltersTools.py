# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np
from BasicTools.Containers.Filters import ElementFilter

def VerifyExclusiveFilters(listOfFilters,mesh):
    """Funtion to check if a list of ElementFilter is exclusive each element is present only one
    listOfFilters : the list of filters to check
    mesh : the mesh to apply the filters

    raise exception is the case the intersection is not empty
    """

    for elemtype,data in mesh.elements.items():
        mask = np.zeros(data.GetNumberOfElements(),dtype=int)-1
        for fnb,f in enumerate(listOfFilters):
            ids = f.GetIdsToTreat(data)
            if np.any(mask[ids] > -1):
                idd = np.where(mask[ids] > -1)[0][0]
                raise(Exception(f"Filter {fnb} incompatiblewith filter {idd} "))
            mask[ids] = fnb

def ListOfElementFiltersFromETagList(taglist, mesh=None):
    """Function to construct a list of filter from a list of tags (or a list of tags)

    taglist: is a list containing a string or a list of strings ("tag1",("tag2","tag3"))
    mesh: mesh to pass to the filters

    return : It returns a list of ElementFilter with the associated tags as filters

    """

    listOfFilters = []
    for matname in taglist:
        if type(matname) in [ list, tuple ] :
            listOfFilters.append(ElementFilter(mesh=mesh,tags=matname))
        else:
            listOfFilters.append(ElementFilter(mesh=mesh,tag=matname))

    return listOfFilters

def ListOfElementFiltersFromMask(maskVector, mesh=None):
    """Function to construct a list of filter from a mask vector

    maskVector : a element size vector of ints to determine to which filter each
                elements belongs to. The number of filters is calculated unsing the np.unique

    mesh: mesh to pass to the filters

    return : It returns a list of ElementFilter with the associated tags as filters

    """
    ids = np.unique(maskVector.flatten())

    listOfFilters = []
    for partition_id in ids:
        mask  = maskVector == partition_id
        listOfFilters.append(ElementFilter(mesh=mesh,mask=mask))

    return listOfFilters

#--------------  CheckIntegrity ---------------
def CheckIntegrity_ListOfElementFiltersFromMask(GUI=False):

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByElementFilter


    mesh = CreateSquare(dimensions=[10,10])
    mask = np.asarray(np.random.random(mesh.GetNumberOfElements())*5,dtype=int)

    listOfFilters = ListOfElementFiltersFromMask(mask,mesh)
    VerifyExclusiveFilters(listOfFilters,mesh)

    tmp = TestTempDir.GetTempPath()

    nbelements = 0
    for i,ff in enumerate(listOfFilters):
        new_mesh = ExtractElementsByElementFilter(mesh,ff)
        nbelements += new_mesh.GetNumberOfElements()
        print(f"Number of element in the mesh {i}:" + str(new_mesh.GetNumberOfElements()))
        WriteMeshToXdmf(tmp + f"CI_ListOfElementFiltersFromMask{i}.xdmf", new_mesh)

    WriteMeshToXdmf(tmp + "CI_ListOfElementFiltersFromMask_base.xdmf", mesh,CellFields=[mask],CellFieldsNames=["RANK"])

    print((nbelements,mesh.GetNumberOfElements()))

    if nbelements != mesh.GetNumberOfElements():
        return "Not OK"

    return "OK"

def CheckIntegrity_ListOfElementFiltersFromETagList(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    mesh = CreateSquare()
    tags = ("2D",("X0","X1"),("Y0","Y1"))

    listOfFilters = ListOfElementFiltersFromETagList(tags)

    return "OK"


def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_ListOfElementFiltersFromETagList,
    CheckIntegrity_ListOfElementFiltersFromMask
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
