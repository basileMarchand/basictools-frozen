# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from typing import List, Union, Optional
from numpy.typing import ArrayLike

from BasicTools.Containers.Filters import ElementFilter, FilterOP
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh

def VerifyExclusiveFilters(listOfFilters: List[Union[ElementFilter,FilterOP]], mesh:UnstructuredMesh) -> bool :
    """Funtion to check if a list of ElementFilter is exclusive.
       (each element is present at the most in one filter)

    Parameters
    ----------
    listOfFilters : List[Union[ElementFilter,FilterOP]]
        The list of filters to check
    mesh : UnstructuredMesh
        Mesh to evaluate the filters

    Returns
    -------
    bool
        True if the filters are exclusive
    """

    for elemtype,data in mesh.elements.items():
        mask = np.zeros(data.GetNumberOfElements(),dtype=int)-1
        for fnb,f in enumerate(listOfFilters):
            ids = f.GetIdsToTreat(data)
            if np.any(mask[ids] > -1):
                idd = np.where(mask[ids] > -1)[0][0]
                print(f"Filter {fnb} incompatible with filter {idd} ")
                return 0
            mask[ids] = fnb

    return True

def ListOfElementFiltersFromETagList(taglist: List[Union[str,List[str]]], mesh: Optional[UnstructuredMesh] = None) -> List[ElementFilter]:
    """Function to construct a list of filters from a list of tags (or a list of tags)

    Parameters
    ----------
    taglist : List[Union[str,List[str]]]
        A list containing a string or a list of strings ("tag1",("tag2","tag3"))
    mesh : Optional[UnstructuredMesh], optional
        Mesh to pass to the filters

    Returns
    -------
    List[ElementFilter]
        List of ElementFilter with the associated tags as filters
    """

    listOfFilters = []
    for matname in taglist:
        if type(matname) in [ list, tuple ] :
            listOfFilters.append(ElementFilter(mesh=mesh,tags=matname))
        else:
            listOfFilters.append(ElementFilter(mesh=mesh,tag=matname))

    return listOfFilters


def ListOfElementFiltersFromMask(maskVector: ArrayLike, mesh: Optional[UnstructuredMesh]=None) -> List[ElementFilter]:
    """Function to construct a list of filter from a mask vector

    Parameters
    ----------
    maskVector : ArrayLike
        A element size vector of ints to determine to which filter each elements belongs to.
        The number of filters is calculated unsing the np.unique

    mesh : Optional[UnstructuredMesh], optional
        Mesh to pass to the filters

    Returns
    -------
    List[ElementFilter]
        List of ElementFilter with the associated tags as filters

    """
    ids = np.unique(maskVector.flatten())

    listOfFilters = []
    for partition_id in ids:
        mask  = maskVector == partition_id
        listOfFilters.append(ElementFilter(mesh=mesh,mask=mask))

    return listOfFilters

def FilterToETag(mesh: UnstructuredMesh, elementFilter: Union[ElementFilter,FilterOP], tagname: str) -> None:
    """Create a Element tag with the name tagname using the elementFilter

    Parameters
    ----------
    mesh : UnstructuredMesh
        Mesh to work on
    elementFilter : Union[ElementFilter,FilterOP]
        The element filter to use to select the elements
    tagname : str
        the name of tag to create

    Returns:
    --------
    the mesh with a new tag with, if the tag alreade exist a exception is raised
    """

    elementFilter.mesh = mesh
    for name, data, ids in elementFilter:
        data.tags.CreateTag(tagname).SetIds(ids)

#--------------  CheckIntegrity ---------------
def CheckIntegrity_FilterToETag(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    from BasicTools.Containers.ElementNames import Quadrangle_4
    mesh = CreateSquare(dimensions=[10,10])
    ef = ElementFilter(mesh=mesh,elementType=Quadrangle_4)
    FilterToETag(mesh,ef, "quads")
    if len(mesh.elements[Quadrangle_4].tags["quads"].GetIds()) != 9*9:
        raise # pragma: no cover
    return "ok"

def CheckIntegrity_VerifyExclusiveFilters(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    from BasicTools.Containers.ElementNames import Quadrangle_4
    mesh = CreateSquare(dimensions=[10,10])

    efI = ElementFilter(mesh=mesh,elementType=Quadrangle_4)
    mask = np.ones(mesh.GetNumberOfElements(),dtype=bool)
    efII = ElementFilter(mesh=mesh,mask=mask)
    if VerifyExclusiveFilters([efI,efII],mesh) == 1:
        raise # pragma: no cover

    return "ok"

def CheckIntegrity_ListOfElementFiltersFromMask(GUI=False):

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByElementFilter


    mesh = CreateSquare(dimensions=[10,10])
    mask = np.asarray(np.random.random(mesh.GetNumberOfElements())*5,dtype=int)

    listOfFilters = ListOfElementFiltersFromMask(mask,mesh)
    if VerifyExclusiveFilters(listOfFilters,mesh) == 0:
        raise # pragma: no cover

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
        return "Not OK" # pragma: no cover

    return "OK"

def CheckIntegrity_ListOfElementFiltersFromETagList(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    mesh = CreateSquare()
    tags = ("2D",("X0","X1"),("Y0","Y1"))

    listOfFilters = ListOfElementFiltersFromETagList(tags)

    return "OK"


def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_FilterToETag,
    CheckIntegrity_VerifyExclusiveFilters,
    CheckIntegrity_ListOfElementFiltersFromETagList,
    CheckIntegrity_ListOfElementFiltersFromMask
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok": # pragma: no cover
            return "error in "+str(f) + " res"
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(True)) # pragma: no cover
