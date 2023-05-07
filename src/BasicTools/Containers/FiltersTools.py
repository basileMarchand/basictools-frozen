# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from typing import List, Union, Optional
from numpy.typing import ArrayLike

from BasicTools.Containers.Filters import ElementFilter, FilterOP, FrozenFilter, PartialElementFilter
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh

def GetListOfPartialElementFilter(elementFilter, nbPartitions:int, frozen=True) -> List[PartialElementFilter]:
    """Generate a list of PartialElementFilter for a number of partitions

    Parameters
    ----------
    elementFilter : ElementFilter
        _description_
    nbPartitions : int
        Number of partitions
    frozen : bool, optional
        if True the filters are frozen (evaluated at the moment of creation), no more modification is allowed, by default True

    Returns
    -------
    List[PartialElementFilter]
        a list of partialElementFilter or frozen filter. the union of all the element of the list matches the original elementFilter
    """
    if frozen:
        elementFilterFrozen = elementFilter.GetFrozenFilter()
        res = [ FrozenFilter(elementFilter.mesh) for f in range(nbPartitions) ]
        #print("toto")
        for name, data, ids  in elementFilterFrozen:
            localIds = np.array_split(ids, nbPartitions)
            for r, lids in zip(res,localIds):
                r.SetIdsToTreatFor(data,lids)
    else:
        res = [PartialElementFilter(elementFilter,nbPartitions, i )  for i in range(nbPartitions) ]

    return res

def VerifyExclusiveFilters(listOfFilters: List[Union[ElementFilter,FilterOP]], mesh: UnstructuredMesh) -> bool :
    """Function to check if a list of ElementFilter is exclusive.
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

    for elemName, data in mesh.elements.items():
        mask = np.zeros(data.GetNumberOfElements(),dtype=int)-1
        for fnb,f in enumerate(listOfFilters):
            ids = f.GetIdsToTreat(data)
            if np.any(mask[ids] > -1):
                idd = np.where(mask[ids] > -1)[0][0]
                print(f"Filter {fnb} incompatible with filter {idd} ")
                return 0
            mask[ids] = fnb

    return True

def ListOfElementFiltersFromETagList(tagList: List[Union[str,List[str]]], mesh: Optional[UnstructuredMesh] = None) -> List[ElementFilter]:
    """Function to construct a list of filters from a list of tags (or a list of tags)

    Parameters
    ----------
    tagList : List[Union[str,List[str]]]
        A list containing a string or a list of strings. Example ("tag1",("tag2","tag3")) -> list with 2 ElementFilters
    mesh : Optional[UnstructuredMesh], optional
        Mesh to pass to the filters

    Returns
    -------
    List[ElementFilter]
        List of ElementFilter with the associated tags as filters
    """

    listOfFilters = []
    for matName in tagList:
        if not isinstance(matName, (list, tuple)):
            matName = [matName]
        listOfFilters.append(ElementFilter(mesh=mesh, tags=matName))

    return listOfFilters

def ListOfElementFiltersFromMask(maskVector: ArrayLike, mesh: Optional[UnstructuredMesh]=None) -> List[ElementFilter]:
    """Function to construct a list of filter from a mask vector

    Parameters
    ----------
    maskVector : ArrayLike
        A element size vector of ints to determine to which filter each elements belongs to.
        The number of filters is calculated using np.unique

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
    The tag is added the mesh.

    Parameters
    ----------
    mesh : UnstructuredMesh
        Mesh to work on
    elementFilter : Union[ElementFilter,FilterOP]
        The element filter to use to select the elements
    tagname : str
        the name of tag to create

    if the tag already exists, an exception is raised.
    """

    elementFilter.mesh = mesh
    for name, data, ids in elementFilter:
        data.tags.CreateTag(tagname).SetIds(ids)

def ReadElementFilter(string: str) -> ElementFilter:
    """Function to read from a string all the parameter of a ElementFiler

    Parameters
    ----------
    string : str
        a string in the form of "key(names) & key(options) && ..."
        possible keywords are :
            Tags  -> Tags(bulk1,bulk2)
            nTags -> nTags(points1,points2)
            Dim   -> Dim(3)  (only one int)
            Exprs ->  Exprs(x-1,y+3, x**2 +y**2 - 5)
            nTagsOn : to select the nTags  Treatment
            ExprsOn : to select the Exprs Treatment

    Examples
    --------
    >>> ReadElementFilter("Tags(Inside) & Dim(2) & nTags(Toto,Tata) & Exprs(x-1, y+3, x**2 + y**2 - 5)")

    Please read the documentation of ElementFilter

    Returns
    -------
    ElementFilter
        An initialized instance of an ElementFilter.
    """

    tokens = [ token.strip() for token in string.split("&")]
    res  = ElementFilter()
    if len(string) == 0:
        return res

    #sanity check
    testString = string
    for k in ["nTagsOn","nTags", "Tags", "Dim", "ExprsOn","Exprs"] :
        if testString.count(k)  > 1:
            raise Exception(f"Keyword '{k}' can't be used more than once")
        testString = testString.replace(k,"keyword")

    for t in tokens:
        #sanity check
        nbOpen = t.count("(")
        nbClose = t.count(")")
        if nbOpen != nbClose or nbOpen != 1 or nbClose != 1:
            raise Exception(f"Error parsing token {t}, missing & ")

        keyword = t.split("(")[0]
        ds = t.split("(")[1].split(")")[0]

        data = [ d.strip() for d in ds.split(",")]

        if keyword == "Tags":
            res.SetTags(data)
            continue

        if keyword == "nTags":
            res.SetNTags(data)
            continue

        if keyword == "Dim":
            if len(data) != 1 :
                raise Exception("Cant treat more than one dim")
            res.SetDimensionality(int(data[0]))
            continue

        if keyword == "Exprs":
            from BasicTools.Containers.SymExpr import SymExprWithPos
            res.SetZones([SymExprWithPos(t) for t in data])
            continue

        if keyword == "nTagsOn":
            res.SetNTagsTreatment(data[0])
            continue
        if keyword == "ExprsOn":
            res.SetZoneTreatment(data[0])
            continue

        raise(Exception(f"Cant read the expression '{t}', please check your syntax"))

    return res

#--------------  CheckIntegrity ---------------
def CheckIntegrityReadElementFilter(GUI=False):
    print(ReadElementFilter("Tags(Inside) & Dim(2) & nTags(Toto,Tata) & Exprs(x-1,y+3, x**2 +y**2 - 5)  "))
    print(ReadElementFilter("Tags(Inside) & Dim(3)"))
    print(ReadElementFilter("Dim(3) "))
    print(ReadElementFilter("Exprs(x-3) & ExprsOn(allnodes) & nTagsOn(allnodes) "))
    print(ReadElementFilter(" Tags(toto,tata,titi)"))
    print(ReadElementFilter(""))
    return "ok"

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

    nbElements = 0
    for i,ff in enumerate(listOfFilters):
        new_mesh = ExtractElementsByElementFilter(mesh,ff)
        nbElements += new_mesh.GetNumberOfElements()
        print(f"Number of element in the mesh {i}:" + str(new_mesh.GetNumberOfElements()))
        WriteMeshToXdmf(tmp + f"CI_ListOfElementFiltersFromMask{i}.xdmf", new_mesh)

    WriteMeshToXdmf(tmp + "CI_ListOfElementFiltersFromMask_base.xdmf", mesh,CellFields=[mask],CellFieldsNames=["RANK"])

    print((nbElements,mesh.GetNumberOfElements()))

    if nbElements != mesh.GetNumberOfElements():
        return "Not OK" # pragma: no cover

    return "OK"

def CheckIntegrity_ListOfElementFiltersFromETagList(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateSquare
    mesh = CreateSquare()
    tags = ("2D",("X0","X1"),("Y0","Y1"))

    listOfFilters = ListOfElementFiltersFromETagList(tags)

    return "OK"

def CheckIntegrity(GUI=False):
    toTest= [
        CheckIntegrityReadElementFilter,
        CheckIntegrity_FilterToETag,
        CheckIntegrity_VerifyExclusiveFilters,
        CheckIntegrity_ListOfElementFiltersFromETagList,
        CheckIntegrity_ListOfElementFiltersFromMask
    ]
    for f in toTest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok": # pragma: no cover
            return "error in "+str(f) + " res"
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True)) # pragma: no cover
