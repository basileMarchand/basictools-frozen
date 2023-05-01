"""Ansys CDB format writer
"""

import numpy as np
from BasicTools.Containers.Filters import ElementFilter, ElementCounter


def ExportElementTagsInCDBFormat(mesh, tagnames, filename=None):
    """Function to export multiple elements tags to Ansys CDB format

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh from which elements ids for each tag are searched
    tagnames : list[str]
        list of element tags
    filename : str, optional
        filename in CDB format where the element ids for each tag are exported, by default None
    """
    ft = False
    for tname in tagnames:
        ExportElementTagInCDBFormat(mesh, tname, filename=filename,append=ft)
        ft = True

def ExportElementTagInCDBFormat(mesh, tagname, filename=None, append=False):
    """
    Functions to export one element tag as Ansys CDB format. It uses the originals ids

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh from which elements ids for each tag are searched
    tagname : str
        element tag
    filename : str, optional
        filename in CDB format where the element ids for the provided tag are exported, by default None
    append : str, optional
        selects if the info is appended to the file (a+) or overwritten (w or False), by default False
    """
    l = []
    ef = ElementFilter(mesh, tag = tagname)
    for elemtype, data,ids  in ef:
        l.extend(data.originalIds[ids])

    if len(l) == 0:
        print(f"Empty tag {tagname}")
        return

    if filename is None:
        filename = tagname+ ".cdb"

    if append:
        mode = "a+"
    else:
        mode = "w"

    fh = open(filename, mode)
    NB_ELEM_DANS_GROUPE = len(l)
    fh.write("CMBLOCK,{NOM_GROUPE},ELEM,{NB_ELEM_DANS_GROUPE}\n".format(NOM_GROUPE=tagname, NB_ELEM_DANS_GROUPE=NB_ELEM_DANS_GROUPE))
    fh.write("(8i10)\n")

    l = np.array(l)
    stop = 0
    for cpt in range(len(l)//8):
        start = cpt*8
        stop = (cpt+1)*8
        if stop >= len(l):
            stop = len(l)-1
        np.savetxt(fh,l[start:stop][np.newaxis], fmt="%10i",delimiter="")

    np.savetxt(fh,l[stop:][np.newaxis], fmt="%10i",delimiter="")
    fh.close()


def CheckIntegrity():
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube()
    mesh.GenerateManufacturedOriginalIDs()

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    ExportElementTagsInCDBFormat(mesh, mesh.elements.GetTagsNames(), filename=tempdir+"/CheckIntegrity_AnsysWriter.cdb")

    return "ok"
