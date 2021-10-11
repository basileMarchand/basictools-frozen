
import numpy as np
from BasicTools.Containers.Filters import ElementFilter, ElementCounter


def ExportElementTagsInCDBFormat(mesh, tagnames, filename=None):

    ft = False
    for tname in tagnames:
        ExportElementTagInCDBFormat(mesh, tname, filename=filename,append=ft)
        ft = True

def ExportElementTagInCDBFormat(mesh, tagname, filename=None,append=False):
    """
    Functions to export elements Tags as Ansys CDB format.
    It uses the originals ids. 
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




