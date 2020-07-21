# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh

#missing element in meshio
#EN.Wedge_15
#EN.Pyramid_13

meshIONameByElementName = {}
meshIONameByElementName[EN.Bar_2] = "line"
meshIONameByElementName[EN.Triangle_3] = "triangle"
meshIONameByElementName[EN.Quadrangle_4] = "quad"
meshIONameByElementName[EN.Tetrahedron_4] = "tetra"
meshIONameByElementName[EN.Hexaedron_8] = "hexahedron"
meshIONameByElementName[EN.Wedge_6] = "wedge"
meshIONameByElementName[EN.Pyramid_5] = "pyramid"
meshIONameByElementName[EN.Bar_3] = "line3"
meshIONameByElementName[EN.Triangle_6] = "triangle6"
meshIONameByElementName[EN.Quadrangle_9] = "quad9"
meshIONameByElementName[EN.Tetrahedron_10] = "tetra10"
meshIONameByElementName[EN.Hexaedron_27] = "hexahedron27"
meshIONameByElementName[EN.Wedge_18] = "wedge18"
#meshIONameByElementName[] = "pyramid14"
meshIONameByElementName[EN.Point_1] = "vertex"
meshIONameByElementName[EN.Quadrangle_8] = "quad8"
meshIONameByElementName[EN.Hexaedron_20] = "hexahedron20"
#meshIONameByElementName[] = "triangle10"
#meshIONameByElementName[] = "triangle15"
#meshIONameByElementName[] = "triangle21"
#meshIONameByElementName[] = "line4"
#meshIONameByElementName[] = "line5"
#meshIONameByElementName[] = "line6"
#meshIONameByElementName[] = "tetra20"
#meshIONameByElementName[] = "tetra35"
#meshIONameByElementName[] = "tetra56"
#meshIONameByElementName[] = "quad16"
#meshIONameByElementName[] = "quad25"
#meshIONameByElementName[] = "quad36"
#meshIONameByElementName[] = "triangle28"
#meshIONameByElementName[] = "triangle36"
#meshIONameByElementName[] = "triangle45"
#meshIONameByElementName[] = "triangle55"
#meshIONameByElementName[] = "triangle66"
#meshIONameByElementName[] = "quad49"
#meshIONameByElementName[] = "quad64"
#meshIONameByElementName[] = "quad81"
#meshIONameByElementName[] = "quad100"
#meshIONameByElementName[] = "quad121"
#meshIONameByElementName[] = "line7"
#meshIONameByElementName[] = "line8"
#meshIONameByElementName[] = "line9"
#meshIONameByElementName[] = "line10"
#meshIONameByElementName[] = "line11"
#meshIONameByElementName[] = "tetra84"
#meshIONameByElementName[] = "tetra120"
#meshIONameByElementName[] = "tetra165"
#meshIONameByElementName[] = "tetra220"
#meshIONameByElementName[] = "tetra286"
#meshIONameByElementName[] = "wedge40"
#meshIONameByElementName[] = "wedge75"
#meshIONameByElementName[] = "hexahedron64"
#meshIONameByElementName[] = "hexahedron125"
#meshIONameByElementName[] = "hexahedron216"
#meshIONameByElementName[] = "hexahedron343"
#meshIONameByElementName[] = "hexahedron512"
#meshIONameByElementName[] = "hexahedron729"
#meshIONameByElementName[] = "hexahedron1000"
#meshIONameByElementName[] = "wedge126"
#meshIONameByElementName[] = "wedge196"
#meshIONameByElementName[] = "wedge288"
#meshIONameByElementName[] = "wedge405"
#meshIONameByElementName[] = "wedge550"

elementNameByMeshIOName = {}

for key,meshioname in meshIONameByElementName.items():
    elementNameByMeshIOName[meshioname] = key

def MeshToMeshIO(mesh, TagsAsFields=False):
    import meshio

    points = mesh.nodes

    mesh.ComputeGlobalOffset()

    cells = []
    for name, data in mesh.elements.items():
        cells.append( (meshIONameByElementName[name],data.connectivity) )

    cell_sets = {}
    for name in mesh.GetNamesOfElemTags():
        cell_sets[name] = mesh.GetElementsInTag(name)

    point_sets = {}
    for tag in mesh.nodesTags:
        point_sets[tag.name] = tag.GetIds()

    mesh = meshio.Mesh(points, cells, point_data = mesh.nodeFields.copy(), cell_data = mesh.elemFields.copy(), point_sets=point_sets, cell_sets=cell_sets )

    return mesh

def MeshIOToMesh(mesh, TagsAsFields=False):
    res = UnstructuredMesh()

    res.nodes = mesh.points
    res.nodeFields = mesh.point_data

    for tagname, tagdata in mesh.point_sets.items():
        res.nodesTags.CreateTag(tagname).SetIds(tagdata)

    cpt = 0
    for name, data in mesh.cells:
        elems = res.GetElementsOfType( elementNameByMeshIOName[name] )
        elems.connectivity = data
        nbelems = data.shape[0]
        elems.cpt = nbelems
        for tagname,tagdata in mesh.cell_sets.items():
            #ids = list(filter(lambda x : x >= cpt and x < nbelems +cpt ,tagdata))
            ids = [t for t in tagdata if t >= cpt and t < nbelems +cpt ]
            if len(ids) :
                elems.tags.CreateTag(tagname).SetIds(ids)
        cpt += nbelems

    res.GenerateManufacturedOriginalIDs()
    res.PrepareForOutput()

    return res

    #cell_data = mesh.elemFields

readers = {}
def InitAllReaders():
    import meshio
    from meshio import _helpers as HH
    for ext, filetype in HH.extension_to_filetype.items():
        def Getinit(filetype):
            def __init__(self):
                super(type(self),self).__init__()
                self.canHandleTemporal = False
                self.meshIOInternalfiletype = filetype
            return __init__

        def GetSetFileName():
            def SetFileName(self,fileName):
                self.filename = fileName
            return SetFileName

        def GetRead():
            def Read(self,out=None):
                mesh_io  = meshio.Mesh.read(self.filename, self.meshIOInternalfiletype)
                return MeshIOToMesh(mesh_io)
            return Read

        wrapperClassName = "MeshIO_"+ ext[1:]+"_Reader"
        obj = type(wrapperClassName,
              (object,),
              {"__init__":Getinit(filetype),
               "SetFileName": GetSetFileName(),
               "Read": GetRead(),}
               )

        globals()[wrapperClassName] = obj
        globals()["readers"][ext] = obj

def AddReadersToBasicToolsFactory():
    from BasicTools.IO.IOFactory import RegisterReaderClass
    for ext,obj in globals()["readers"].items():
        RegisterReaderClass(ext,obj, withError= False)

writers = {}
def InitAllWriters():
    from meshio import _helpers as HH

    for ext, filetype in HH.extension_to_filetype.items():

        def GetInit(filetype):
            def __init__(self):
                super(type(self),self).__init__()
                self.canHandleTemporal = False
            return __init__

        def GetSetFileName():
            def SetFileName(self,fileName):
                self.filename = fileName
            return SetFileName

        def GetWrite(meshIOInternalfiletype):
            def Write(self,mesh=None, PointFieldsNames=None, PointFields=None, CellFieldsNames=None, CellFields=None):
                meshio  = MeshToMeshIO(mesh)

                if PointFieldsNames is not None:
                    for n,d in zip(PointFieldsNames, PointFields):
                        meshio.point_data[n] = d

                if CellFieldsNames is not None:
                    for n,d in zip(CellFieldsNames, CellFields):
                        mesh.cell_data[n] = d

                meshio.write(self.filename,meshIOInternalfiletype)
            return Write

        def Open(self):
            pass

        def Close(self):
            pass

        def SetBinary(self):
            pass


        wrapperClassName = "MeshIO_"+ ext[1:]+"_Writer"
        obj = type(wrapperClassName,
              (object,),
              {"__init__":GetInit(filetype),
               "SetFileName": GetSetFileName(),
               "Write": GetWrite(filetype),
               "Open": Open,
               "Close": Close,}
               )

        globals()[wrapperClassName] = obj
        globals()["writers"][ext] = obj

def AddWritersToBasicToolsFactory():
    from BasicTools.IO.IOFactory import RegisterWriterClass
    for ext,obj in globals()["writers"].items():
        RegisterWriterClass(ext,obj, withError= False)

def CheckIntegrity():

    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("MESHIO_NO_FAIL"): return "ok"

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,2,3],[0,1,0]], [[0,1,2]])
    res.nodesTags.CreateTag('point3').SetIds([3])
    elements = res.GetElementsOfType(EN.Bar_2)
    elements.AddNewElement([1,2],1)

    print(res)
    iomesh = MeshToMeshIO(res)
    print("----------")
    res2 = MeshIOToMesh(iomesh)
    print(res2)
    print("--------------------------------")
    iomesh2 = MeshToMeshIO(res2)
    print(iomesh )
    print("----------")
    print(iomesh2 )

    InitAllReaders()
    print(readers)
    #AddReadersToBasicToolsFactory()

    InitAllWriters()
    print(writers)
    #AddWritersToBasicToolsFactory()

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity()) #pragma: no cover
