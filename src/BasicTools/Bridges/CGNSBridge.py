# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from CGNS.MAP import load, save
import CGNS.PAT.cgnskeywords as CGK
import CGNS.PAT.cgnsutils as CGU
import CGNS.PAT.cgnslib as CGL

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh, ElementsContainer
import BasicTools.Containers.ElementNames as EN

CGNSNameToBasicTools = {}
CGNSNameToBasicTools[CGK.NODE_s] = EN.Point_1
# 1D
CGNSNameToBasicTools[CGK.BAR_2_s] = EN.Bar_2
CGNSNameToBasicTools[CGK.BAR_3_s] = EN.Bar_3
# 2D
CGNSNameToBasicTools[CGK.TRI_3_s] = EN.Triangle_3
CGNSNameToBasicTools[CGK.TRI_6_s] = EN.Triangle_6
CGNSNameToBasicTools[CGK.QUAD_4_s] = EN.Quadrangle_4
CGNSNameToBasicTools[CGK.QUAD_8_s] = EN.Quadrangle_8
CGNSNameToBasicTools[CGK.QUAD_9_s] = EN.Quadrangle_9
# 3D
CGNSNameToBasicTools[CGK.TETRA_4_s] = EN.Tetrahedron_4
CGNSNameToBasicTools[CGK.TETRA_10_s] = EN.Tetrahedron_10
CGNSNameToBasicTools[CGK.PYRA_5_s] = EN.Pyramid_5
CGNSNameToBasicTools[CGK.PYRA_13_s] = EN.Pyramid_13
CGNSNameToBasicTools[CGK.PENTA_6_s] = EN.Wedge_6
CGNSNameToBasicTools[CGK.PENTA_15_s] = EN.Wedge_15
CGNSNameToBasicTools[CGK.PENTA_18_s] = EN.Wedge_18
CGNSNameToBasicTools[CGK.HEXA_8_s] = EN.Hexaedron_8
CGNSNameToBasicTools[CGK.HEXA_20_s] = EN.Hexaedron_20
CGNSNameToBasicTools[CGK.HEXA_27_s] = EN.Hexaedron_27
# when adding a new element please add the Pertmutaion
#https://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
#example  CGNSToBasicToolsPermutation[CGK.TETRA_4_s] = [0 3 2 1]
CGNSToBasicToolsPermutation = {}
CGNSToBasicToolsPermutation[CGK.HEXA_20_s] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17 ,18, 19, 12, 13, 14, 15]
CGNSToBasicToolsPermutation[CGK.HEXA_27_s] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17 ,18, 19, 12, 13, 14, 15, 24, 22, 21, 23, 20, 25, 26]


CGNSNumberToBasicTools = { CGK.ElementType_l.index(k):v for k,v in CGNSNameToBasicTools.items() }

BasicToolsToCGNSNames = { y:x for x,y in CGNSNameToBasicTools.items() }
BasicToolsToCGNSNumber = { y:x for x,y in CGNSNumberToBasicTools.items() }


def GetCGNSNumberToBasicTools(cgnsNumber):
    res = CGNSNumberToBasicTools.get(cgnsNumber,None)
    if res is None:
        raise Exception(f"Elements of type : {CGK.ElementType_l[cgnsNumber]} ({cgnsNumber}) not coded yet")
    return res


def __ReadIndex(pyTree):
    a = __ReadIndexArray(pyTree)
    b =  __ReadIndexRange(pyTree)
    return np.hstack( (a,b) )

def __ReadIndexArray(pyTree):
    indexArrayPaths = CGU.getAllNodesByTypeSet(pyTree, ['IndexArray_t'])
    res =[]
    for indexArrayPath in indexArrayPaths:
        data = CGU.getNodeByPath(pyTree, indexArrayPath)
        res.extend(data[1].ravel())
    return np.array(res, dtype=PBasicIndexType).ravel()

def __ReadIndexRange(pyTree):
    indexRangePaths = CGU.getAllNodesByTypeSet(pyTree, ['IndexRange_t'])
    res =[]
    #if len(indexRangePaths) == 0:
    #    return np.zeros((0), dtype = PBasicIndexType)

    for indexRangePath in indexRangePaths:
        indexRange = CGU.getNodeByPath(pyTree, indexRangePath)[1]
        begin = indexRange[0]#[:,0]
        end = indexRange[1]#[:,1]
        res.extend(np.arange(begin, end+1).ravel())

    return np.array(res, dtype = PBasicIndexType).ravel()

def CGNSToMesh(pyTree, baseNumberOrName= 0, zoneNumberOrName = 0)-> UnstructuredMesh:

    res = UnstructuredMesh()

    if type(baseNumberOrName) is int:
        basepath = CGU.getAllNodesByTypeSet(pyTree,[CGK.CGNSBase_ts])[baseNumberOrName][1:]
    else:
        basepath = baseNumberOrName

    base = CGU.getNodeByPath(pyTree,basepath)

    if type(zoneNumberOrName) is int:
        zonepath = CGU.getAllNodesByTypeSet(base,[CGK.Zone_t])[zoneNumberOrName]
    else:
        zonepath = zoneNumberOrName

    zonePyTree = CGU.getNodeByPath(pyTree,"/"+ zonepath)

    gridCoordinatesPath = CGU.getAllNodesByTypeSet(zonePyTree, [CGK.GridCoordinates_ts])[0]

    gx = CGU.getNodeByPath(zonePyTree,gridCoordinatesPath+'/CoordinateX')[1]
    res.nodes = np.empty((gx.shape[0],3), dtype= PBasicFloatType)
    res.nodes[:,0] = gx
    res.nodes[:,1] = CGU.getNodeByPath(zonePyTree,gridCoordinatesPath+'/CoordinateY')[1]
    res.nodes[:,2] = CGU.getNodeByPath(zonePyTree,gridCoordinatesPath+'/CoordinateZ')[1]
    res.originalIDNodes = np.arange(1, res.nodes.shape[0]+1, dtype= PBasicIndexType)

    #elements
    elementsPaths = CGU.getAllNodesByTypeSet(zonePyTree, [CGK.Elements_ts])

    for elementsPath in elementsPaths:
        cgnsElements = CGU.getNodeByPath(zonePyTree,elementsPath)
        cgnsUserName = "Elements_" + cgnsElements[0]
        cgnsElemType = cgnsElements[1][0]
        originalIds = __ReadIndexRange(cgnsElements)

        if cgnsElemType == 20: # we are in mixed mode:
            cpt = 0
            elementConnectivity = CGU.getNodeByPath(cgnsElements,'/ElementConnectivity')[1]
            for oid in originalIds:
                basicToolsElemType = GetCGNSNumberToBasicTools(elementConnectivity[cpt])
                cpt +=1
                nbp = EN.numberOfNodes[basicToolsElemType]
                coon = elementConnectivity[cpt:cpt+nbp] - 1
                if CGNSToBasicToolsPermutation.get(cgnsElemType,None) is not None:
                    coon = coon[CGNSToBasicToolsPermutation[cgnsElemType]]
                cpt += nbp
                elems: ElementsContainer = res.elements.GetElementsOfType(basicToolsElemType)
                elcpt = elems.AddNewElements(coon[None,:], [oid] )
                elems.tags.CreateTag(cgnsUserName,False).AddToTag(elcpt-1)
        else:
            basicToolsElemType = GetCGNSNumberToBasicTools(cgnsElemType)
            elementConnectivity = np.asarray(CGU.getNodeByPath(cgnsElements,'/ElementConnectivity')[1], dtype=PBasicIndexType).reshape( (-1, EN.numberOfNodes[basicToolsElemType]  ) )-1
            elems: ElementsContainer = res.elements.GetElementsOfType(basicToolsElemType)

            if CGNSToBasicToolsPermutation.get(cgnsElemType,None) is not None:
                elementConnectivity = elementConnectivity[:,CGNSToBasicToolsPermutation[cgnsElemType]]
            elcpt = elems.cpt

            nelcpt = elems.AddNewElements(elementConnectivity, originalIds)

            elems.tags.CreateTag(cgnsUserName,False).AddToTag(np.arange(elcpt,nelcpt))

    datasPaths = CGU.getAllNodesByTypeSet(zonePyTree, [CGK.FlowSolution_ts])
    for dataPath in datasPaths :
        datas = CGU.getNodeByPath(zonePyTree,dataPath)
        gl = CGU.getAllNodesByTypeSet(datas, [CGK.GridLocation_ts])

        store  = res.nodeFields
        if len(gl) > 0:
            if "".join(np.array(CGU.getNodeByPath(datas,gl[0])[1], dtype=str) )== CGK.CellCenter_s:
                store  = res.elemFields
            if "".join(np.array(CGU.getNodeByPath(datas,gl[0])[1], dtype=str) )== CGK.Vertex_s:
                store  = res.nodeFields

        fieldPaths = CGU.getAllNodesByTypeSet(datas, [CGK.DataArray_ts])
        for fieldPath in fieldPaths:
            fieldData = CGU.getNodeByPath(datas,fieldPath)
            dataName = fieldData[0]
            data = fieldData[1]
            if  dataName == "OriginalIds" and store is res.nodeFields:
                res.originalIDNodes =   np.asarray(data, dtype=PBasicIndexType,order='C')
            elif dataName == "OriginalIds" and store is res.elemFields:
               res.SetElementsOriginalIDs(np.asarray(data, dtype=PBasicIndexType))
            else:
                res.nodeFields[dataName] = np.asarray(data)

    ZoneBCPaths = CGU.getAllNodesByTypeSet(zonePyTree, [CGK.ZoneBC_ts])
    for ZoneBCPath in ZoneBCPaths:
        ZoneBC = CGU.getNodeByPath(zonePyTree,ZoneBCPath)

        BCPaths = CGU.getAllNodesByTypeSet(ZoneBC, [CGK.BC_ts])
        for BCPath in BCPaths:
            BCNode = CGU.getNodeByPath(ZoneBC,BCPath)
            BCName = str(BCNode[0])


            indices = __ReadIndex(BCNode)
            if len(indices) == 0:
                continue
            if "Point" in BCNode[2][0][0]:
                res.nodesTags.CreateTag(BCName).SetIds(indices-1)
            else:
                res.AddElementToTagUsingOriginalId(indices, BCName)

    res.PrepareForOutput()
    return res

def MeshToCGNS( mesh: UnstructuredMesh, outputPyTree = None,  baseNumberOrName= 0, zoneNumberOrName = 0) :

    if type(baseNumberOrName) is int:
        baseName = f"Base_{baseNumberOrName}"
    else:
        baseName = f"{baseNumberOrName}"

    if type(zoneNumberOrName) is int:
        zoneName = f"Zone_{baseNumberOrName}"
    else:
        zoneName = f"{baseNumberOrName}"

    if outputPyTree is None:
        outputPyTree=CGL.newCGNSTree()
        physicalDim = mesh.GetPointsDimensionality()
        topologicalDim = mesh.GetElementsDimensionality()
        base = CGL.newBase(outputPyTree,baseName,physicalDim,topologicalDim) # name, physical dim, topological dim

        s = np.array([mesh.GetNumberOfNodes(),mesh.GetNumberOfElements(),0],dtype=np.int32)
        grid = CGL.newZone(base,zoneName,s,CGK.Unstructured_s)
    else:
        grid = CGU.getNodeByPath(outputPyTree,f"/{baseName}/{zoneName}/")[1]

    #add nodes
    gridCoordinates = CGL.newGridCoordinates(grid,"GridCoordinates")
    for i,coord in enumerate(["X", "Y", "Z"]):
        da = CGL.newDataArray(gridCoordinates, "Coordinate"+coord, np.copy(mesh.nodes[:,i]) )


    #nodes tags
    if len(mesh.nodesTags):
        zbg = CGL.newZoneBC(grid)
        for tag in mesh.nodesTags:
           tt = CGL.newBC(zbg,tag.name,pttype=CGK.IndexArray_ts)
           tt[2][0][0] = "PointList"
           tt[2][0][1] = tag.GetIds()+1

    for name, data in mesh.elements.items():
        #elem tags
        nbelem = data.GetNumberOfElements()
        if nbelem == 0:
            continue
        CGL.newElements(grid,"Elements_"+ BasicToolsToCGNSNames[name],BasicToolsToCGNSNumber[name], np.array((data.globaloffset+1, data.globaloffset+nbelem) ), econnectivity= (data.connectivity+1).ravel() )


    flowSolution = CGL.newFlowSolution(grid, "PointData", gridlocation=CGK.Vertex_s)
    da = CGL.newDataArray(flowSolution, "OriginalIds", mesh.originalIDNodes )
    for name, pointData in mesh.nodeFields.items():
        da = CGL.newDataArray(flowSolution, name,pointData)

    ## elem fields
    flowSolutionCell = CGL.newFlowSolution(grid, "CellData", gridlocation=CGK.CellCenter_s)
    da = CGL.newDataArray(flowSolutionCell, "OriginalIds", mesh.GetElementsOriginalIDs() )
    for name, cellData in mesh.elemFields.items():
        CGL.newDataArray(flowSolutionCell, name, np.copy(cellData))

    return outputPyTree

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import GetUniqueTempFile
    from BasicTools.Actions.OpenInParaView import OpenInParaView
    import BasicTools.TestData as BasicToolsTestData


    import BasicTools.IO.GeofReader as GR
    res = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")

    res.VerifyIntegrity()
    resII = CGNSToMesh(MeshToCGNS(res))
    resII.VerifyIntegrity()

    from BasicTools.Containers.MeshTools import IsClose
    #IsClose(res,resII)
    #for the moment the element tags are not exported

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
