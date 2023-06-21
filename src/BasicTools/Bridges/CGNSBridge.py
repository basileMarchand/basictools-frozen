# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh, ElementsContainer
import BasicTools.Containers.ElementNames as EN

CGNSNameToBasicTools = {}
CGNSNameToBasicTools["NODE"] = EN.Point_1
# 1D
CGNSNameToBasicTools["BAR_2"] = EN.Bar_2
CGNSNameToBasicTools["BAR_3"] = EN.Bar_3
# 2D
CGNSNameToBasicTools["TRI_3"] = EN.Triangle_3
CGNSNameToBasicTools["TRI_6"] = EN.Triangle_6
CGNSNameToBasicTools["QUAD_4"] = EN.Quadrangle_4
CGNSNameToBasicTools["QUAD_8"] = EN.Quadrangle_8
CGNSNameToBasicTools["QUAD_9"] = EN.Quadrangle_9
# 3D
CGNSNameToBasicTools["TETRA_4"] = EN.Tetrahedron_4
CGNSNameToBasicTools["TETRA_10"] = EN.Tetrahedron_10
CGNSNameToBasicTools["PYRA_5"] = EN.Pyramid_5
CGNSNameToBasicTools["PYRA_13"] = EN.Pyramid_13
CGNSNameToBasicTools["PENTA_6"] = EN.Wedge_6
CGNSNameToBasicTools["PENTA_15"] = EN.Wedge_15
CGNSNameToBasicTools["PENTA_18"] = EN.Wedge_18
CGNSNameToBasicTools["HEXA_8"] = EN.Hexaedron_8
CGNSNameToBasicTools["HEXA_20"] = EN.Hexaedron_20
CGNSNameToBasicTools["HEXA_27"] = EN.Hexaedron_27
# when adding a new element please add the Pertmutaion
#https://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
#example  CGNSToBasicToolsPermutation["TETRA_4"] = [0 3 2 1]
CGNSToBasicToolsPermutation = {}
CGNSToBasicToolsPermutation["HEXA_20"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17 ,18, 19, 12, 13, 14, 15]
CGNSToBasicToolsPermutation["HEXA_27"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17 ,18, 19, 12, 13, 14, 15, 24, 22, 21, 23, 20, 25, 26]

CGNS_ElementType_l = ['Null', 'UserDefined', 'NODE', 'BAR_2', 'BAR_3', 'TRI_3', 'TRI_6', 'QUAD_4', 'QUAD_8', 'QUAD_9', 'TETRA_4', 'TETRA_10', 'PYRA_5', 'PYRA_14', 'PENTA_6', 'PENTA_15', 'PENTA_18', 'HEXA_8', 'HEXA_20', 'HEXA_27', 'MIXED', 'PYRA_13', 'NGON_n', 'NFACE_n', 'BAR_4', 'TRI_9', 'TRI_10', 'QUAD_12', 'QUAD_16', 'TETRA_16', 'TETRA_20', 'PYRA_21', 'PYRA_29', 'PYRA_30', 'PENTA_24', 'PENTA_38', 'PENTA_40', 'HEXA_32', 'HEXA_56', 'HEXA_64', 'BAR_5', 'TRI_12', 'TRI_15', 'QUAD_P4_16', 'QUAD_25', 'TETRA_22', 'TETRA_34', 'TETRA_35', 'PYRA_P4_29', 'PYRA_50', 'PYRA_55', 'PENTA_33', 'PENTA_66', 'PENTA_75', 'HEXA_44', 'HEXA_98', 'HEXA_125']

CGNSNumberToBasicTools = { CGNS_ElementType_l.index(k):v for k,v in CGNSNameToBasicTools.items() }

BasicToolsToCGNSNames = { y:x for x,y in CGNSNameToBasicTools.items() }
BasicToolsToCGNSNumber = { y:x for x,y in CGNSNumberToBasicTools.items() }


def GetCGNSNumberToBasicTools(cgnsNumber):
    res = CGNSNumberToBasicTools.get(cgnsNumber,None)
    if res is None:
        raise Exception(f"CGNC elements of number : {cgnsNumber} not coded yet")
    return res


def __ReadIndex(pyTree):
    a = __ReadIndexArray(pyTree)
    b =  __ReadIndexRange(pyTree)
    return np.hstack( (a,b) )

def __ReadIndexArray(pyTree):
    indexArrayPaths = GetAllNodesByTypeSet(pyTree, 'IndexArray_t')
    res =[]
    for indexArrayPath in indexArrayPaths:
        data = GetNodeByPath(pyTree, indexArrayPath)[0]
        res.extend(data[1].ravel())
    return np.array(res, dtype=PBasicIndexType).ravel()

def __ReadIndexRange(pyTree):
    indexRangePaths = GetAllNodesByTypeSet(pyTree, 'IndexRange_t')
    res =[]
    #if len(indexRangePaths) == 0:
    #    return np.zeros((0), dtype = PBasicIndexType)

    for indexRangePath in indexRangePaths:
        indexRange = GetNodeByPath(pyTree, indexRangePath)[0][1]
        begin = indexRange[0]#[:,0]
        end = indexRange[1]#[:,1]
        res.extend(np.arange(begin, end+1).ravel())

    return np.array(res, dtype = PBasicIndexType).ravel()

def GetRecursiveObjects(pyTree, filter, state= None):
    res = list()
    state, lres = filter(state, pyTree)
    res.extend(lres)
    for subTree in pyTree[2]:
        lstate, lres = GetRecursiveObjects(subTree,filter,state)
        res.extend(lres)
    return state, res


def GetAllNodesByTypeSet(pyTree, typeset):

    def fil(state, data):
        if state is None:
            state = ""
        state += "/" + data[0]
        if data[3] == typeset:
            return state, [state]
        else:
            return state, []

    return GetRecursiveObjects(pyTree, fil)[1]

def GetNodeByPath(pyTree, path):

    def fil(state, data):
        if len(state)==1 and data[0] == state[0]:
            return state[1:], [data]
        else:
            return state[1:], []

    return GetRecursiveObjects(pyTree, fil, path.strip("/").split("/"))[1]

def CGNSToMesh(pyTree, baseNumberOrName= 0, zoneNumberOrName = 0)-> UnstructuredMesh:

    res = UnstructuredMesh()

    if type(baseNumberOrName) is int:
        basepath = GetAllNodesByTypeSet(pyTree,"CGNSBase_t")[baseNumberOrName]
    else:
        basepath = baseNumberOrName

    base = GetNodeByPath(pyTree,basepath)[0]

    if type(zoneNumberOrName) is int:
        zonepath = GetAllNodesByTypeSet(base,"Zone_t")[zoneNumberOrName]
    else:
        zonepath = zoneNumberOrName

    zonePyTree = GetNodeByPath(base,"/"+ zonepath)[0]


    gridCoordinatesPath = GetAllNodesByTypeSet(zonePyTree, "GridCoordinates_t")[0]

    gx = GetNodeByPath(zonePyTree,gridCoordinatesPath+'/CoordinateX')[0][1]
    res.nodes = np.empty((gx.shape[0],3), dtype= PBasicFloatType)
    res.nodes[:,0] = gx
    res.nodes[:,1] = GetNodeByPath(zonePyTree,gridCoordinatesPath+'/CoordinateY')[0][1]
    res.nodes[:,2] = GetNodeByPath(zonePyTree,gridCoordinatesPath+'/CoordinateZ')[0][1]
    res.originalIDNodes = np.arange(1, res.nodes.shape[0]+1, dtype= PBasicIndexType)

    #elements
    elementsPaths = GetAllNodesByTypeSet(zonePyTree, "Elements_t")

    for elementsPath in elementsPaths:
        cgnsElements = GetNodeByPath(zonePyTree,elementsPath)[0]
        cgnsUserName = "Elements_" + cgnsElements[0]
        cgnsElemType = cgnsElements[1][0]
        originalIds = __ReadIndexRange(cgnsElements)

        if cgnsElemType == 20: # we are in mixed mode:
            cpt = 0
            elementConnectivity = GetNodeByPath(cgnsElements,cgnsUserName+'/ElementConnectivity')[0][1]
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
            elementConnectivity = np.asarray(GetNodeByPath(cgnsElements,cgnsUserName+'/ElementConnectivity')[0][1], dtype=PBasicIndexType).reshape( (-1, EN.numberOfNodes[basicToolsElemType]  ) )-1
            elems: ElementsContainer = res.elements.GetElementsOfType(basicToolsElemType)

            if CGNSToBasicToolsPermutation.get(cgnsElemType,None) is not None:
                elementConnectivity = elementConnectivity[:,CGNSToBasicToolsPermutation[cgnsElemType]]
            elcpt = elems.cpt

            nelcpt = elems.AddNewElements(elementConnectivity, originalIds)

            elems.tags.CreateTag(cgnsUserName,False).AddToTag(np.arange(elcpt,nelcpt))

    datasPaths = GetAllNodesByTypeSet(zonePyTree, "FlowSolution_t")
    for dataPath in datasPaths :
        datas = GetNodeByPath(zonePyTree,dataPath)[0]
        gl = GetAllNodesByTypeSet(datas, "GridLocation_t")

        store  = res.nodeFields
        if len(gl) > 0:
            if "".join(np.array(GetNodeByPath(datas,gl[0])[0][1], dtype=str) )== "CellCenter":
                store  = res.elemFields
            if "".join(np.array(GetNodeByPath(datas,gl[0])[0][1], dtype=str) )== "Vertex":
                store  = res.nodeFields

        fieldPaths = GetAllNodesByTypeSet(datas, "DataArray_t")
        for fieldPath in fieldPaths:
            fieldData = GetNodeByPath(datas,fieldPath)[0]
            dataName = fieldData[0]
            data = fieldData[1]
            if  dataName == "OriginalIds" and store is res.nodeFields:
                res.originalIDNodes =   np.asarray(data, dtype=PBasicIndexType,order='C')
            elif dataName == "OriginalIds" and store is res.elemFields:
               res.SetElementsOriginalIDs(np.asarray(data, dtype=PBasicIndexType))
            else:
                res.nodeFields[dataName] = np.asarray(data)

    ZoneBCPaths = GetAllNodesByTypeSet(zonePyTree, "ZoneBC_t")
    for ZoneBCPath in ZoneBCPaths:
        ZoneBC = GetNodeByPath(zonePyTree,ZoneBCPath)[0]

        BCPaths = GetAllNodesByTypeSet(ZoneBC, "BC_t")
        for BCPath in BCPaths:
            BCNode = GetNodeByPath(ZoneBC,BCPath)[0]
            BCName = str(BCNode[0])

            indices = __ReadIndex(BCNode)
            if len(indices) == 0:
                continue

            gl = GetAllNodesByTypeSet(BCNode, "GridLocation_t")
            if "".join(np.array(GetNodeByPath(BCNode,gl[0])[0][1], dtype=str) )== "CellCenter":
                #res.AddElementToTagUsingOriginalId(indices-1, BCName)
                #res.AddElementsToTag(indices-1, BCName)
                print(f"Reading Element selections is not supported. skiping {BCName}")
                pass
            if "".join(np.array(GetNodeByPath(BCNode,gl[0])[0][1], dtype=str) )== "Vertex":
                res.nodesTags.CreateTag(BCName).SetIds(indices-1)

    res.PrepareForOutput()
    return res

def NewDataArray(father, name, data):
    res =  [name, data, [], 'DataArray_t']
    father[2].append(res)
    return res

def NewBC(father,name,data, pttype, onPoints=True):
    res =  [name, np.array([b'F', b'a', b'm', b'i', b'l', b'y', b'S', b'p', b'e', b'c', b'i', b'f', b'i', b'e', b'd'], dtype='|S1'), [], 'BC_t']

    plist = ['PointList', data[None,:], [], pttype]
    res[2].append(plist)

    family = ['FamilyName', np.array([b'N', b'u', b'l', b'l'], dtype='|S1'), [], 'FamilyName_t']
    res[2].append(family)

    if onPoints:
        glocation = ['GridLocation', np.array([b'V', b'e', b'r', b't', b'e', b'x'], dtype='|S1'), [], 'GridLocation_t']
    else:
        glocation = ['GridLocation', np.array([b'C', b'e', b'l', b'l', b'C', b'e', b'n', b't', b'e', b'r'], dtype='|S1'), [], 'GridLocation_t']
    res[2].append(glocation)

    father[2].append(res)
    return res

def NewElements(father, data):
    name = data.elementType
    nbelem = data.GetNumberOfElements()
    res = ["Elements_"+ BasicToolsToCGNSNames[name], np.array([BasicToolsToCGNSNumber[name],  0], dtype=np.int32), [['ElementRange', np.array((data.globaloffset+1, data.globaloffset+nbelem) ), [], 'IndexRange_t'], ['ElementConnectivity', (data.connectivity+1).ravel(), [], 'DataArray_t']], 'Elements_t']
    father[2].append(res)
    return res

def NewFlowSolution(father,name,gridlocation ):
    res = [name, None,[['GridLocation', np.array([ c for c in gridlocation], dtype='|S1'), [], 'GridLocation_t']] , 'FlowSolution_t']
    father[2].append(res)
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
        outputPyTree=['CGNSTree', None, [], 'CGNSTree_t']
        version = ['CGNSLibraryVersion', np.array([3.4], dtype=np.float32), [], 'CGNSLibraryVersion_t']
        outputPyTree[2].append(version)
        physicalDim = mesh.GetPointsDimensionality()
        topologicalDim = mesh.GetElementsDimensionality()
        nbElemTags = len(mesh.GetNamesOfElemTags())
        nbNodesTags = len(mesh.nodesTags)
        totalNumberOfTags = nbElemTags +nbNodesTags

        base = [baseName, np.array([physicalDim, topologicalDim], dtype=np.int32), [], 'CGNSBase_t']  # name, physical dim, topological dim
        outputPyTree[2].append(base)

        family = ['Bulk', np.array([b'B', b'u', b'l', b'k'], dtype='|S1'), [], 'FamilyName_t']
        base[2].append(family)

        s = np.array([[mesh.GetNumberOfNodes(),mesh.GetNumberOfElements(),totalNumberOfTags]],dtype=np.int32)
        grid = [zoneName, s, [['ZoneType', np.array([b'U', b'n', b's', b't', b'r', b'u', b'c', b't', b'u', b'r', b'e',b'd'], dtype='|S1'), [], 'ZoneType_t']], 'Zone_t']
        base[2].append(grid)
        zone_family = ['FamilyName', np.array([b'B', b'u', b'l', b'k'], dtype='|S1'), [], 'FamilyName_t']
        grid[2].append(zone_family)

    else:
        grid = GetNodeByPath(outputPyTree,f"/{baseName}/{zoneName}/")[0][1]

    #add nodes
    gridCoordinates = ['GridCoordinates', None, [], 'GridCoordinates_t']
    grid[2].append(gridCoordinates)

    for i,coord in enumerate(["X", "Y", "Z"]):
        da = NewDataArray(gridCoordinates, "Coordinate"+coord, np.copy(mesh.nodes[:,i]) )

    #nodes tags
    if len(mesh.nodesTags):
        zbg = ['Points_Selections', None, [], 'ZoneBC_t']
        grid[2].append(zbg)
        for tag in mesh.nodesTags:
           NewBC(zbg,tag.name,tag.GetIds()+1,pttype="IndexArray_t")

    for name, data in mesh.elements.items():
        #elem tags
        nbelem = data.GetNumberOfElements()
        if nbelem == 0:
            continue
        NewElements(grid, data)

    ezbg = ['Elements_Selections', None, [], 'ZoneBC_t']
    grid[2].append(ezbg)

    for name in mesh.GetNamesOfElemTags():
        ids = mesh.GetElementsInTag(name)+1
        NewBC(ezbg,name,ids, pttype="IndexArray_t", onPoints=False)

    flowSolution = NewFlowSolution(grid, "PointData", gridlocation="Vertex")
    NewDataArray(flowSolution, "OriginalIds", mesh.originalIDNodes )
    if len(mesh.nodeFields) > 0:
        for name, pointData in mesh.nodeFields.items():
            if pointData.dtype.char == "U":
                print(f"skipping nodeFields '{name}' because if of type: {pointData.dtype}"  )
                continue
            NewDataArray(flowSolution, name,pointData)

    ## elem fields
    flowSolutionCell = NewFlowSolution(grid, "CellData", gridlocation="CellCenter")
    NewDataArray(flowSolutionCell, "OriginalIds", mesh.GetElementsOriginalIDs() )
    if len(mesh.elemFields) > 0:
        for name, cellData in mesh.elemFields.items():
            if cellData.dtype.char == "U":
                print(f"skipping elemFields '{name}' because if of type: {cellData.dtype}"  )
                continue
            NewDataArray(flowSolutionCell, name, np.copy(cellData))


    return outputPyTree

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import GetUniqueTempFile
    from BasicTools.Actions.OpenInParaView import OpenInParaView
    import BasicTools.TestData as BasicToolsTestData


    import BasicTools.IO.GeofReader as GR
    res = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    res = CreateCube()


    res.VerifyIntegrity()
    resII = CGNSToMesh(MeshToCGNS(res))
    resII.VerifyIntegrity()

    from BasicTools.Containers.MeshTools import IsClose
    #IsClose(res,resII)
    #for the moment the element tags are not exported

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
