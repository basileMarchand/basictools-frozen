# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.DofNumbering import ComputeDofNumbering
import BasicTools.Containers.ElementNames as EN
from BasicTools.NumpyDefs import PBasicIndexType

from scipy.sparse import coo_matrix, csr_matrix
from BasicTools.FE.IntegrationsRules import LagrangeIsoParam
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.Containers.Filters import ElementFilter


def GetElementaryMatrixForFormulation(elemName, wform, unknownNames, space=LagrangeSpaceP1,geoFactor=None):
    # Explicitly specify signature to cleanly display default argument values
    # in sphinx autodoc generated documentation
    """
    GetElementaryMatrixForFormulation(elemName, wform, unknownNames, \
            space=BasicTools.FE.Spaces.FESpaces.LagrangeSpaceP1)
    """
    from BasicTools.FE.Integration import IntegrateGeneral

    from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
    mesh = UnstructuredMesh()

    mesh.nodes = np.asarray(space[elemName].posN,dtype=float)
    if geoFactor is not None:
        mesh.nodes = mesh.nodes*geoFactor
    mesh.originalIDNodes = np.arange(0,mesh.GetNumberOfNodes(),dtype=PBasicIndexType)

    elements = mesh.GetElementsOfType(elemName)
    elements.connectivity = np.arange(space[elemName].GetNumberOfShapeFunctions(),dtype=PBasicIndexType)

    elements.connectivity.shape = (1,space[elemName].GetNumberOfShapeFunctions())
    elements.GetTag("3D").AddToTag(0)

    elements.originalIds = np.arange(0,1,dtype=PBasicIndexType)
    elements.cpt = elements.connectivity.shape[0]

    mesh.PrepareForOutput()
    numbering = ComputeDofNumbering(mesh,space,)

    unkownFields = []
    for name in unknownNames:
        unkownFields.append(FEField(name,mesh,space,numbering))

    M,f = IntegrateGeneral(mesh=mesh,wform=wform, unkownFields= unkownFields,constants={},fields=[])

    return M.toarray(), f

def PrepareFEComputation(mesh, elementFilter = None, numberOfComponents = None):
    """
    Prepares a finite element computation with Lagrange isoparametric finite
    elements by computing FESpace, numberings, offset and NGauss

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    elementFilter : ElementFilter, optional
        Filter on which the function is applied
    numberOfComponents : int
        the number of components

    Returns
    -------
    FESpace
    numberings
    offset : list of int
    NGauss : int
    """
    dim = mesh.GetDimensionality()
    if elementFilter == None:
        elementFilter = ElementFilter(mesh)
        elementFilter.SetDimensionality(dim)

    if numberOfComponents == None:
        numberOfComponents = dim

    NGauss = 0
    spaces = LagrangeSpaceGeo

    for name, data, ids in elementFilter:
        p,w =  LagrangeIsoParam[name]
        NGauss += len(ids)*len(w)

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)
    numberings = [numbering]*numberOfComponents

    offset = []
    totaldofs = 0
    for n in numberings:
        offset.append(totaldofs)
        totaldofs += n["size"]

    return spaces, numberings, offset, NGauss

def ComputeL2ScalarProducMatrix(mesh, numberOfComponents, elementFilter = None):
    """
    Computes the L2 scalar product used to compute the correlations
    between the primal solution snapshots. The numberOfComponents
    depends on the solution type: 3 for solid mechanics in 3D, or 1 for
    thermal in any dimension. (Lagrange isoparametric finite elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    numberOfComponents : int
        the number of components
    elementFilter : ElementFilter, optional
        Filter on which the function is applied

    Returns
    -------
    scipy.sparse.csr_matrix
        the sparse matrix of the L2 scalar product
    """

    from BasicTools.FE.Integration import IntegrateGeneral

    dim = mesh.GetDimensionality()

    if elementFilter == None:
        elementFilter = ElementFilter(mesh)
        elementFilter.SetDimensionality(dim)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, elementFilter, numberOfComponents)

    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField


    T = GetField("T",numberOfComponents)
    Tt = GetTestField("T",numberOfComponents)

    wf = T.T*Tt

    constants = {}
    fields  = []

    if numberOfComponents == 1:
        unkownFields = [FEField("T",mesh=mesh,space=spaces,numbering=numberings[0])]
    else:
        unkownFields = [FEField("T_"+str(i),mesh=mesh,space=spaces,numbering=numberings[i]) for i in range(numberOfComponents)]

    K,_ = IntegrateGeneral(mesh=mesh,wform=wf, constants=constants, fields=fields,
                    unkownFields = unkownFields, elementFilter = elementFilter)

    return K

def ComputeH10ScalarProductMatrix(mesh, numberOfComponents):
    """
    Computes the H10 scalar product used to compute the correlations
    between the primal solution snapshots. The numberOfComponents
    depends on the solution type: 3 for solid mechanics in 3D, or 1 for
    thermal in any dimension (Lagrange isoparametric finite elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    numberOfComponents : int
        the number of components of the primal variable snapshots

    Returns
    -------
    scipy.sparse.csr_matrix
        the sparse matrix of the H10 scalar product
    """

    nbNodes = mesh.GetNumberOfNodes()
    dim     = mesh.GetDimensionality()

    ff = ElementFilter(mesh)
    ff.SetDimensionality(dim)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff, numberOfComponents)

    ev = []
    ei = []
    ej = []


    for name,data,ids in ff:
        p,w = LagrangeIsoParam[name]

        lenNumbering = len(numberings[0][name][0,:])
        ones = np.ones(lenNumbering,dtype=int)
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)
        for el in ids:

            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumberings = [numberings[j][name][el,:]+offset[j] for j in range(numberOfComponents)]

            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                BxByBzI = Jinv(space_ipvalues.valdphidxi[ip])

                for j in range(numberOfComponents):
                    ev.extend((w[ip]*Jdet)*np.tensordot(BxByBzI,BxByBzI, axes=(0,0)).ravel())

                    for i in leftNumberings[j]:
                        ei.extend(i*ones)
                        ej.extend(leftNumberings[j].ravel())

    mat = coo_matrix((ev, (ei,ej)), shape=(numberOfComponents*nbNodes,numberOfComponents*nbNodes)).tocsr()

    return mat

def ComputeInterpolationMatrix_FE_GaussPoint(mesh, feSpace, integrationRule,feNumbering=None,ipNumbering=None, elementFilter=None):
    from BasicTools.FE.Integration import IntegrateGeneral

    if elementFilter is None:
        elementFilter = ElementFilter(mesh)
        dim = mesh.GetDimensionality()
        elementFilter.SetDimensionality(dim)
    else:
        pass

    if feNumbering is None:
        numberingLeft = ComputeDofNumbering(mesh,Space=feSpace,elementFilter=elementFilter)
    else:
        numberingLeft  = feNumbering

    leftField = FEField(name="P1",numbering=numberingLeft,mesh=mesh,space=feSpace)

    from BasicTools.FE.Spaces.IPSpaces import GenerateSpaceForIntegrationPointInterpolation
    gaussSpace = GenerateSpaceForIntegrationPointInterpolation(integrationRule)
    if ipNumbering is None:
        numberingRight = ComputeDofNumbering(mesh,Space=gaussSpace,elementFilter=elementFilter)
    else:
        numberingRight = ipNumbering

    rightField = FEField(name="Gauss'",numbering=numberingRight,mesh=mesh,space=gaussSpace)

    from BasicTools.FE.SymWeakForm import GetField,GetTestField
    LF = GetField("P1",1)
    RF = GetTestField("Gauss",1)

    symForm = LF.T*RF
    interpMatrixMatrix,_ = IntegrateGeneral(mesh=mesh,constants={},fields=[],wform=symForm, unkownFields= [leftField],testFields=[rightField],onlyEvaluation=True,integrationRule=integrationRule,elementFilter=elementFilter)
    return interpMatrixMatrix

def ComputeJdetAtIntegPoint(mesh, elementSets = None, relativeDimension = 0):
    """
    Computes determinant of the Jacobian of the transformation of the
    transformation between the reference element and the mesh element, at
    the integration points. (Lagrange isoparametric finite elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    elementSets : list if strings
        sets of elements on which the det-Jacobians are computed
    relativeDimension : int (0, -1 or -2)
        difference between the dimension of the elements on which the function
        is computed and the dimensionality of the mesh

    Returns
    -------
    np.ndarray
        of size (NGauss,)
    """
    ff = ElementFilter(mesh)

    dimension = mesh.GetDimensionality() + relativeDimension
    ff.SetDimensionality(dimension)

    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for elementSet in elementSets:
            if elementSet:
                ff.AddTag(elementSet)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff)

    jDet = np.zeros(NGauss)

    countTotal = 0
    for name,data,ids in ff:

        p,w =  LagrangeIsoParam[name]
        NGaussperEl = len(w)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            space_ipvalues = spaces[name].SetIntegrationRule(p,w)
            for ip in range(NGaussperEl):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)

                jDet[countTotal] = Jdet

                countTotal += 1

    return jDet

def ComputePhiAtIntegPoint(mesh, elementSets = None, relativeDimension = 0):
    """
    Computes the value of the finite element shape functions at the integration
    points. (Lagrange isoparametric finite elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    elementSets : list if strings
        sets of elements on which the det-Jacobians are computed
    relativeDimension : int (0, -1 or -2)
        difference between the dimension of the elements on which the function
        is computed and the dimensionality of the mesh

    Returns
    -------
    np.ndarray
        of size (NGauss,)
    coo_matrix
        of size (NGauss, nbNodes)
    """
    ff = ElementFilter(mesh)

    dimension = mesh.GetDimensionality() + relativeDimension
    ff.SetDimensionality(dimension)

    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for elementSet in elementSets:
            if elementSet:
                ff.AddTag(elementSet)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff)

    nbNodes = mesh.GetNumberOfNodes()

    phiAtIntegPointIndices = []
    phiAtIntegPointValues = []
    row = []

    integrationWeights = np.zeros(NGauss)

    countTotal = 0
    for name,data,ids in ff:

        p,w =  LagrangeIsoParam[name]
        NGaussperEl = len(w)
        lenNumbering = len(numberings[0][name][0,:])

        ones = np.ones(lenNumbering,dtype=int)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            space_ipvalues = spaces[name].SetIntegrationRule(p,w)
            for ip in range(NGaussperEl):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)

                integrationWeights[countTotal] = w[ip]*Jdet

                leftNumbering = numberings[0][name][el,:]
                left = space_ipvalues.valN[ip]

                phiAtIntegPointIndices.extend(leftNumbering)
                phiAtIntegPointValues.extend(left)
                row.extend(countTotal*ones)
                countTotal += 1

    phiAtIntegPointMatrix = coo_matrix((phiAtIntegPointValues, (row, phiAtIntegPointIndices)), shape=(NGauss, nbNodes))

    return integrationWeights, phiAtIntegPointMatrix

def ComputePhiAtIntegPointFromElFilter(mesh, elFilter):
    """
    Computes the shape functions on the integration
    points and the integration weights associated with the integration scheme

    Parameters
    ----------
    mesh : UnstructuredMesh
    elFilter : elementFilter

    Returns
    -------
    np.ndarray
        integrationWeights

    coo_matrix of size (NGauss, nbNodes)
        phiAtIntegPointMatrix
    """

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)

    nbNodes = mesh.GetNumberOfNodes()

    phiAtIntegPointIndices = []
    phiAtIntegPointValues = []
    row = []
    integrationWeights = []

    countTotal = 0
    for name,data,ids in elFilter:

        p,w =  LagrangeIsoParam[name]
        NGaussperEl = len(w)
        lenNumbering = len(numbering[name][0,:])

        ones = np.ones(lenNumbering,dtype=int)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            space_ipvalues = LagrangeSpaceGeo[name].SetIntegrationRule(p,w)
            for ip in range(NGaussperEl):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)

                integrationWeights.append(w[ip]*Jdet)

                leftNumbering = numbering[name][el,:]
                left = space_ipvalues.valN[ip]

                phiAtIntegPointIndices.extend(leftNumbering)
                phiAtIntegPointValues.extend(left)
                row.extend(countTotal*ones)
                countTotal += 1

    NGauss = len(integrationWeights)
    integrationWeights = np.array(integrationWeights)

    phiAtIntegPointMatrix = coo_matrix((phiAtIntegPointValues, (row, phiAtIntegPointIndices)), shape=(NGauss, nbNodes))

    return integrationWeights, phiAtIntegPointMatrix

def ComputeGradPhiAtIntegPoint(mesh, elementSets = None, relativeDimension = 0):
    """
    Computes the components of the gradient of the shape functions on the
    integration points and the integration weights associated with the
    integration scheme

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    elementSets : list of strings
        element sets, support for the shape functions
    relativeDimension : int
        0, -1, or -2: the dimension of the element set relative to the dimension of the mesh

    Returns
    -------
    np.ndarray
        integrationWeights

    list (length dimensionality of the mesh) coo_matrix of size (NGauss, nbNodes)
        gradPhiAtIntegPoint
    """

    ff = ElementFilter(mesh)

    dimension = mesh.GetDimensionality() + relativeDimension
    ff.SetDimensionality(dimension)

    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for elementSet in elementSets:
            if elementSet:
                ff.AddTag(elementSet)


    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff)


    nbNodes = mesh.GetNumberOfNodes()
    meshDimension = mesh.GetDimensionality()

    GradPhiAtIntegPointIndices = []
    GradPhiAtIntegPointValues = [[] for i in range(meshDimension)]
    row = []

    integrationWeights = np.zeros(NGauss)

    countTotal = 0
    for name,data,ids in ff:

        p,w =  LagrangeIsoParam[name]
        NGaussperEl = len(w)
        lenNumbering = len(numberings[0][name][0,:])

        ones = np.ones(lenNumbering,dtype=int)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            space_ipvalues = spaces[name].SetIntegrationRule(p,w)
            for ip in range(NGaussperEl):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                BxByBzI = Jinv(space_ipvalues.valdphidxi[ip])

                integrationWeights[countTotal] = w[ip]*Jdet

                leftNumbering = numberings[0][name][el,:]

                GradPhiAtIntegPointIndices.extend(leftNumbering)
                for k in range(meshDimension):
                    GradPhiAtIntegPointValues[k].extend(BxByBzI[k])

                row.extend(countTotal*ones)
                countTotal += 1

    GradPhiAtIntegPointMatrix = [coo_matrix((GradPhiAtIntegPointValues[k], (row, GradPhiAtIntegPointIndices)), shape=(NGauss, nbNodes)) for k in range(meshDimension)]

    return integrationWeights, GradPhiAtIntegPointMatrix

def ComputeNormalsAtIntegPoint(mesh, elementSets):

    """
    Computes the normals at the elements from the sets elementSets in the
    ambiant space (i.e. if mesh is of dimensionality 3, elementSets should
    contain 2D elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    elementSets : list of strings
        sets of elements on which the function is computed

    Returns
    -------
    np.ndarray
        of size (dimensionality, NGauss)
    """
    ff = ElementFilter(mesh)

    dimension = mesh.GetDimensionality() - 1
    ff.SetDimensionality(dimension)


    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for elementSet in elementSets:
            if elementSet:
                ff.AddTag(elementSet)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff)

    normalsAtIntegPoint = np.empty((mesh.GetDimensionality(), NGauss))

    countTotal = 0
    for name,data,ids in ff:

        p,w =  LagrangeIsoParam[name]
        NGaussperEl = len(w)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            space_ipvalues = spaces[name].SetIntegrationRule(p,w)
            for ip in range(NGaussperEl):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                normal = space_ipvalues.GetNormal(Jack)
                normalsAtIntegPoint[:,countTotal] = normal
                countTotal += 1

    return normalsAtIntegPoint



def ComputeIntegrationPointsTags(mesh, dimension = None):
    """
    Returns a list of lists (of str) at each integration point (Lagrange
    isoparametric finite elements). Each list contains all the tags of the
    element of dimension "dimension" containg the considered integration points

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    elementSets : list if strings
        sets of elements on which the function is computed

    Returns
    -------
    list of lists (of str)
        of length numberOfIntegrationPoints
    """

    if dimension == None:
        dimension = mesh.GetDimensionality()

    ff = ElementFilter(mesh)
    ff.SetDimensionality(dimension)

    _, _, _, NGauss = PrepareFEComputation(mesh, ff, dimension)

    listOfTags = [[] for i in range(NGauss)]

    totalIntPointOffset = 0
    elementOffset = 0
    for name,data,ids in ff:
        _,w = LagrangeIsoParam[name]
        elNbeOfIntPoints = len(w)
        for tag in data.tags:
            for id in tag.GetIds():
                for intPoint in range((id - elementOffset)*elNbeOfIntPoints,(id - elementOffset + 1)*elNbeOfIntPoints):
                    listOfTags[intPoint].append(tag.name)

        totalIntPointOffset += elNbeOfIntPoints*data.GetNumberOfElements()
        elementOffset += data.GetNumberOfElements()

    return listOfTags


def CellDataToIntegrationPointsData(mesh, scalarFields, elementSet = None, relativeDimension = 0):
    """
    Change the representation of scalarFields from data constant by cell
    (elements) to data at integration points. (Lagrange isoparametric finite
    elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    scalarFields : np.ndarray of size (nbe of fields, nbe of elements) or dict
        with "nbe of fields" as keys and np.ndarray of size (nbe of elements,) or floats
        as values fields whose representation is changed by the function
    elementSet : elementSet defining the elements on which the function is
        computed. If None, takes all the elements of considered dimension
    relativeDimension : int (0, -1 or -2)
        difference between the dimension of the elements on which the function
        is computed and the dimensionality of the mesh

    Returns
    -------
    np.ndarray
        of size (nbe of fields, numberOfIntegrationPoints)
    """
    dimension = mesh.GetDimensionality() + relativeDimension

    ff = ElementFilter(mesh)
    ff.SetDimensionality(dimension)
    if elementSet != None:
        ff.AddTag(elementSet)

    _, _, _, NGauss = PrepareFEComputation(mesh, ff, dimension)

    if type(scalarFields) == dict:
        keymap = list(scalarFields.keys())
    elif type(scalarFields) == np.ndarray:
        keymap = np.arange(scalarFields.shape[0])
    else:
        raise("scalarFields is neither an np.ndarray nor a dict")

    numberOfFields = len(keymap)
    typeField = []
    for f in range(numberOfFields):
        if type(scalarFields[keymap[f]]) == float:
            typeField.append("scalar")
        else:
            typeField.append("vector")

    integrationPointsData = np.empty((numberOfFields, NGauss))

    countEl = 0
    countIp = 0
    for name,data,ids in ff:

        p,w = LagrangeIsoParam[name]
        numberOfIPperEl = len(w)

        for el in ids:

            for f in range(numberOfFields):
                if typeField[f] == "vector":
                    integrationPointsData[f,countIp:countIp+numberOfIPperEl] = scalarFields[keymap[f]][countEl]
                else:
                    integrationPointsData[f,countIp:countIp+numberOfIPperEl] = scalarFields[keymap[f]]

            countEl += 1
            countIp += numberOfIPperEl

    return integrationPointsData


def IntegrationPointsToCellData(mesh, scalarFields):
    """
    Change the representation of scalarFields from data at integration points
    to data constant by cell (elements) - taking the mean of values at the
    integration points in each elements. (Lagrange isoparametric finite
    elements)

    Parameters
    ----------
    mesh : UnstructuredMesh
        mesh on which the function is applied
    scalarFields : np.ndarray of size (nbe of fields, nbe of elements) or dict
        with "nbe of fields" as keys and np.ndarray of size (nbe of elements,) as
        values fields whose representation in changed by the function

    Returns
    -------
    np.ndarray
        of size (nbe of elements,)
    """

    if type(scalarFields) == dict:
        keymap = list(scalarFields.keys())
    elif type(scalarFields) == np.ndarray:
        keymap = np.arange(scalarFields.shape[0])
    else:
        raise("scalarFields is neither an np.ndarray nor a dict")

    numberOfFields = len(keymap)

    from BasicTools.FE.Fields import IPField as IPF

    cellData = []

    ef = ElementFilter(mesh=mesh,dimensionality =mesh.GetDimensionality())

    for f in range(numberOfFields):
        iPField = IPF.IPField("", mesh, rule=LagrangeIsoParam)
        iPField.Allocate()
        iPField.SetDataFromNumpy(scalarFields[keymap[f]], ef)
        cellData.append(iPField.GetCellRepresentation())

    return np.array(cellData)


def CheckIntegrity(GUI=False):
    from BasicTools.FE.SymPhysics import MecaPhysics

    mecaPhysics = MecaPhysics()
    wform = mecaPhysics.GetBulkFormulation(1.0,0.3)
    res,rhs = GetElementaryMatrixForFormulation(EN.Hexaedron_8,wform, unknownNames =mecaPhysics.GetPrimalNames() )

    import BasicTools.TestData as BasicToolsTestData
    from BasicTools.IO import GeofReader as GR
    mesh = GR.ReadGeof(BasicToolsTestData.GetTestDataPath()+"cube2.geof")
    ComputeL2ScalarProducMatrix(mesh, 1)
    ComputeL2ScalarProducMatrix(mesh, 3)
    ComputeH10ScalarProductMatrix(mesh, 1)
    ComputeH10ScalarProductMatrix(mesh, 3)
    ComputePhiAtIntegPoint(mesh)
    ComputeGradPhiAtIntegPoint(mesh)
    ComputeIntegrationPointsTags(mesh, 3)
    ComputeNormalsAtIntegPoint(mesh, ["x0"])
    length = len(mesh.elements["quad4"].tags["x0"].GetIds())
    fields = {'a': np.ones(length), 'b':1.}
    CellDataToIntegrationPointsData(mesh, fields, "x0", relativeDimension=-1)
    res = CellDataToIntegrationPointsData(mesh, np.ones((2,216)))
    nGauss = res.shape[1]
    IntegrationPointsToCellData(mesh, np.ones((2,nGauss)))

    from BasicTools.FE.IntegrationsRules import LagrangeP1
    mesh.elements["hex8"].tags.CreateTag("Transfert").SetIds([0,1] )
    elementFilter = ElementFilter(mesh,tag="Transfert")
    data = ComputeInterpolationMatrix_FE_GaussPoint(mesh=mesh,feSpace=LagrangeSpaceP1,integrationRule=LagrangeP1,elementFilter=elementFilter)

    print(np.dot(data.toarray(),np.arange(12)))

    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity(GUI=True))

    print("Done")

