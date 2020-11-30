# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

from BasicTools.FE.Integration import IntegrateGeneral
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.DofNumbering import ComputeDofNumbering
import BasicTools.Containers.ElementNames as EN

from scipy.sparse import coo_matrix, csr_matrix
from BasicTools.FE.IntegrationsRules import LagrangeIsoParam
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.Containers import Filters



def GetElementaryMatrixForFormulation(elemName, wform, unknownNames, space=LagrangeSpaceP1,geoFactor=None):
    # Explicitly specify signature to cleanly display default argument values
    # in sphinx autodoc generated documentation
    """
    GetElementaryMatrixForFormulation(elemName, wform, unknownNames, \
            space=BasicTools.FE.Spaces.FESpaces.LagrangeSpaceP1)
    """
    from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
    mesh = UnstructuredMesh()

    mesh.nodes = np.asarray(space[elemName].posN,dtype=float)
    if geoFactor is not None:
        mesh.nodes = mesh.nodes*geoFactor
    mesh.originalIDNodes = np.arange(0,mesh.GetNumberOfNodes(),dtype=np.int)

    elements = mesh.GetElementsOfType(elemName)
    elements.connectivity = np.arange(space[elemName].GetNumberOfShapeFunctions(),dtype=np.int)

    elements.connectivity.shape = (1,space[elemName].GetNumberOfShapeFunctions())
    elements.GetTag("3D").AddToTag(0)

    elements.originalIds = np.arange(0,1,dtype=np.int)
    elements.cpt = elements.connectivity.shape[0]

    mesh.PrepareForOutput()
    numbering = ComputeDofNumbering(mesh,space,)

    unkownFields = []
    for name in unknownNames:
        unkownFields.append(FEField(name,mesh,space,numbering))

    M,f = IntegrateGeneral(mesh=mesh,wform=wform, unkownFields= unkownFields,constants={},fields=[])

    return M.toarray(), f



def PrepareFEComputation(mesh, elementFilter = None, numberOfComponents = None):


    dim = mesh.GetDimensionality()
    if elementFilter == None:
        elementFilter = Filters.ElementFilter(mesh)
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



def ComputeL2ScalarProducMatrix(mesh, numberOfComponents):


    nbNodes = mesh.GetNumberOfNodes()
    dim     = mesh.GetDimensionality()

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dim)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff, numberOfComponents)

    ev = []
    ei = []
    ej = []


    for name,data,ids in ff:
        p,w =  LagrangeIsoParam[name]
        lenNumbering = len(numberings[0][name][0,:])
        #replace lenNumbering by nbsf = spaces[name].GetNumberOfShapeFunctions() ?
        ones = np.ones(lenNumbering,dtype=int)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            space_ipvalues = spaces[name].SetIntegrationRule(p,w)
            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                #Jack, Jdet, Jinv = spaces[name].GetJackAndDetI(ip,xcoor)
                for j in range(numberOfComponents):
                    leftNumbering = numberings[j][name][el,:] + offset[j]
                    left = space_ipvalues.valN[ip]
                    ev.extend(((w[ip]*Jdet)*np.outer(left, left)).ravel())
                    for i in leftNumbering:
                        ei.extend(i*ones)
                        ej.extend(leftNumbering.ravel())

    return coo_matrix((ev, (ei,ej)), shape=(numberOfComponents*nbNodes,numberOfComponents*nbNodes)).tocsr()



def ComputeH10ScalarProductMatrix(mesh, numberOfComponents):

    nbNodes = mesh.GetNumberOfNodes()
    dim     = mesh.GetDimensionality()

    ff = Filters.ElementFilter(mesh)
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

    if elementFilter is None:
        elementFilter = Filters.ElementFilter(mesh)
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



def ComputePhiAtIntegPoint(mesh, elementSets = None, relativeDimension = 0):
    """
    Computes the shape functions on the integration
    points and the integration weights associated with the integration scheme

    Parameters
    ----------
    mesh : UnstructuredMesh
    elementSets : list of strings
        element sets, support for the shape functions
    relativeDimension : int
        0, -1, or -2: the dimension of the element set relative to the dimension of the mesh


    Returns
    -------
    np.ndarray
        integrationWeights

    coo_matrix of size (NGauss, nbNodes)
        phiAtIntegPointMatrix
    """
    
    ff = Filters.ElementFilter(mesh)

    dimension = mesh.GetDimensionality() + relativeDimension
    ff.SetDimensionality(dimension)

    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for set in elementSets:
            if set:
                ff.AddTag(set)
                

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



def ComputeGradPhiAtIntegPoint(mesh, elementSets = None, relativeDimension = 0):
    """
    Computes the components of the gradient of the shape functions on the integration
    points and the integration weights associated with the integration scheme

    Parameters
    ----------
    mesh : UnstructuredMesh
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
    
    
    ff = Filters.ElementFilter(mesh)

    dimension = mesh.GetDimensionality() + relativeDimension
    ff.SetDimensionality(dimension)

    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for set in elementSets:
            if set:
                ff.AddTag(set)
                

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
    vector is the size of the number of nodes of "set"
    """
    
    ff = Filters.ElementFilter(mesh)

    dimension = mesh.GetDimensionality() - 1
    ff.SetDimensionality(dimension)


    if elementSets:
        assert type(elementSets) == list, "elementSets should be a list of elementSets"
        for set in elementSets:
            if set:
                ff.AddTag(set)


                
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

    
    if dimension == None:
        dimension = mesh.GetDimensionality()

    ff = Filters.ElementFilter(mesh)
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



def CellDataToIntegrationPointsData(mesh, set, scalarFields, relativeDimension = 0):
    """
    scalarFields is an np.ndarray of size (nbe of fields, number of dofs of the mesh)
    or a dict with "nbe of fields" keys
    """

    dimension = mesh.GetDimensionality() + relativeDimension

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension)
    ff.AddTag(set)

    _, _, _, NGauss = PrepareFEComputation(mesh, ff, dimension)
    
    if type(scalarFields) == dict:
        keymap = list(scalarFields.keys())
    elif type(scalarFields) == np.ndarray:
        keymap = np.arange(scalarFields.shape[0])
    else:
        raise("scalarFields is neither an np.ndarray nor a dict")
    
    numberOfFields = len(keymap)
    integrationPointsData = np.empty((numberOfFields, NGauss))

    countEl = 0
    countIp = 0
    for name,data,ids in ff:

        p,w = LagrangeIsoParam[name]
        numberOfIPperEl = len(w)

        for el in ids:
            
            for f in range(numberOfFields):
                integrationPointsData[f,countIp:countIp+numberOfIPperEl] = scalarFields[keymap[f]][countEl]
                
            countEl += 1
            countIp += numberOfIPperEl

    return integrationPointsData







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
    vector = np.ones(len(mesh.elements["quad4"].tags["x0"].GetIds()))
    ComputeNormalsAtIntegPoint(mesh, ["x0"])
    length = len(mesh.elements["quad4"].tags["x0"].GetIds())
    CellDataToIntegrationPointsData(mesh, "x0", np.ones(length))    

    #TODELETE
    #ComputeMecaIntegrator(mesh)
    #IntegrateVectorNormalComponentOnSurface(mesh, "x0", vector)
    #scalarFields = np.ones((3,mesh.GetNumberOfNodes()))
    #IntegrateOrderOneTensorOnSurface(mesh, "x0", scalarFields)
    #IntegrateOrderTwoTensorOnSurface(mesh, "x0", scalarFields)
    #IntegrateCentrifugalEffect(mesh, {'ALLELEMENT':1.}, np.array([1.,0.,0.]), np.array([0.,0.,0.]))    

    from BasicTools.FE.IntegrationsRules import LagrangeP1
    mesh.elements["hex8"].tags.CreateTag("Transfert").SetIds([0,1] )
    elementFilter = Filters.ElementFilter(mesh,tag="Transfert")
    data = ComputeInterpolationMatrix_FE_GaussPoint(mesh=mesh,feSpace=LagrangeSpaceP1,integrationRule=LagrangeP1,elementFilter=elementFilter)

    print(np.dot(data.toarray(),np.arange(12)))

    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity(GUI=True))

    print("Done")






"""
#TODELETE
def ComputeMecaIntegrator(mesh, elementSet = "ALLELEMENT"):
    #elementSet must tag element of same dimensionality as mesh

    nbNodes = mesh.GetNumberOfNodes()
    dimension = mesh.GetDimensionality()

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension)
    if elementSet != "ALLELEMENT":
        ff.AddTag(elementSet)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff, dimension)

    integrationWeights = np.zeros(NGauss)

    #convInd = {1:[0], 2:[0, 2, 2, 1], 3:[0, 3, 4, 3, 1, 5, 4, 5, 2]}
    convInd = {1:[0], 2:[0, 2, 2, 1], 3:[0, 3, 5, 3, 1, 4, 5, 4, 2]}
    nbeInd  = {1:1, 2:3, 3:6}

    row = []
    col = []
    dat = []

    count = 0
    for name,data,ids in ff:
        p,w = LagrangeIsoParam[name]

        nbsf = spaces[name].GetNumberOfShapeFunctions()
        ones = np.ones(dimension*nbsf,dtype=int)

        B = np.zeros((dimension*nbsf, nbeInd[dimension]), dtype=np.float)
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        for el in ids:

            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumbering = np.concatenate([numberings[j][name][el,:]+offset[j] for j in range(dimension)])

            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                BxByBzI = Jinv(space_ipvalues.valdphidxi[ip])

                integrationWeights[count] = w[ip]*Jdet

                for i in range(nbsf):
                    for k in range(dimension):
                        for l in range(dimension):
                            B[i+l*nbsf, convInd[dimension][k*dimension + l]] = BxByBzI[k,i]

                dat.extend((B.T).ravel())
                for i in range(nbeInd[dimension]):
                    row.extend(leftNumbering.ravel())
                    col.extend(ones*(nbeInd[dimension]*count+i))

                count += 1

    dat = np.array(dat)
    row = np.array(row)
    col = np.array(col)

    integrator = coo_matrix((dat, (row, col)), shape=(dimension*nbNodes,nbeInd[dimension]*NGauss)).tocsr()

    return integrationWeights, integrator"""



"""
#TODELETE
def IntegrateVectorNormalComponentOnSurface(mesh, set, vector):
    '''
    vector is the size of the number of nodes of "set"
    '''

    nbNodes = mesh.GetNumberOfNodes()
    dimension = mesh.GetDimensionality()

    res = np.zeros(dimension*nbNodes)


    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension-1)
    ff.AddTag(set)

    spaces, numberings, offset, _ = PrepareFEComputation(mesh, ff, dimension)

    count = 0
    for name,data,ids in ff:

        p,w = LagrangeIsoParam[name]
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        for el in ids:
            pressureValue = vector[count]; count += 1
            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumbering = np.concatenate([numberings[j][name][el,:]+offset[j] for j in range(dimension)])

            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                normal = space_ipvalues.GetNormal(Jack)
                left = space_ipvalues.valN[ip]
                cartesian_product = np.dot(normal.reshape((normal.shape[0],1)),left.reshape((1,left.shape[0])))

                res[leftNumbering] += ((w[ip]*Jdet*pressureValue)*cartesian_product).ravel()

    return res
"""



"""
#TODELETE
def IntegrateOrderOneTensorOnSurface(mesh, set, scalarFields):
    '''
    scalarFields is of size (nbe of fields, number of dofs of the mesh)
    '''


    dimension = mesh.GetDimensionality()

    nFields = scalarFields.shape[0]

    res = np.zeros(nFields)

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension-1)
    ff.AddTag(set)

    spaces, numberings, offset, _ = PrepareFEComputation(mesh, ff, dimension)

    for name,data,ids in ff:
        p,w = LagrangeIsoParam[name]
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        for el in ids:

            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumbering = numberings[0][name][el,:] + offset[0]

            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                left = space_ipvalues.valN[ip]

                res += w[ip]*Jdet*np.dot(scalarFields[:,leftNumbering],left)
                #contrib = (w[ip]*Jdet)
                #for i in range(nFields):
                #  res[i] += contrib*np.dot(scalarFields[i][leftNumbering],left)

    return res
"""



"""
#TODELETE
def IntegrateOrderTwoTensorOnSurface(mesh, set, scalarFields):
    '''
    scalarFields is of size (nbe of fields, number of dofs of the mesh)
    '''


    dimension = mesh.GetDimensionality()

    nFields = scalarFields.shape[0]

    res = np.zeros((nFields,nFields))


    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension-1)
    ff.AddTag(set)

    spaces, numberings, offset, _ = PrepareFEComputation(mesh, ff, dimension)

    for name,data,ids in ff:
        p,w = LagrangeIsoParam[name]
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        for el in ids:

            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumbering = numberings[0][name][el,:] + offset[0]

            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                left = space_ipvalues.valN[ip]

                tempt = np.dot(scalarFields[:,leftNumbering],left)
                res += w[ip]*Jdet*np.outer(tempt, tempt)

                '''contrib = w[ip]*Jdet
                for i in range(nFields):
                    contribI = contrib*np.dot(scalarFields[i][leftNumbering],left)
                    for j in range(nFields):
                        res[i,j] += contribI*np.dot(scalarFields[j][leftNumbering],left)'''

    return res
"""


"""
#TODELETE
def IntegrateCentrifugalEffect(mesh, density, rotationAxis, rotationCenter):

    '''
    density: dictionary with keys: elementSet and value: float
    '''

    nbNodes = mesh.GetNumberOfNodes()
    dimension = mesh.GetDimensionality()

    res = np.zeros(dimension*nbNodes)

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff, dimension)

    _, phiAtIntegPoint = ComputePhiAtIntegPoint(mesh)

    integrationPointsPosition = phiAtIntegPoint.dot(mesh.nodes)

    densityTags = set(list(density.keys()))

    count = 0
    for name,data,ids in ff:

        elementTags = {}
        for el in ids:
            elementTags[el] = []
        for tag in densityTags:
            for el in mesh.GetElementsInTag(tag):
                elementTags[el].append(tag)

        p,w =  LagrangeIsoParam[name]
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        for el in ids:

            localTags = elementTags[el]
            if not localTags:
                localTags.append("ALLELEMENT")
            if len(localTags) > 1:
                raise("more than one behavior law associate with element "+str(el))
            localDensity = density[localTags[0]]

            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumbering = np.concatenate([numberings[j][name][el,:]+offset[j] for j in range(dimension)])

            for ip in range(len(w)):
                ipPositionFromRotationCenter = integrationPointsPosition[count,:] - rotationCenter

                count += 1

                length = np.vdot(ipPositionFromRotationCenter,rotationAxis)
                ipProjectionFromRotationCenter = length*rotationAxis
                r = ipPositionFromRotationCenter - ipProjectionFromRotationCenter

                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                left = space_ipvalues.valN[ip]
                cartesian_product = np.dot(r.reshape((r.shape[0],1)),left.reshape((1,left.shape[0])))

                res[leftNumbering] += localDensity*((w[ip]*Jdet)*cartesian_product).ravel()

    return res
"""




"""
#TODELETE
def ComputeFEInterpGradMatAtGaussPoint(mesh):

    nbNodes = mesh.GetNumberOfNodes()
    dim = mesh.GetDimensionality()

    spaces, _, _, NGauss0 = PrepareFEComputation(mesh)

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dim)

    FEInterpAtIntegPointIndices = []
    FEInterpAtIntegPointGradMatrix = [[] for i in range(dim)]
    row = []

    integrationWeights = np.zeros(NGauss0)

    countElementType = 0
    countTotal = 0
    for name,data,ids in ff:

        NnodeperEl = EN.numberOfNodes[name]
        p,w = LagrangeIsoParam[name]
        NGaussperEl = len(w)
        NGauss = data.GetNumberOfElements()*NGaussperEl
        nbElements = data.GetNumberOfElements()

        FEInterpAtIntegPointIndices.append(np.zeros((NGauss,NnodeperEl),dtype=np.int32))
        for i in range(dim):
            FEInterpAtIntegPointGradMatrix[i].append(np.zeros((NGauss,NnodeperEl)))
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        count = 0
        for i in range(nbElements):
            xcoor = np.array([mesh.nodes[data.connectivity[i,j]] for j in range(NnodeperEl)])
            for j in range(NGaussperEl):

                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(j,xcoor)
                BxByBzI = Jinv(space_ipvalues.valdphidxi[j])

                FEInterpAtIntegPointIndices[countElementType][count,:] = data.connectivity[i,:]
                for k in range(dim):
                    FEInterpAtIntegPointGradMatrix[k][countElementType][count,:] = BxByBzI[k]

                integrationWeights[countTotal] = w[j]*Jdet

                count += 1
                countTotal += 1

        row.append(np.concatenate([[i for j in range(NnodeperEl)] for i in range(NGauss)]))

        countElementType += 1

    FEInterpAtIntegPointIndices = np.concatenate([ind.flatten() for ind in FEInterpAtIntegPointIndices])
    FEInterpAtIntegPointGradMatrix  = [np.concatenate([mat.flatten() for mat in FEInterpAtIntegPointGradMatrix[k]]) for k in range(dim)]

    row = np.concatenate(row)

    FEInterpAtIntegPointGradMatrix = [coo_matrix((FEInterpAtIntegPointGradMatrix[k], (row, FEInterpAtIntegPointIndices)), shape=(NGauss0, nbNodes)) for k in range(dim)]


    return integrationWeights, FEInterpAtIntegPointGradMatrix
"""



"""
#TODELETE
def ComputeFEInterpMatAtGaussPoint(mesh):

    nbNodes = mesh.GetNumberOfNodes()
    dim = mesh.GetDimensionality()

    spaces, _, _, NGauss0 = PrepareFEComputation(mesh)

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dim)

    FEInterpAtIntegPointIndices = []
    FEInterpAtIntegPointMatrix = []
    row = []

    countElementType = 0
    for name,data,ids in ff:

        NnodeperEl = EN.numberOfNodes[name]
        p,w = LagrangeIsoParam[name]
        NGaussperEl = len(w)
        NGauss = data.GetNumberOfElements()*NGaussperEl
        nbElements = data.GetNumberOfElements()

        FEInterpAtIntegPointIndices.append(np.zeros((NGauss,NnodeperEl),dtype=np.int32))
        FEInterpAtIntegPointMatrix.append(np.zeros((NGauss,NnodeperEl)))
        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        count = 0
        for i in range(nbElements):
          xcoor = np.array([mesh.nodes[data.connectivity[i,j]] for j in range(NnodeperEl)])
          for j in range(NGaussperEl):

            Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(j,xcoor)

            FEInterpAtIntegPointIndices[countElementType][count,:] = data.connectivity[i,:]
            FEInterpAtIntegPointMatrix[countElementType][count,:] = space_ipvalues.valN[j]

            count += 1

        row.append(np.concatenate([[i for j in range(NnodeperEl)] for i in range(NGauss)]))

        countElementType += 1

    FEInterpAtIntegPointIndices = np.concatenate([ind.flatten() for ind in FEInterpAtIntegPointIndices])
    FEInterpAtIntegPointMatrix = np.concatenate([ind.flatten() for ind in FEInterpAtIntegPointMatrix])

    row = np.concatenate(row)

    FEInterpAtIntegPointMatrix = coo_matrix((FEInterpAtIntegPointMatrix, (row, FEInterpAtIntegPointIndices)), shape=(NGauss0, nbNodes))

    return FEInterpAtIntegPointMatrix
"""

"""
#TODELETE
def ComputeNodeToIntegThermalOperatorRadiation(mesh, elementSets):
    #elementSet must tag element of one less dimensionality than mesh

    nbNodes = mesh.GetNumberOfNodes()
    dimension = mesh.GetDimensionality()

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension-1)
    for set in elementSets:
        ff.AddTag(set)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff, dimension-1)

    integrationWeightsRadiation = np.zeros(NGauss)

    if len(elementSets) > 0:

        row = []
        col = []
        dat = []

        count = 0
        for name,data,ids in ff:
            p,w = LagrangeIsoParam[name]

            lenNumbering = len(numberings[0][name][0,:])
            ones = np.ones(lenNumbering,dtype=int)
            space_ipvalues = spaces[name].SetIntegrationRule(p,w)

            for el in ids:
                xcoor = mesh.nodes[data.connectivity[el],:]
                leftNumbering = numberings[0][name][el,:] + offset[0]

                for ip in range(len(w)):
                    Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                    integrationWeightsRadiation[count] = w[ip]*Jdet

                    left = space_ipvalues.valN[ip]
                    dat.extend(left.ravel())

                    row.extend(leftNumbering.ravel())
                    col.extend(ones*count)

                    count += 1

        dat = np.array(dat)
        row = np.array(row)
        col = np.array(col)

        nodeToIntegThermalRadiationOperator = coo_matrix((dat, (row, col)), shape=(nbNodes,NGauss)).tocsr()

    else:

        nodeToIntegThermalRadiationOperator = csr_matrix((nbNodes, NGauss), dtype=float)

    return integrationWeightsRadiation, nodeToIntegThermalRadiationOperator
"""

"""
#TODELETE
def ComputeNodeToIntegThermalOperator(mesh):

    nbNodes = mesh.GetNumberOfNodes()
    dimension = mesh.GetDimensionality()

    ff = Filters.ElementFilter(mesh)
    ff.SetDimensionality(dimension)

    spaces, numberings, offset, NGauss = PrepareFEComputation(mesh, ff, dimension)

    integrationWeights = np.zeros(NGauss)

    row = []
    col = []
    dat = []

    row2 = []
    col2 = []
    dat2 = []

    count = 0
    for name,data,ids in ff:
        p,w = LagrangeIsoParam[name]

        lenNumbering = len(numberings[0][name][0,:])
        ones = np.ones(lenNumbering,dtype=int)

        space_ipvalues = spaces[name].SetIntegrationRule(p,w)

        for el in ids:
            xcoor = mesh.nodes[data.connectivity[el],:]
            leftNumbering = numberings[0][name][el,:]

            for ip in range(len(w)):
                Jack, Jdet, Jinv = space_ipvalues.GetJackAndDetI(ip,xcoor)
                BxByBzI = Jinv(space_ipvalues.valdphidxi[ip])

                integrationWeights[count] = w[ip]*Jdet

                left = space_ipvalues.valN[ip]
                #prod = np.tensordot(BxByBzI,BxByBzI, axes=(0,0))

                dat.extend(left.ravel())
                row.extend(leftNumbering.ravel())
                col.extend(ones*count)

                for i in range(dimension):
                    dat2.extend(BxByBzI[i,:])
                    row2.extend(leftNumbering)
                    col2.extend(ones*(dimension*count+i))

                count += 1

    dat = np.array(dat)
    row = np.array(row)
    col = np.array(col)

    dat2 = np.array(dat2)
    row2 = np.array(row2)
    col2 = np.array(col2)

    nodeToIntegThermalOperator = coo_matrix((dat, (row, col)), shape=(nbNodes,NGauss)).tocsr()
    nodeToIntegGradThermalOperator = coo_matrix((dat2, (row2, col2)), shape=(nbNodes,NGauss*dimension)).tocsr()


    return integrationWeights, nodeToIntegThermalOperator, nodeToIntegGradThermalOperator
"""