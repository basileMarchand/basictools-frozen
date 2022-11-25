# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np
import warnings
import networkx
from scipy.sparse.linalg import eigsh
from scipy import sparse, optimize

from BasicTools.Containers.UnstructuredMeshInspectionTools import GetDualGraphNodeToElement, GetDualGraph
import BasicTools.Containers.UnstructuredMeshModificationTools as UMMT
import BasicTools.Containers.UnstructuredMeshInspectionTools as UMIT
from BasicTools.Containers import UnstructuredMesh
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType
import BasicTools.Containers.ElementNames as ElementNames
import BasicTools.Containers.Filters as Filters
import BasicTools.Containers.Tags as T



def PartitionMesh(inmesh, nbSubdomains, driver="Metis"): # TBD: add other drivers
    '''Generates a subdomain vector providing labels ordered for each element of the input mesh
    based on a third party graph partitionning solver. Current working options are: Metis'''

    elemGraph, usedCells = GetElementGraph(inmesh)

    if driver=="Metis":
        import pymetis
        # pymetis cannot interprete -1s is nparrays --> transforming into a list of list
        meshGraph = [np.delete(np.asarray(line, dtype=int), np.where(line==-1)) for line in elemGraph]
        cuts, subdomainMap = pymetis.part_graph(nbSubdomains, adjacency=meshGraph)
    elif driver=="Scotch":
        raise ValueError("Scotch not implemented yet")
    elif driver=="Chaco":
        raise ValueError("Chaco not implemented yet")
    else:
        raise ValueError("""The driver option you have selected is not implemented or incorrectly spelled!
            Choose between Scotch, Metis, or Chaco.""")

    #elementsBySubdomain = [np.argwhere(np.array(maskVector) == sd).ravel() for sd in [i for i in range(nbSubdomains)]]
    return subdomainMap

def GetElementGraph(inmesh, maxNumConnections = 200): # TBD: Test on 1M elts mesh, make use of dimensionality
    '''Generates the graph connecting elements through faces (in 3D) or edges (in 2D).
    Also provides the array of cells used.'''

    meshGraph, dump = GetDualGraphNodeToElement(inmesh) # build node connectivity table first
    elemGraph = np.zeros((inmesh.GetNumberOfElements(), maxNumConnections), dtype = PBasicIndexType) - 1
    usedCells = np.zeros(inmesh.GetNumberOfElements(), dtype = PBasicIndexType)
    maxDim = max([ElementNames.dimension[groupName] for groupName, elemGroup in inmesh.elements.items()]) # check largest dimensionality
    if min([ElementNames.dimension[groupName] for groupName, elemGroup in inmesh.elements.items()]) < maxDim:
        warnings.warn("Input mesh has more than one dimensionality: Only the largest will be used to build element graph.")

    for groupName, elemGroup, ids in Filters.ElementFilter(inmesh,dimensionality=maxDim): # for each element group of largest dimensionality
        neighboringElems = np.zeros((elemGroup.connectivity.shape[1], meshGraph.shape[1]), dtype = PBasicIndexType)
        nbNodesPerFace = []
        for face in ElementNames.faces[groupName]: # evaluate face connectivity
            nbNodesPerFace.append(len(face[1]))
        for element in ids: # for each element in the filtered groups
            elementNodes = elemGroup.connectivity[element,:]
            neighboringElems.fill(0)
            for j, node in enumerate(elementNodes):
                neighboringElems[j] = meshGraph[node]
            unique, counts = np.unique(neighboringElems, return_counts=True) # count how many time each other element is connected to a node
            yindex = 0
            for j in range(len(unique)):
                if counts[j] in nbNodesPerFace: # if nb of connections correspond to face connectivity, add it to graph
                    elemGraph[element,yindex] = unique[j]
                    usedCells[element] += 1
                    yindex += 1

    maxsize = np.max(np.sum(elemGraph>=0,axis=1)) # crop output data
    elemGraph = elemGraph[:,0:maxsize]
    return elemGraph, usedCells



def InitializeGraphPointsFromMeshPoints(inMesh):
    '''Initializes a networkx graph with nodes consistant with the number of nodes of an UnstructuredMesh.
    This enables further edge addition compatible with the connectivity of the elements of the UnstructuredMesh.

    Parameters
    ----------
    inMesh : UnstructuredMesh
        input mesh

    Returns
    -------
    networkx.Graph
        initialized graph
    '''
    G = networkx.Graph()
    for i in range(inMesh.GetNumberOfNodes()):
        G.add_node(i)
    return G


def InitializeGraphPointsFromMeshElements(inMesh, dimensionality = None):
    '''Initializes a networkx graph with nodes consistant with the number of elements of an UnstructuredMesh.
    This enables further edge addition compatible with a chosen global numbering of the elements of the UnstructuredMesh.

    Parameters
    ----------
    inMesh : UnstructuredMesh
        input mesh
    dimensionality : int
        dimension of the elements considered to initalize the graph

    Returns
    -------
    networkx.Graph
        initialized graph
    '''
    if dimensionality == None:
        dimensionality = inMesh.GetDimensionality()

    G = networkx.Graph()

    elFilter = Filters.ElementFilter(inMesh, dimensionality = dimensionality)
    count = 0
    for _, _, ids in elFilter:
        for _ in range(len(ids)):
            G.add_node(count)
            count += 1
    return G


def GetExtendedElementGraph(inMesh, dimensionality = None):
    '''Creates a networkx graph from the element connectivity on an UnstrucutredMesh in the following sense:
    an element is linked to another in the graph if they share a vertex.

    Parameters
    ----------
    inMesh : UnstructuredMesh
        input mesh
    dimensionality : int
        dimension of the elements considered to initalize the graph

    Returns
    -------
    networkx.Graph
        initialized graph
    '''
    if dimensionality == None:
        dimensionality = inMesh.GetDimensionality()

    G = InitializeGraphPointsFromMeshElements(inMesh, dimensionality = dimensionality)

    elFilter = Filters.ElementFilter(inMesh, dimensionality = dimensionality)
    for _, data, ids in elFilter:
        for idd in range(len(ids)):
            inMesh.GetNodalTag("current").SetIds(data.connectivity[idd])
            ff = Filters.ElementFilter(inMesh, dimensionality = dimensionality, nTags=["current"], nTagsTreatment="leastonenode")
            for _, _, ids2 in ff:
                for el in ids2:
                    G.add_edge(idd, el)
    return G



def ComputeMeshLaplacianEigenmaps(mesh, dimensionality = None, nEigenmaps = 10, distFunc = None, normalizedLaplacian = False):
    '''Computes the Laplacian engenmaps of a mesh

    Parameters
    ----------
    mesh : UnstructuredMesh
        input mesh
    dimensionality : int
        dimension of the elements considered to initalize the graph
    nEigenmaps : int
        number of computed eigenmaps (less or equal to the number of nodes of mesh)
    distFunc : func
        function applied to the lengh of the edges of the mesh, and attached of the
        corresponding edge of the graph of the mesh
    normalizedLaplacian : bool
        if "True", the normalized Laplacian matrix is taken

    Returns
    -------
    ndarray(nEigenmaps)
        eigenvalues
    ndarray(numberOfNodes, nEigenmaps)
        Laplacian eigenmaps
    '''

    if dimensionality == None:
        dimensionality = mesh.GetDimensionality()

    elFilter = Filters.ElementFilter(mesh, dimensionality = dimensionality)

    if distFunc == None:
        def distFunc(x):
            return x

    edges = {}
    for name, data ,ids in elFilter:
        for idd in ids:
            for face in ElementNames.faces[name]:
                edge = np.sort(data.connectivity[idd][face[1]])
                length = np.linalg.norm(mesh.nodes[edge[1]]-mesh.nodes[edge[0]])
                edges[(edge[0], edge[1])] = distFunc(length)


    G = InitializeGraphPointsFromMeshPoints(mesh)
    for edge, length in edges.items():
        G.add_edge(edge[0], edge[1], weight = length)

    if normalizedLaplacian == False:
        L = networkx.laplacian_matrix(G).asfptype()
    else:
        L = networkx.normalized_laplacian_matrix(G).asfptype()

    if nEigenmaps > L.shape[0] - 1:
        nEigenmaps = L.shape[0] - 1
        print("reducing number of eigenmaps to mesh number of nodes minus 1")

    w, v = eigsh(L, k = nEigenmaps, which = "SM")

    return w, v


def RenumberMeshForParametrization(inMesh, inPlace = True, boundaryOrientation = "direct", fixedBoundaryPoints = None):
    """
    Only for linear triangle meshes
    Renumbering: the last nodeIDs correspond to the boundary, ordered


    Parameters
    ----------
    inMesh : UnstructuredMesh
        input triangular to be renumbered
    inPlace : bool
        if "True", inMesh is modified
        if "False", inMesh is let unmodified, and a new mesh is produced
    boundaryOrientation : str
        if "direct, the boundary of the parametrisation is constructed in the direct trigonometric order
        if "indirect", the boundary of the parametrisation is constructed in the indirect trigonometric orderc order
    fixedBoundaryPoints : list
        list containing lists of two np.ndarrays. Each 2-member list is used to identify one
    point on the boundary: the first array contains the specified components, and the second the

    Returns
    -------
    UnstructuredMesh
        renumbered mesh
    ndarray(1) of ints
        renumbering of the nodes of the returned renumbered mesh, with respect to inMesh
    int
        number of node of the boundary of inMesh

    """
    # assert mesh of linear triangles
    for name, data in inMesh.elements.items():
        name == ElementNames.Triangle_3

    if inPlace == True:
        mesh = inMesh
    else:
        import copy
        mesh = copy.deepcopy(inMesh)

    # Retrieve the elements of the line boundary
    skin = UMMT.ComputeSkin(mesh, md = 2)
    skin.ComputeBoundingBox()


    # Create a path linking nodes of the line boundary, starting with the node with smallest coordinates
    # and going in the direction increasing the value of the second coordinate the least

    bars = skin.elements[ElementNames.Bar_2].connectivity
    nodeGraph = UMIT.GetDualGraph(skin)[0]
    assert nodeGraph.shape[1] == 2, "skin is not a mesh of bars"

    indicesBars = np.sort(np.unique(bars.flatten()))

    if fixedBoundaryPoints == None:

        indicesNodesXmin = inMesh.nodes[indicesBars,0] == inMesh.nodes[np.argmin(inMesh.nodes[indicesBars,0]),0]
        nodesXmin = inMesh.nodes[indicesBars[indicesNodesXmin], :]

        indicesNodesmin = nodesXmin[:,1] == nodesXmin[np.argmin(nodesXmin[:,1]),1]
        nodesmin = nodesXmin[indicesNodesmin, :]

        if inMesh.GetDimensionality() == 3:
            indicesNodesmin = nodesmin[:,2] == nodesmin[np.argmin(nodesmin[:,2]),2]
            nodesmin = nodesmin[indicesNodesmin, :]

        indexInBars = np.where((inMesh.nodes[indicesBars,:] == nodesmin).all(axis=1))[0]
        assert indexInBars.shape == (1,)
        indexInBars = indexInBars[0]
        assert (inMesh.nodes[indicesBars[indexInBars],:] == nodesmin).all()

        pMin = indicesBars[indexInBars]
        print("starting walking along line boundary at point... =", str(nodesmin))

    else:
        inds, point = fixedBoundaryPoints[0][0], fixedBoundaryPoints[0][1]
        indexInBars = (np.linalg.norm(np.subtract(inMesh.nodes[indicesBars,:][:,inds], point), axis = 1)).argmin()
        pMin = indicesBars[indexInBars]
        print("starting walking along line boundary at point... =", str(inMesh.nodes[pMin,:]))


    p1 = p1init = pMin
    p2_candidate = [nodeGraph[pMin][0], nodeGraph[pMin][1]]

    if fixedBoundaryPoints == None:
        # choose direction
        p2 = p2_candidate[np.argmin(np.array([inMesh.nodes[p2_candidate[0],1], inMesh.nodes[p2_candidate[1],1]]))]

    else:
        # choose direction from second point set on boundary
        inds = fixedBoundaryPoints[1][0]
        delta_fixedBoundaryPoints = fixedBoundaryPoints[1][1] - fixedBoundaryPoints[0][1]
        delta_fixedBoundaryPoints /= np.linalg.norm(delta_fixedBoundaryPoints)

        delta_candidate = np.array([inMesh.nodes[p2c,inds] - inMesh.nodes[pMin,inds] for p2c in p2_candidate])
        delta_candidate[0] /= np.linalg.norm(delta_candidate[0])
        delta_candidate[1] /= np.linalg.norm(delta_candidate[1])

        error_delta_candidate = []
        error_delta_candidate.append(np.subtract(delta_candidate[0], delta_fixedBoundaryPoints))
        error_delta_candidate.append(np.subtract(delta_candidate[1], delta_fixedBoundaryPoints))

        p2 = p2_candidate[np.linalg.norm(error_delta_candidate, axis = 1).argmin()]

    print("... walking toward point =", str(inMesh.nodes[p2,:]))

    path = [p1, p2]
    while p2 != p1init:
        p2save = p2
        p2 = nodeGraph[p2][nodeGraph[p2]!=p1][0]
        p1 = p2save
        path.append(p2)
    path = path[:-1]

    if boundaryOrientation == "indirect":
        path = path[::-1]

    # Renumber the node, keeping at the end the continuous path along the line boundary
    N = mesh.GetNumberOfNodes()
    nBoundary = len(path)

    initOrder = np.arange(N)
    interiorNumberings = np.delete(initOrder, path)

    renumb = np.hstack((interiorNumberings, path))

    assert len(renumb) == N

    invRenumb = np.argsort(renumb)

    mesh.nodes = mesh.nodes[renumb,:]
    for _, data in mesh.elements.items():
        data.connectivity = invRenumb[data.connectivity]
    mesh.ConvertDataForNativeTreatment()

    return mesh, renumb, nBoundary



def FloaterMeshParametrization(inMesh, nBoundary, outShape = "circle", boundaryOrientation = "direct", fixedInteriorPoints = None, fixedBoundaryPoints = None):
    """
    STILL LARGELY EXPERIMENTAL

    Only for linear triangular meshes

    Computes the Mesh Parametrization algorithm [1] proposed by Floater,
    in the case of target parametrization fitted to the unit 2D circle (R=1) or square (L=1).
    Adapted for ML need: the outShape's boundary is sampled following the curvilinear abscissa along
    the boundary on inMesh (only for outShape = "circle" for the moment)

    Parameters
    ----------
    inMesh : UnstructuredMesh
        Renumbered triangular mesh to parametrize
    nBoundary : int
        number nodes on the line boundary
    outShape : str
        if "circle", the boundary of inMesh is mapped into the unit circle
        if "square", the boundary of inMesh is mapped into the unit square
    boundaryOrientation : str
        if "direct, the boundary of the parametrisation is constructed in the direct trigonometric order
        if "indirect", the boundary of the parametrisation is constructed in the indirect trigonometric order
    fixedInteriorPoints : dict
        with one key, and corresponding value, a list: [ndarray(n), ndarray(n,2)],
        with n the number of interior points to be fixed; the first ndarray is the index of the considered
        interior point, the second ndarray is the corresponding prescribed positions
        if key is "mean", the interior points are displaced by the mean of the prescribed positions
        if key is "value", the interior points are displaced by the value of the prescribed positions
    fixedBoundaryPoints: list
        list of lists: [ndarray(2), ndarray(2)], helping definining a point in inMesh; the first ndarray is the component
        of a point on the boundary, and the second array is the value of corresponding component. Tested for triangular meshes
        in the 3D space.

    Returns
    -------
    UnstructuredMesh
        parametrization of mesh
    dict
        containing two keys: "minEdge" and "maxEdge", with values floats containing the minimal
        and maximal edged length of the parametrized mesh

    Notes
    -----
        mesh mush be a renumbered UnstructuredMesh of triangles (either in
        a 2D or 3D ambiant space), with a line boundary (no closed surface in 3D).
        outShape = "circle" is more robust in the sense that is inMesh has a 2D square-like,
        for triangles may ended up flat with  outShape = "square"

    References
    ----------
        [1] M. S. Floater. Parametrization and smooth approximation of surface
        triangulations, 1997. URL: https://www.sciencedirect.com/science/article/abs/pii/S0167839696000313
    """
    import copy
    mesh = copy.deepcopy(inMesh)

    N = mesh.GetNumberOfNodes()
    n = N - nBoundary

    u = np.zeros((mesh.nodes.shape[0],2))

    if outShape == "square":
        print("!!! Warning, the implmentation outShape == 'square' is *very* experimental !!!")
        if boundaryOrientation == "indirect":
            raise NotImplementedError("Cannot use 'square' outShape with 'indirect' boundaryOrientation")
        if fixedInteriorPoints != None:
            raise NotImplementedError("Cannot use 'square' outShape with fixedInteriorPoints not None")
        if fixedBoundaryPoints != None:
            raise NotImplementedError("Cannot use 'square' outShape with fixedBoundaryPoints not None")

        # Set the boundary on the parametrization on the unit square
        L = nBoundary//4
        r = nBoundary%4

        u[n:n+L,0] = np.linspace(1/L,1,L)
        u[n:n+L,1] = 0.
        u[n+L:n+2*L,0] = 1.
        u[n+L:n+2*L,1] = np.linspace(1/L,1,L)
        u[n+2*L:n+3*L,0] = np.linspace(1-1/L,0,L)
        u[n+2*L:n+3*L,1] = 1.
        u[n+3*L:n+4*L+r,0] = 0.
        u[n+3*L:n+4*L+r,1] = np.linspace(1-1/(L+r),0,(L+r))

    elif outShape == "circle":
        # Set the boundary on the parametrization on the unit circle

        lengthAlongBoundary = [0]
        cumulativeLength = 0.
        indices = np.arange(n+1, N)
        for i in indices:
            p1 = mesh.nodes[i-1,:]
            p2 = mesh.nodes[i,:]
            cumulativeLength += np.linalg.norm(p2-p1)
            lengthAlongBoundary.append(cumulativeLength)
        lengthAlongBoundary = np.array(lengthAlongBoundary)

        if fixedBoundaryPoints != None:
            fixedRanksOnBoundary = [0]
            nFixedPointsOnBoundary = 1
            for fixedBoundaryPoint in fixedBoundaryPoints[1:]:
                inds, point = fixedBoundaryPoint[0], fixedBoundaryPoint[1]
                #indexInBars = np.where((inMesh.nodes[n:,:][:,inds] == point).all(axis=1))[0]
                indexInBars = (np.linalg.norm(np.subtract(inMesh.nodes[n:,:][:,inds], point), axis = 1)).argmin()

                fixedRanksOnBoundary.append(indexInBars)
                nFixedPointsOnBoundary += 1
            fixedRanksOnBoundary.append(-1)

            angles = []
            deltaAngle = 2*np.pi/nFixedPointsOnBoundary
            #print("deltaAngle =", deltaAngle)
            for k in range(nFixedPointsOnBoundary):

                deltaLengthAlongBoundary = lengthAlongBoundary[fixedRanksOnBoundary[k]:fixedRanksOnBoundary[k+1]]-lengthAlongBoundary[fixedRanksOnBoundary[k]]
                deltaUnitLengthAlongBoundary = deltaLengthAlongBoundary/(lengthAlongBoundary[fixedRanksOnBoundary[k+1]]-lengthAlongBoundary[fixedRanksOnBoundary[k]])
                res = (k+deltaUnitLengthAlongBoundary)*deltaAngle
                angles = np.hstack((angles, res))

            angles = np.hstack((angles, 2.*np.pi))

        else:
            angles = np.linspace(2*np.pi/nBoundary, 2*np.pi, nBoundary)
            #angles = 2*np.pi*lengthAlongBoundary/cumulativeLength

        if boundaryOrientation == "direct":
            for i, a in enumerate(angles):
                u[n+i,0] = np.cos(a)
                u[n+i,1] = np.sin(a)
        else:
            for i, a in enumerate(angles):
                u[n+i,0] = np.cos(a)
                u[n+i,1] = -np.sin(a)

    else:
        raise NotImplementedError("outShape"+str(outShape)+" not implemented")

    # Compute a node graphe for the mesh
    edges = set()
    elFilter = Filters.ElementFilter(mesh, dimensionality = 2, elementTypes = [ElementNames.Triangle_3])
    for name, data ,ids in elFilter:
        for face in ElementNames.faces[name]:
            for idd in ids:
                edge = np.sort(data.connectivity[idd][face[1]])
                edges.add((edge[0], edge[1]))

    G2 = InitializeGraphPointsFromMeshPoints(mesh)
    for edge in edges:
        G2.add_edge(edge[0], edge[1])

    # Compute the weights of each node of the mesh (number of edges linked to each node): the inverse of the degrees
    Ad = networkx.adjacency_matrix(G2)

    weights = np.zeros(N)
    for i in range(N):
        weights[i] = 1./np.sum(Ad[i,:])

    # Construct the sparse linear system to solve to find the position of the interior points in the parametrization


    A = sparse.eye(n).tolil()
    RHSmat = sparse.lil_matrix((n, N))
    for edge in edges:
        for edg in [(edge[0], edge[1]), (edge[1], edge[0])]:
            if edg[0] < n and edg[1] < n:
                A[edg[0], edg[1]] = -weights[edg[0]]
            elif edg[0] < n:
                RHSmat[edg[0], edg[1]] = weights[edg[0]]

    RHS = RHSmat.dot(u)
    A = A.tocsr()

    # update the position of the interior points
    res = sparse.linalg.spsolve(A, RHS)
    u[:n,:] = res


    if fixedInteriorPoints != None:
        mesh.nodes = u
        mesh.ConvertDataForNativeTreatment()

        displacement = None
        mask = None

        if "mean" in fixedInteriorPoints:

            meanPos = np.mean(u[fixedInteriorPoints["mean"][0],:], axis=0)

            if displacement == None:
                displacement = -np.tile(meanPos,(fixedInteriorPoints["mean"][1].shape[0],1))
            else:
                displacement = np.vstack((displacement, -np.tile(meanPos,(fixedInteriorPoints["mean"][1].shape[0],1))))

            if mask == None:
                mask = fixedInteriorPoints["mean"][0]
            else:
                mask = np.hstack((mask, fixedInteriorPoints["mean"][0]))

        if "value" in fixedInteriorPoints:

            if displacement == None:
                displacement = fixedInteriorPoints["value"][1] - u[fixedInteriorPoints["value"][0],:]
            else:
                displacement = np.vstack((displacement, fixedInteriorPoints["value"][1] - u[fixedInteriorPoints["value"][0],:]))

            if mask == None:
                mask = fixedInteriorPoints["value"][0]
            else:
                mask = np.hstack((mask, fixedInteriorPoints["value"][0]))


        if displacement is not None and mask is not None:
            displacement = np.vstack((displacement, np.zeros((N-n,2))))
            mask = np.hstack((mask, np.arange(n,N)))

            from BasicTools.Containers import UnstructuredMeshModificationTools as UMMT
            new_nodes = UMMT.Morphing(mesh, displacement, mask, radius= 1.)

            mesh.nodes = new_nodes

    else:
        mesh.nodes = u
        mesh.ConvertDataForNativeTreatment()

    infos = {}
    endgeLengths = []
    for edge in edges:
        endgeLengths.append(np.linalg.norm(mesh.nodes[edge[1],:]-mesh.nodes[edge[0],:]))

    infos = {"minEdge":np.min(endgeLengths), "maxEdge":np.max(endgeLengths)}

    return mesh, infos


def FloaterMesh3DParametrizationStarDomain(inMesh: UnstructuredMesh, boundaryTag: str, shape = None):
    """
    Naive implementation only for linear tetrahedral meshes and 3D star domains

    Parameters
    ----------
    inMesh : UnstructuredMesh
        input mesh of a 3D star domain
    boundaryTag : str
        elementtag for the boundary of the mesh
    shape : dict
        with key "type":str, "center":None or ndarray(3) of floats , "dimension":ndarray(3) of floats
        if shape == "None", the boundary of inMesh is mapped into a ball of radius 1, centered at the origin
        if "type"="cuboid", inMesh is mapped into a cuboid
        if "type"="ellipsoid", inMesh is mapped into a ellipsoid
        "center" defined the center of symmetry of the obtained parametrized mesh
        "dimension" defined the diameter of the obtained parametrized mesh along each direction

    Returns
    -------
    UnstructuredMesh
        parametrization of mesh
    ndarray(1) of ints
        renumbering of the nodes of the returned parametrized mesh, with respect to inMesh
    dict
        containing two keys: "minEdge" and "maxEdge", with values floats containing the minimal
        and maximal edged length of the parametrized mesh

    Notes
    -----
        mesh mush be a UnstructuredMesh of linear tetrahedrons, with an elementSet defining the boundary (the skin)
    """

    import copy
    mesh = copy.deepcopy(inMesh)

    # Retrieve the nodes of the surface boundary
    nf = Filters.NodeFilter(mesh, etag=boundaryTag)
    skinNodesIds = nf.GetIdsToTreat()


    N = mesh.GetNumberOfNodes()
    n = N - len(skinNodesIds)

    skinNodesIds = np.array(skinNodesIds)
    otherIds = np.delete(np.arange(N), nf.GetIdsToTreat())
    renumb2 = np.hstack((otherIds, skinNodesIds))
    invRenumb2 = np.argsort(renumb2)

    mesh.nodes = mesh.nodes[renumb2,:]
    for _, data in mesh.elements.items():
        data.connectivity = invRenumb2[data.connectivity]
    mesh.ConvertDataForNativeTreatment()


    # Set the boundary of the parametrization
    u = np.zeros((mesh.nodes.shape[0],3))

    mesh.ComputeBoundingBox()

    if shape == None:
        shape = {"type":"ellipsoid", "center":[0.,0.,0.], "dimension":[1.,1.,1.]}

    if "center" not in shape or shape["center"] == None:
        shape["center"] = (mesh.boundingMin + mesh.boundingMax)/2.

    if shape["type"] == "ellipsoid":

        for i in range(n, N):
            vect = mesh.nodes[i,:] - (mesh.boundingMin + mesh.boundingMax)/2.
            u[i,:] = np.multiply(shape["dimension"], vect)/np.linalg.norm(vect) + shape["center"]

    elif shape["type"] == "cuboid":

        delta = mesh.boundingMax - mesh.boundingMin

        ubound = mesh.nodes[n:,:] - ((mesh.boundingMin + mesh.boundingMax)/2.)[np.newaxis, :]
        coefs = np.divide(shape["dimension"], delta)

        nodes = np.multiply(coefs[np.newaxis, :], ubound) + shape["center"]

        u[n:,:] = nodes


    # Compute a node graphe for the mesh
    elFilter = Filters.ElementFilter(mesh, dimensionality = 3, elementTypes = [ElementNames.Tetrahedron_4])

    edges = set()
    for name, data ,ids in elFilter:
        for face in ElementNames.faces2[name]:
            for idd in ids:
                edge = np.sort(data.connectivity[idd][face[1]])
                edges.add((edge[0], edge[1]))


    G2 = InitializeGraphPointsFromMeshPoints(mesh)
    for edge in edges:
        G2.add_edge(edge[0], edge[1])

    # Compute the weights of each node of the mesh (number of edges linked to each node): the inverse of the degrees
    Ad = networkx.adjacency_matrix(G2)

    weights = np.zeros(N)
    for i in range(N):
        weights[i] = 1./np.sum(Ad[i,:])

    # Construct the sparse linear system to solve to find the position of the interior points in the parametrization
    A = sparse.eye(n).tolil()
    RHSmat = sparse.lil_matrix((n, N))

    for edge in edges:
        for edg in [(edge[0], edge[1]), (edge[1], edge[0])]:
            if edg[0] < n and edg[1] < n:
                A[edg[0], edg[1]] = -weights[edg[0]]
            elif edg[0] < n:
                RHSmat[edg[0], edg[1]] = weights[edg[0]]

    RHS = RHSmat.dot(u)
    A = A.tocsr()

    # update the position of the interior points
    res = sparse.linalg.spsolve(A, RHS)
    u[:n,:] = res

    mesh.nodes = u
    mesh.ConvertDataForNativeTreatment()

    infos = {}
    endgeLengths = []
    for edge in edges:
        endgeLengths.append(np.linalg.norm(mesh.nodes[edge[1],:]-mesh.nodes[edge[0],:]))

    infos = {"minEdge":np.min(endgeLengths), "maxEdge":np.max(endgeLengths)}


    return mesh, renumb2, infos



def FloaterMesh3DParametrization(inMesh: UnstructuredMesh, boundaryTag: str, inPlace: bool = True, fixedInteriorPoints = None, fixedBoundaryPoints = None):
    """
    STILL LARGELY EXPERIMENTAL

    Only for linear tetrahedral meshes

    Computes a 3D mesh parametrization by cutting a 3D mesh using a plane and applying the 2D mesh parametrization
    to the boundaries of the half-meshes on disks, projecting the meshes on two half-spheres, and then morphing the
    interior points.

    Parameters
    ----------
    inMesh : UnstructuredMesh
        Renumbered triangular mesh to parametrize
    boundaryTag : str
        element tag for the surface triangles of half of the boundary of inMesh (whichever one)
    inPlace : bool
        if "True", inMesh is modified
        if "False", inMesh is let unmodified, and a new mesh is produced
    boundaryOrientation : str
        if "direct, the boundary of the parametrisation is constructed in the direct trigonometric order
        if "indirect", the boundary of the parametrisation is constructed in the indirect trigonometric order
    fixedInteriorPoints : dict
        with two keys: "skinIds1", and "skinIds2", and dict as value:
            with one key, and corresponding value, a list: [ndarray(n), ndarray(n,2)],
            with n the number of interior points to be fixed; the first ndarray is the index of the considered
            interior point, the second ndarray is the corresponding prescribed positions
            if key is "mean", the interior points are displaced by the mean of the prescribed positions
            if key is "value", the interior points are displaced by the value of the prescribed positions
    fixedBoundaryPoints: list
        list of lists: [ndarray(2), ndarray(2)], helping definining a point in inMesh; the first ndarray is the component
        of a point on the boundary, and the second array is the value of corresponding component. Tested for triangular meshes
        in the 3D space.

    Returns
    -------
    UnstructuredMesh
        parametrization of the 3D mesh
    ndarray(1) of ints
        renumbering of the nodes of the returned parametrized mesh, with respect to inMesh
    dict
        containing 4 keys: "minEdge" and "maxEdge", with values floats containing the minimal
        and maximal edged length of the parametrized mesh; and "paramMesh1" and "paramMesh2",
        containing the 2D mesh parametrization (disks) of each half-boundary on inMesh

    Notes
    -----
        inMesh must have its skin partitionned in two, with boundaryTag being a element tag for
        surface triangles of one element of this partition
    """

    if fixedInteriorPoints == None:
        fixedInteriorPoints = {'skinIds1':None, 'skinIds2':None}


    if inPlace == True:
        meshTemp = inMesh
    else:
        import copy
        meshTemp = copy.deepcopy(inMesh)

    skinIds1 = meshTemp.elements['tri3'].GetTag(boundaryTag).GetIds()
    UMMT.ComputeSkin(meshTemp, md = 3, inPlace = True)
    skinIds2 = np.delete(meshTemp.elements['tri3'].GetTag('Skin').GetIds(), skinIds1)
    meshTemp.elements['tri3'].tags = T.Tags()
    meshTemp.elements['tri3'].GetTag('skinIds1').SetIds(skinIds1)
    meshTemp.elements['tri3'].GetTag('skinIds2').SetIds(skinIds2)

    # Treat 1st half-sphere
    ef = Filters.ElementFilter(meshTemp, dimensionality = 2, tag='skinIds1')
    skinHS1Mesh = UMIT.ExtractElementsByElementFilter(meshTemp, ef)
    UMMT.CleanLonelyNodes(skinHS1Mesh)

    renumSkinHS1Mesh, renumbHS1, nBoundaryHS1 = RenumberMeshForParametrization(skinHS1Mesh, inPlace = False, boundaryOrientation = "direct", fixedBoundaryPoints = fixedBoundaryPoints)
    parametrizedRenumSkinHS1Mesh, _ = FloaterMeshParametrization(renumSkinHS1Mesh, nBoundaryHS1, outShape = "circle", boundaryOrientation = "direct", fixedInteriorPoints = fixedInteriorPoints['skinIds1'], fixedBoundaryPoints = fixedBoundaryPoints)

    # Treat 2nd half-shpere
    ef2 = Filters.ElementFilter(meshTemp, dimensionality = 2, tag='skinIds2')
    skinHS2Mesh = UMIT.ExtractElementsByElementFilter(meshTemp, ef2)
    UMMT.CleanLonelyNodes(skinHS2Mesh)

    renumSkinHS2Mesh, renumbHS2, nBoundaryHS2 = RenumberMeshForParametrization(skinHS2Mesh, inPlace = False, boundaryOrientation = "direct", fixedBoundaryPoints = fixedBoundaryPoints)
    parametrizedRenumSkinHS2Mesh, _ = FloaterMeshParametrization(renumSkinHS2Mesh, nBoundaryHS2, outShape = "circle", boundaryOrientation = "indirect", fixedInteriorPoints = fixedInteriorPoints['skinIds2'], fixedBoundaryPoints = fixedBoundaryPoints)

    assert nBoundaryHS1 == nBoundaryHS2, "boundary between both boundarySet not compatible (different size)"

    np.testing.assert_allclose(skinHS1Mesh.originalIDNodes[renumbHS1[-nBoundaryHS1:]], skinHS2Mesh.originalIDNodes[renumbHS2[-nBoundaryHS2:]])

    rankInteriorHS1InMesh = skinHS1Mesh.originalIDNodes[renumbHS1[:-nBoundaryHS1]]
    rankInteriorHS2InMesh = skinHS2Mesh.originalIDNodes[renumbHS2[:-nBoundaryHS2]]
    rankLineBoundaryInMesh = skinHS1Mesh.originalIDNodes[renumbHS1[-nBoundaryHS1:]]
    np.testing.assert_allclose(meshTemp.nodes[rankInteriorHS1InMesh,:], renumSkinHS1Mesh.nodes[:-nBoundaryHS1,:])
    np.testing.assert_allclose(meshTemp.nodes[rankInteriorHS2InMesh,:], renumSkinHS2Mesh.nodes[:-nBoundaryHS2,:])
    np.testing.assert_allclose(meshTemp.nodes[rankLineBoundaryInMesh,:], renumSkinHS1Mesh.nodes[-nBoundaryHS1:,:])

    N = meshTemp.GetNumberOfNodes()
    rankSInMesh = np.union1d(rankInteriorHS1InMesh, rankInteriorHS2InMesh)
    rankSInMesh = np.union1d(rankSInMesh, rankLineBoundaryInMesh)
    otherIds = np.delete(np.arange(N),rankSInMesh)

    renumb = np.hstack((otherIds, rankInteriorHS1InMesh))
    renumb = np.hstack((renumb, rankInteriorHS2InMesh))
    renumb = np.hstack((renumb, rankLineBoundaryInMesh))

    invRenumb = np.argsort(renumb)
    meshTemp.nodes = meshTemp.nodes[renumb,:]
    for _, data in meshTemp.elements.items():
        data.connectivity = invRenumb[data.connectivity]
    meshTemp.ConvertDataForNativeTreatment()

    # Set the boundary of the parametrization
    n = N-len(rankInteriorHS1InMesh)-len(rankInteriorHS2InMesh)-nBoundaryHS1
    u = np.zeros((N,3))

    for i in range(n, n+len(rankInteriorHS1InMesh)):
        point = parametrizedRenumSkinHS1Mesh.nodes[i-n,:]
        l = np.linalg.norm(point)
        theta = np.pi*(1-l)/2
        phi = np.arccos(point[0]/l)
        u[i,0] = np.cos(theta)*np.cos(phi)
        u[i,1] = np.sign(point[1])*np.cos(theta)*np.sin(phi)
        u[i,2] = np.sin(theta)

    for i in range(n+len(rankInteriorHS1InMesh), n+len(rankInteriorHS1InMesh)+len(rankInteriorHS2InMesh)+nBoundaryHS1):
        point = parametrizedRenumSkinHS2Mesh.nodes[i-n-len(rankInteriorHS1InMesh),:]
        l = np.linalg.norm(point)
        if i<n+len(rankInteriorHS1InMesh)+len(rankInteriorHS2InMesh):
            theta = np.pi*(1-l)/2
        else:
            theta = 0.
        phi = np.arccos(point[0]/l)
        u[i,0] = np.cos(theta)*np.cos(phi)
        u[i,1] = -np.sign(point[1])*np.cos(theta)*np.sin(phi)
        u[i,2] = -np.sin(theta)


    # Compute a node graphe for the mesh
    edges = set()
    elFilter = Filters.ElementFilter(meshTemp, dimensionality = 3, elementTypes = [ElementNames.Tetrahedron_4])
    for name, data, ids in elFilter:
        for face in ElementNames.faces2[name]:
            for idd in ids:
                edge = np.sort(data.connectivity[idd][face[1]])
                edges.add((edge[0], edge[1]))

    G2 = InitializeGraphPointsFromMeshPoints(meshTemp)
    for edge in edges:
        G2.add_edge(edge[0], edge[1])

    # Compute the weights of each node of the mesh (number of edges linked to each node): the inverse of the degrees
    Ad = networkx.adjacency_matrix(G2)

    weights = np.zeros(N)
    for i in range(N):
        if 1./np.sum(Ad[i,:])>0:
            weights[i] = 1./np.sum(Ad[i,:])
        else:
            weights[i] = 0.

    # Construct the sparse linear system to solve to find the position of the interior points in the parametrization
    A = sparse.eye(n).tolil()
    RHSmat = sparse.lil_matrix((n, N))

    for edge in edges:
        for edg in [(edge[0], edge[1]), (edge[1], edge[0])]:
            if edg[0] < n and edg[1] < n:
                A[edg[0], edg[1]] = -weights[edg[0]]
            elif edg[0] < n:
                RHSmat[edg[0], edg[1]] = weights[edg[0]]

    RHS = RHSmat.dot(u)
    A = A.tocsr()

    # update the position of the interior points
    res = sparse.linalg.spsolve(A, RHS)
    u[:n,:] = res

    meshTemp.nodes = u
    meshTemp.ConvertDataForNativeTreatment()

    infos = {}
    endgeLengths = []
    for edge in edges:
        endgeLengths.append(np.linalg.norm(meshTemp.nodes[edge[1],:]-meshTemp.nodes[edge[0],:]))

    infos = {"minEdge":np.min(endgeLengths), "maxEdge":np.max(endgeLengths), "paramMesh1":parametrizedRenumSkinHS1Mesh, "paramMesh2":parametrizedRenumSkinHS2Mesh}


    return meshTemp, renumb, infos



def ComputeStretchMetric(originMesh, parametrizedMesh):
    """
    Computes the texture stretch metric between an inital triangle mesh
    and its parametrization, as defined in [1].

    Parameters
    ----------
    originMesh : UnstructuredMesh
        Renumbered triangular mesh to parametrize
    parametrizedMesh : int
        Parametrization of originMesh

    Returns
    -------
    np.ndarray, of size(nTriangles)
        texture stretch metric
    np.ndarray, of size(nTriangles)
        areas of the triangles of the mesh

    Notes
    -----
        mesh mush be a renumbered UnstructuredMesh of triangles (either in
        a 2D or 3D ambiant space), with a line boundary (no closed surface in 3D)

    References
    ----------
        [1] P.V. Sander, J. Snyder, S.J. Gortler, H. Hoppe. Texture mapping progressive meshes, 2001.
        URL: https://hhoppe.com/tmpm.pdf
    """
    assert(originMesh.GetNumberOfNodes() == parametrizedMesh.GetNumberOfNodes())
    connectivity = originMesh.elements[ElementNames.Triangle_3].connectivity

    np.testing.assert_almost_equal(connectivity, parametrizedMesh.elements[ElementNames.Triangle_3].connectivity)

    nTriangles = connectivity.shape[0]
    stretchMetric = np.zeros(nTriangles)
    area = np.zeros(nTriangles)

    q = np.zeros((3,2))

    for i in range(nTriangles):
        nodeIndices = connectivity[i]
        q[:,:] = originMesh.nodes[nodeIndices,:]
        s = parametrizedMesh.nodes[nodeIndices,0]
        t = parametrizedMesh.nodes[nodeIndices,1]

        Ax2 = (s[1]-s[0])*(t[2]-t[0])-(s[2]-s[0])*(t[1]-t[0])
        S_s = (q[0,:]*(t[1]-t[2]) + q[1,:]*(t[2]-t[0]) + q[2,:]*(t[0]-t[1]))/Ax2
        S_t = (q[0,:]*(s[1]-s[2]) + q[1,:]*(s[2]-s[0]) + q[2,:]*(s[0]-s[1]))/Ax2

        a = np.dot(S_s, S_s)
        b = np.dot(S_s, S_t)
        c = np.dot(S_t, S_t)

        rad = np.sqrt((a-c)**2+4*b**2)

        Gamma = np.sqrt((a+c+rad)/2)
        gamma = np.sqrt((a+c-rad)/2)

        stretchMetric[i] = np.sqrt((Gamma**2+gamma**2)/2)
        area[i] = 0.5*np.linalg.norm(np.cross(q[1,:]-q[0,:], q[2,:]-q[0,:]))

    return stretchMetric, area



def ComputeNodalAveragedStretchMetric(originMesh, parametrizedMesh):
    """
    Computes a nodal-averaged stretch metric between an inital triangle mesh
    and its parametrization, as defined in [1].

    Parameters
    ----------
    originMesh : UnstructuredMesh
        Renumbered triangular mesh to parametrize
    parametrizedMesh : int
        Parametrization of originMesh

    Returns
    -------
    np.ndarray, of size(nNodes)
        texture stretch metric

    Notes
    -----
        mesh mush be a renumbered UnstructuredMesh of triangles (either in
        a 2D or 3D ambiant space), with a line boundary (no closed surface in 3D)

    References
    ----------
        [1] S. Yoshizawa, A. Belyaev, H.-P. Seidel. A fast and simple stretch-minimizing mesh parameterization, 2004.
        URL, http://www2.riken.jp/brict/Yoshizawa/Papers/smi04ybs.pdf
    """
    stretchMetric, area = ComputeStretchMetric(originMesh, parametrizedMesh)

    NodeToElementConnectivity, _ = GetDualGraphNodeToElement(originMesh)

    nNodes = originMesh.GetNumberOfNodes()
    nodalStretchMetric = np.zeros(nNodes)

    for i in range(nNodes):
        tris = NodeToElementConnectivity[i][NodeToElementConnectivity[i]>0]
        stret2 = np.square(stretchMetric[tris])
        are = area[tris]
        nodalStretchMetric[i] = np.sqrt(np.dot(are, stret2) / np.sum(are))

    return nodalStretchMetric



def ComputeMeanEdgeAtPoints(mesh):
    """
    Computes a notion of node density at each vertex,
    defined by the mean length of the edged connected to this point in the mesh

    Parameters
    ----------
    originMesh : mesh
        input mesh whose node density is searched
    parametrizedMesh : int
        Parametrization of originMesh

    Returns
    -------
    np.ndarray, of size(nNodes)
         mean length of the edged connected to each point
    """
    N = mesh.GetNumberOfNodes()
    meanEdgeAtPoints = np.zeros(N)
    nodalGraph = GetDualGraph(mesh)[0]

    for i in range(N):
        linkedNodes = nodalGraph[i][nodalGraph[i]>0]
        d = len(linkedNodes)

        for j in linkedNodes:
            meanEdgeAtPoints[i] += np.linalg.norm(mesh.nodes[i] - mesh.nodes[j])/d

    return meanEdgeAtPoints


# INTEGRITY CHECKS # -----------------------------------------------------------------------------------


def CreateMeshForCheckIntegrity():
    from BasicTools.Containers import UnstructuredMeshCreationTools as UMCT
    from BasicTools.Containers import ConstantRectilinearMesh as CRM

    myMesh = CRM.ConstantRectilinearMesh(2)
    myMesh.SetDimensions([3, 3])
    myMesh.SetSpacing([1, 1])

    return UMCT.CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex = True)



def Create3DMeshForCheckIntegrity():
    from BasicTools.Containers import UnstructuredMeshCreationTools as UMCT
    from BasicTools.Containers import ConstantRectilinearMesh as CRM

    myMesh = CRM.ConstantRectilinearMesh(3)
    myMesh.SetDimensions([5, 5, 5])
    myMesh.SetSpacing([1, 1, 1])

    return UMCT.CreateMeshFromConstantRectilinearMesh(myMesh, ofSimplex = True)


def CheckIntegrity_PartitionMesh_Metis(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    dg, nused = GetElementGraph(res)

    # must convert dg in UnstructTuredMeshFormat for next line to run
    #mv = PartitionMesh(dg, 2, driver="Metis")
    print("test to complete")
    return "ok"


def CheckIntegrity_InitializeGraphPointsFromMeshPoints(GUI=False):

    mesh = CreateMeshForCheckIntegrity()
    G = InitializeGraphPointsFromMeshPoints(mesh)

    np.testing.assert_almost_equal(G.nodes, np.arange(9))
    #print("1 =", G.nodes)

    return "ok"


def CheckIntegrity_InitializeGraphPointsFromMeshElements(GUI=False):

    mesh = CreateMeshForCheckIntegrity()
    G = InitializeGraphPointsFromMeshElements(mesh)
    np.testing.assert_almost_equal(G.nodes, np.arange(8))
    #print("2 =", G.nodes)

    return "ok"



def CheckIntegrity_GetExtendedElementGraph(GUI=False):

    mesh = CreateMeshForCheckIntegrity()
    G = GetExtendedElementGraph(mesh)

    np.testing.assert_almost_equal(G.nodes, np.arange(8))
    #print("3 =", G.nodes)

    return "ok"


def CheckIntegrity_ComputeMeshLaplacianEigenmaps(GUI=False):

    mesh = CreateMeshForCheckIntegrity()
    w, v = ComputeMeshLaplacianEigenmaps(mesh, nEigenmaps = 3)

    refW = np.array([-1.2490009e-15,  1.0000000e+00,  2.1559265e+00])
    np.testing.assert_almost_equal(w, refW)

    refV = np.array([[3.33333333e-01,  2.06105878e-16,  3.59572691e-01],
    [ 3.33333333e-01, -2.88675135e-01,  4.50452206e-02],
    [ 3.33333333e-01, -5.77350269e-01, -5.77775050e-01]])
    np.testing.assert_almost_equal(v[:3,:], refV)


    return "ok"



def CheckIntegrity_FloaterMeshParametrization(GUI=False):

    mesh = CreateMeshForCheckIntegrity()

    meshRenumb, renumb, nBoundary = RenumberMeshForParametrization(mesh, inPlace = False)
    meshParam, infos = FloaterMeshParametrization(meshRenumb, nBoundary)

    print("infos FloaterMeshParametrization:", infos)

    np.testing.assert_almost_equal(infos['minEdge'], 0.7653668647301795)
    np.testing.assert_almost_equal(infos['maxEdge'], 1.4142135623730951)

    refNodes = np.array([[ 0.00000000e+00, -1.30659844e-17],
    [ 7.07106781e-01,  7.07106781e-01],
    [ 6.12323400e-17,  1.00000000e+00],
    [-7.07106781e-01,  7.07106781e-01],
    [-1.00000000e+00,  1.22464680e-16],
    [-7.07106781e-01, -7.07106781e-01],
    [-1.83697020e-16, -1.00000000e+00],
    [ 7.07106781e-01, -7.07106781e-01],
    [ 1.00000000e+00, -2.44929360e-16]])

    np.testing.assert_almost_equal(meshParam.nodes, refNodes)

    return "ok"


def CheckIntegrity_FloaterMesh3DParametrizationStarDomain(GUI=False):

    mesh = Create3DMeshForCheckIntegrity()

    from BasicTools.Containers import UnstructuredMeshCreationTools as UMCT
    UMCT.ComputeSkin(mesh, md = 3, inPlace = True)

    meshParam, renumb, infos = FloaterMesh3DParametrizationStarDomain(mesh, boundaryTag = 'Skin')

    print("infos FloaterMesh3DParametrizationStarDomain ellipsoid:", infos)


    np.testing.assert_almost_equal(infos['minEdge'], 0.27477100047357156)
    np.testing.assert_almost_equal(infos['maxEdge'], 0.8023930268619969)

    refNodes = np.array([[-0.36304405, -0.36304405, -0.36304405],
                         [-0.39871955, -0.39871955,  0.01773299],
                         [-0.40280694, -0.40280694,  0.43490912],
                         [-0.39871955,  0.01773299, -0.39871955]])

    np.testing.assert_almost_equal(meshParam.nodes[:4,:], refNodes)


    meshParam, renumb, infos = FloaterMesh3DParametrizationStarDomain(mesh, boundaryTag = 'Skin', shape = {"type":"cuboid", "center":[1.,0.5,0.1], "dimension":[0.1,0.2,0.3]})

    print("infos FloaterMesh3DParametrizationStarDomain cuboid   :", infos)

    np.testing.assert_almost_equal(infos['minEdge'], 0.02499999999999969)
    np.testing.assert_almost_equal(infos['maxEdge'], 0.09354143466934872)

    refNodes = np.array([[0.975, 0.45,  0.025],
                         [0.975, 0.45,  0.1  ],
                         [0.975, 0.45,  0.175],
                         [0.975, 0.5 ,  0.025]])

    np.testing.assert_almost_equal(meshParam.nodes[:4,:], refNodes)

    return "ok"


def CheckIntegrity_FloaterMesh3DParametrization(GUI=False):

    mesh = Create3DMeshForCheckIntegrity()

    from BasicTools.Containers import UnstructuredMeshCreationTools as UMCT
    UMCT.ComputeSkin(mesh, md = 3, inPlace = True)


    mesh.ComputeBoundingBox()

    eF = Filters.ElementFilter(mesh, zone = lambda p: (p[:,1]-0.5*(mesh.boundingMin[1]+mesh.boundingMax[1])), zoneTreatment = "allnodes")
    skinIds1 = eF.GetIdsToTreat(mesh.elements['tri3'])
    mesh.elements['tri3'].GetTag('skinIds1').SetIds(skinIds1)

    meshParam, renumb, infos = FloaterMesh3DParametrization(mesh, boundaryTag = 'skinIds1', inPlace = False)

    print("infos FloaterMesh3DParametrization:", infos)

    np.testing.assert_almost_equal(infos['minEdge'], 0.15893069864588907)
    np.testing.assert_almost_equal(infos['maxEdge'], 0.7653668647301797)

    refNodes = np.array([[ 0.33118664,  0.137182,    0.45861232],
                         [-0.04485619,  0.34618882,  0.45766116],
                         [-0.42191731,  0.48060774,  0.41107971],
                         [ 0.51377967,  0.21281451,  0.01164092]])

    np.testing.assert_almost_equal(meshParam.nodes[:4,:], refNodes)

    positions1 = np.array(
    [[0.2 , 0.],
     [0.  , 0.2],
     [-0.2, 0.],
     [0.  , -0.2]])

    positions2 = np.array(
    [[0.45,0.05]])

    fixedInteriorPoints = {}
    fixedInteriorPoints['skinIds1'] = {}
    fixedInteriorPoints['skinIds1']["mean"] = [np.array([10, 0, 2, 11]), positions1]
    fixedInteriorPoints['skinIds2'] = {}
    fixedInteriorPoints['skinIds2']["value"] = [np.array([6]), positions2]

    meshParam, renumb, infos = FloaterMesh3DParametrization(mesh, boundaryTag = 'skinIds1', inPlace = False, fixedInteriorPoints = fixedInteriorPoints)

    print("infos FloaterMesh3DParametrization :", infos)

    np.testing.assert_almost_equal(infos['minEdge'], 0.06353546953412406)
    np.testing.assert_almost_equal(infos['maxEdge'], 0.8096752560898939)

    refNodes = np.array([[ 0.21487737,  0.04148949,  0.50421315],
                         [-0.16425121,  0.25739021,  0.48133962],
                         [-0.51055236,  0.41595872,  0.41360663],
                         [ 0.48672766,  0.14950934,  0.03759597]])

    np.testing.assert_almost_equal(meshParam.nodes[:4,:], refNodes)

    return "ok"


def CheckIntegrity_ComputeNodalAveragedStretchMetric(GUI=False):

    mesh = CreateMeshForCheckIntegrity()

    meshRenumb, renumb, nBoundary = RenumberMeshForParametrization(mesh, inPlace = False)
    meshParam, infos = FloaterMeshParametrization(meshRenumb, nBoundary)

    res = ComputeNodalAveragedStretchMetric(meshRenumb, meshParam)
    ref = np.array([1.16252822, 1.25928013, 1.48563346, 1.84775907, 1.41421356, 1.25928013, 1.41421356, 1.84775907, 1.41421356])

    np.testing.assert_almost_equal(res, ref)

    return "ok"



def CheckIntegrity_ComputeMeanEdgeAtPoints(GUI=False):

    mesh = CreateMeshForCheckIntegrity()

    meshRenumb, renumb, nBoundary = RenumberMeshForParametrization(mesh, inPlace = False)
    meshParam, infos = FloaterMeshParametrization(meshRenumb, nBoundary)

    res = ComputeMeanEdgeAtPoints(meshParam)
    ref = np.array([1., 0.76536686, 0.9816491, 0.76536686, 0.9816491, 0.76536686, 0.9816491, 0.76536686, 0.9816491])

    np.testing.assert_almost_equal(res, ref)

    return "ok"





def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_PartitionMesh_Metis,
    CheckIntegrity_InitializeGraphPointsFromMeshPoints,
    CheckIntegrity_InitializeGraphPointsFromMeshElements,
    CheckIntegrity_GetExtendedElementGraph,
    CheckIntegrity_ComputeMeshLaplacianEigenmaps,
    CheckIntegrity_FloaterMeshParametrization,
    CheckIntegrity_FloaterMesh3DParametrizationStarDomain,
    CheckIntegrity_FloaterMesh3DParametrization,
    CheckIntegrity_ComputeNodalAveragedStretchMetric,
    CheckIntegrity_ComputeMeanEdgeAtPoints
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover