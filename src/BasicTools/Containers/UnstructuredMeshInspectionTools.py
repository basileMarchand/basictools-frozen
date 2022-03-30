# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

 #   if elementFilter is not None and zone is not None:
 #       raise(Exception("Need only one filter or zone"))
 #   if elementFilter is None and zone is None:
 #       raise(Exception("Need a filter or a zone"))
 #
 #   if elementFilter is not None:
 #       zones = filter.zones


from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType
import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes
from BasicTools.NumpyDefs import PBasicIndexType

def GetElementsFractionInside(field, points, name, elements, ids):

    from scipy.spatial import Delaunay

    def triangle_lob_fraction(index, phis):
        div = phis-phis[index]
        div[index] = phis[index]
        val = phis[index]/(div)
        return np.prod(val)

    def extract_signs(phis):
        signs = np.sign(phis)
        dominant_sign = np.sign(np.sum(signs))
        if dominant_sign == -1:
            signs[signs != 1] = -1
        else:
            signs[signs != -1] = 1
        return signs

    def split_signs(signs):
        minuses = np.nonzero(signs == -1)[0]
        pluses = np.nonzero(signs == 1)[0]
        return minuses, pluses


    # helper fucntion for tets

    def tetrahedron_volumic_fraction(phis):
        assert phis.size == 4
        assert np.count_nonzero(phis) > 0
        minuses, pluses = sort_by_sign(phis)
        if minuses.size == 0:
            return 0.0
        elif minuses.size == 1:
            return tetrahedron_volumic_fraction_dominated(minuses, pluses)
        elif minuses.size == 2:
            return tetrahedron_volumic_fraction_nondominated(minuses, pluses)
        elif minuses.size == 3:
            return 1.0 - tetrahedron_volumic_fraction_dominated(pluses, minuses)
        else:
            return 1.0

    def tetrahedron_volumic_fraction_dominated(phi_in, phis_out):
        return np.prod(phi_in / (phi_in - phis_out))

    def tetrahedron_volumic_fraction_nondominated(phis_in, phis_out):
        phis_in_r = np.reshape(phis_in, (2, 1))
        phis_out_r = np.reshape(phis_out, (1, 2))
        ratios = phis_in_r / (phis_in_r - phis_out_r)
        result = ratios[0, 0] * ratios[0, 1] * (1.0 - ratios[1, 0]) + \
                ratios[1, 0] * ratios[1, 1] * (1.0 - ratios[0, 1]) + \
                ratios[1, 0] * ratios[0, 1]
        return result

    def sort_by_sign(phis):
        signs = np.sign(phis)
        dominant_sign = np.sign(np.sum(signs))
        if dominant_sign == -1:
            minuses = phis[signs != 1]
            pluses = phis[signs == 1]
        else:
            minuses = phis[signs == -1]
            pluses = phis[signs != -1]
        return minuses, pluses

    #-----------------------------------------------------------------------------
    res = np.zeros(len(ids)) -1

    for cpt,id in enumerate(ids):
        coon = elements.connectivity[id,:]
        vals = field[coon]
        if np.all(vals<0)  :
            res[cpt] = 1
        elif np.all(vals>0) :
            res[cpt] = 0
        else:
            elemPoints = points[coon,:]
            simplex = Delaunay(points).simplices
            if ElementNames.dimension[name] == 1 and len(vals) == 2:
                #bar2
                if vals[0] < 0 :
                    res[cpt] = abs(vals[0])/np.sup(abs(vals))
                else:
                    res[cpt] = 1-abs(vals[0])/np.sup(abs(vals))
            elif ElementNames.dimension[name] == 2 and len(vals) == 3:
               # triangles
                signs = extract_signs(vals)
                minuses, pluses = split_signs(signs)

                if minuses.size == 2:
                    res[cpt] =  1.0 - triangle_lob_fraction(pluses[0], vals)

                elif minuses.size == 1:
                    res[cpt] =  triangle_lob_fraction(minuses[0], vals)
            elif ElementNames.dimension[name] == 3 and len(vals) == 4:
                # tets
                res[cpt] = tetrahedron_volumic_fraction(vals)

    return res

def VolumeOfTetrahedrons(inmesh):

    elems =inmesh.GetElementsOfType(ElementNames.Tetrahedron_4)
    conn = elems.connectivity
    a = inmesh.nodes[conn[:,0],:]
    b = inmesh.nodes[conn[:,1],:]
    c = inmesh.nodes[conn[:,2],:]
    d = inmesh.nodes[conn[:,3],:]
    e = np.cross(b-d,c-d)
    f = (a-d)
    res = np.empty(elems.GetNumberOfElements(),dtype=np.float)
    for n in range(elems.GetNumberOfElements()):
        res[n] = np.abs( np.dot(f[n,:],e[n,:])  )

    return res*(1./6.)

def VolumeOfHexaedrons(inmesh):

    elems =inmesh.GetElementsOfType(ElementNames.Hexaedron_8)
    conn = elems.connectivity

    def VolumeInternal(a,b,c,d):
        e = np.cross(b-d,c-d)
        f = (a-d)
        res = np.empty(elems.GetNumberOfElements(),dtype=np.float)
        for n in range(elems.GetNumberOfElements()):
            res[n] = np.abs( np.dot(f[n,:],e[n,:])  )
        return res*(1./6.)

    res = np.zeros(elems.GetNumberOfElements(),dtype=np.float)

    p0 = inmesh.nodes[conn[:,0],:]
    p1 = inmesh.nodes[conn[:,1],:]
    p2 = inmesh.nodes[conn[:,2],:]
    p3 = inmesh.nodes[conn[:,3],:]
    p4 = inmesh.nodes[conn[:,4],:]
    p5 = inmesh.nodes[conn[:,5],:]
    p6 = inmesh.nodes[conn[:,6],:]
    p7 = inmesh.nodes[conn[:,7],:]

    res += VolumeInternal(p0,p6,p5,p1)
    res += VolumeInternal(p0,p6,p1,p2)
    res += VolumeInternal(p0,p6,p2,p3)
    res += VolumeInternal(p0,p6,p3,p7)
    res += VolumeInternal(p0,p6,p7,p4)
    res += VolumeInternal(p0,p6,p4,p5)

    return res

def GetVolume(inmesh) :

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo, ConstantSpaceGlobal
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.FE.Integration import IntegrateGeneral

    numbering = ComputeDofNumbering(inmesh,LagrangeSpaceGeo,fromConnectivity=True)

    wform = GetField("F",1).T*GetTestField("T",1)

    F = FEField("F",inmesh,LagrangeSpaceGeo,numbering)
    F.Allocate(1.)
    gnumbering = ComputeDofNumbering(inmesh,ConstantSpaceGlobal)
    unkownFields = [ FEField("T",mesh=inmesh,space=ConstantSpaceGlobal,numbering=gnumbering) ]
    _,f  = IntegrateGeneral( mesh=inmesh, wform=wform, constants={}, fields=[F], unkownFields=unkownFields)
    return f[0]

def GetDualGraphNodeToElement(inmesh, maxNumConnections=200):
    # generation of the dual graph
    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=PBasicIndexType )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType )

    cpt =0

    for name,elems in inmesh.elements.items():
        for i in range(elems.GetNumberOfElements()):
            coon = elems.connectivity[i,:]
            for j in coon:
                dualGraph[j,usedPoints[j]] =  cpt
                usedPoints[j] += 1
            cpt += 1

    #we crop the output data
    maxsize = np.max(np.sum(dualGraph>=0,axis=1))
    dualGraph = dualGraph[:,0:maxsize]
    return dualGraph,usedPoints

def GetDualGraph(inmesh, maxNumConnections=200):

    # generation of the dual graph
    dualGraph = np.zeros((inmesh.GetNumberOfNodes(),maxNumConnections), dtype=PBasicIndexType )-1
    usedPoints = np.zeros(inmesh.GetNumberOfNodes(), dtype=PBasicIndexType )

    for name,elems in inmesh.elements.items():
        size = elems.GetNumberOfNodesPerElement()
        for i in range(elems.GetNumberOfElements()):
            coon = elems.connectivity[i,:]
            for j in range(size):
                myIndex = coon[j]
                for k in range(size):
                    if k == j:
                        continue
                    dualGraph[myIndex,usedPoints[myIndex]] =  coon[k]
                    usedPoints[myIndex] += 1
                    # we reached the maximun number of connection
                    # normally we have some data duplicated
                    # we try to shrink the vector
                    if usedPoints[myIndex] == maxNumConnections :# pragma: no cover
                        # normaly we shouldn't pas here
                        c = np.unique(dualGraph[myIndex,:])
                        dualGraph[myIndex,0:len(c)] = c
                        usedPoints[myIndex] = len(c)

    maxsize = 0
    # we finish now we try to compact the structure
    for i in range(inmesh.GetNumberOfNodes()):
        c = np.unique(dualGraph[i,0:usedPoints[i]])
        dualGraph[i,0:len(c)] = c
        dualGraph[i,len(c):] = -1
        usedPoints[i] = len(c)
        maxsize = max(len(c),maxsize)

    #we crop the output data
    dualGraph = dualGraph[:,0:maxsize]
    return dualGraph,usedPoints

def ExtractElementsByImplicitZone(inmesh,op,allNodes=True,cellCenter=False):
    inmesh.ComputeGlobalOffset()

    mask = op(inmesh.nodes) <= 0.

    outmesh = type(inmesh)()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = np.copy(inmesh.nodes)
    outmesh.originalIDNodes = np.arange(inmesh.GetNumberOfNodes())

    # keep only the faces with one or zero volumes attached
    for name,data in inmesh.elements.items():
        if allNodes:
            ElemMask = np.sum(mask[data.connectivity],axis=1) == data.GetNumberOfNodesPerElement()
        else:
            ElemMask = np.sum(mask[data.connectivity],axis=1) >= 1
        outmesh.elements[name] = ExtractElementsByMask(data,ElemMask)
    CleanLonelyNodes(outmesh)
    outmesh.PrepareForOutput()
    return outmesh

def ExtractElementsByElementFilter(inmesh,ff):

    inmesh.ComputeGlobalOffset()
    outmesh = type(inmesh)()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = np.copy(inmesh.nodes)
    outmesh.originalIDNodes = np.arange(inmesh.GetNumberOfNodes())
    outmesh.nodesTags = inmesh.nodesTags.Copy()

    ff.mesh = inmesh
    for name,data,ids in ff:
        outmesh.elements[name] = ExtractElementsByMask(data,ids)

    outmesh.PrepareForOutput()
    return outmesh

def ExtractElementsByMask(inelems, _mask):


    outelems = type(inelems)(inelems.elementType)

    newIndex = np.empty(inelems.GetNumberOfElements(),dtype=PBasicIndexType)

    _mask = np.asarray(_mask)
    if _mask.dtype == bool:
        nbels =0
        for i in range(inelems.GetNumberOfElements()):
           newIndex[i] = nbels
           nbels += 1 if _mask[i] else 0
        mask = _mask
    else:
        nbels = len(_mask)
        mask = np.zeros(inelems.GetNumberOfElements(),dtype=bool)
        cpt =0
        for index in _mask:
           newIndex[index ] = cpt
           mask[index] = True
           cpt += 1

    outelems.Allocate(nbels)
    outelems.connectivity = inelems.connectivity[mask,:]
    outelems.originalIds = np.where(mask)[0]

    for tag in inelems.tags  :
       temp = np.extract(mask[tag.GetIds()],tag.GetIds())
       newid = newIndex[temp]
       outelems.tags.CreateTag(tag.name).SetIds(newid)

    return outelems

def ExtractElementByDimensionalityNoCopy(inmesh,dimensionalityFilter):
    outmesh = type(inmesh)()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = inmesh.nodes
    outmesh.originalIDNodes = inmesh.originalIDNodes
    outmesh.nodesTags = inmesh.nodesTags

    if dimensionalityFilter >0 :
        for name,elems in inmesh.elements.items():
            if ElementNames.dimension[name] != dimensionalityFilter:
                continue
            else:
                outmesh.elements[name] = elems
    else:
        for name,elems in inmesh.elements.items():
            if ElementNames.dimension[name] == -dimensionalityFilter:
                continue
            else:
                outmesh.elements[name] = elems

    outmesh.PrepareForOutput()
    return outmesh

def ExtractElementByTags(inmesh,tagsToKeep, allNodes=False,dimensionalityFilter= None, cleanLonelyNodes=True):

    outmesh = UnstructuredMesh()
    outmesh.CopyProperties(inmesh)

    outmesh.nodes = np.copy(inmesh.nodes)
    outmesh.originalIDNodes = np.copy(inmesh.originalIDNodes)

    import copy
    for tag in inmesh.nodesTags:
        outmesh.nodesTags.AddTag(copy.deepcopy(tag) )


    nodalMask = np.zeros(inmesh.GetNumberOfNodes(),dtype = bool)
    for name,elems in inmesh.elements.items():

       #if dimensionalityFilter is not None:
       #    if dimensionalityFilter !=  ElementNames.dimension[name]:
       #        continue

       if (np.any([x in elems.tags for x in tagsToKeep] ) == False) and (dimensionalityFilter is None) :
           if np.any([x in inmesh.nodesTags for x in tagsToKeep]) == False:
               continue# pragma: no cover


       toKeep = np.zeros(elems.GetNumberOfElements(), dtype=bool)
       # check elements tags
       for tagToKeep in tagsToKeep:
           if tagToKeep in elems.tags:
               toKeep[elems.tags[tagToKeep].GetIds()] = True

       # check for nodes tags
       for tagToKeep in tagsToKeep:
           if tagToKeep in inmesh.nodesTags:
             nodalMask.fill(False)
             tag = inmesh.GetNodalTag(tagToKeep)
             nodalMask[tag.GetIds()] = True
             elemMask = np.sum(nodalMask[elems.connectivity],axis=1)
             if allNodes :
                 toKeep[elemMask == elems.GetNumberOfNodesPerElement()] = True
             else:
                 toKeep[elemMask > 0] = True

       # if dimensionality is ok and no tag to keep, we keep all elements
       if dimensionalityFilter is not None:
           if dimensionalityFilter ==  ElementNames.dimension[name]:
               toKeep[:] = True
           #if len(tagsToKeep)  == 0:
           #    toKeep[:] = True

       newIndex = np.empty(elems.GetNumberOfElements(), dtype=PBasicIndexType )
       cpt =0
       for i in range(elems.GetNumberOfElements()):
           newIndex[i] = cpt
           cpt += 1 if toKeep[i] else 0





       outelem = outmesh.GetElementsOfType(name)
       nbTokeep = np.sum(toKeep)
       outelem.Allocate(nbTokeep)
       outelem.connectivity = elems.connectivity[toKeep,:]
       outelem.originalIds = np.where(toKeep)[0]

       for tag in elems.tags  :
           temp = np.extract(toKeep[tag.GetIds()],tag.GetIds())
           newid = newIndex[temp]
           outelem.tags.CreateTag(tag.name,errorIfAlreadyCreated=False).SetIds(newid)

    if cleanLonelyNodes:
        CleanLonelyNodes(outmesh)
    outmesh.PrepareForOutput()
    return outmesh

def EnsureUniquenessElements(mesh):
    cpt = 0
    for name,data in mesh.elements.items():
        dd = dict()
        for el in range(data.GetNumberOfElements()):
            n = tuple(np.sort(data.connectivity[el,:]))
            if len(n) != len(np.unique(n)):
                raise(Exception("("+str(name)+") element " + str(el) +" (global["+str(el+cpt)+"])" + " use a point more than once ("+str(data.connectivity[el,:])+")"  ))
            if n in dd.keys():
                raise(Exception("("+str(name)+") element " + str(el) +" (global["+str(el+cpt)+"])" + " is a duplication of element " + str(dd[n]) +" (global["+str(dd[n]+cpt)+"])" ))
            dd[n] = el
        cpt = data.GetNumberOfElements()

def MeshQualityAspectRatioBeta(mesh):
    #https://cubit.sandia.gov/public/15.2/help_manual/WebHelp/mesh_generation/mesh_quality_assessment/tetrahedral_metrics.htm
    for name,data in mesh.elements.items():
        if ElementNames.dimension[name] == 0:
            pass
        elif ElementNames.dimension[name] == 1:
            pass
        elif ElementNames.dimension[name] == 2:
            pass
        elif ElementNames.dimension[name] == 3:
            if name == ElementNames.Tetrahedron_4:
               mmax = 0
               for el in range(data.GetNumberOfElements()):
                   n = data.connectivity[el,:]
                   nodes = mesh.nodes[n,:]
                   p0 = nodes[0,:]
                   p1 = nodes[1,:]
                   p2 = nodes[2,:]
                   p3 = nodes[3,:]
                   #print("---")

                   a = np.linalg.norm(p1-p0)
                   b = np.linalg.norm(p2-p0)
                   #print(a)
                   #print(b)
                   normal = np.cross(p1-p0,p2-p0)
                   #https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
                   base_area = 0.5*a*b*np.sqrt(1-(np.dot(p1-p0,p2-p0)/(a*b))**2)

                   normal /= np.linalg.norm(normal)

                   # distance to the oposite point
                   d = np.dot(normal,p3-p0)

                   volume = base_area*d
                   #print(base_area)
                   #print(d)
                   inscribed_sphere_radius = d/3
                   #https://math.stackexchange.com/questions/2820212/circumradius-of-a-tetrahedron

                   c = np.linalg.norm(p3-p0)

                   A = np.linalg.norm(p3-p2)
                   B = np.linalg.norm(p1-p3)
                   C = np.linalg.norm(p2-p1)

                   #print("-")
                   #print((a*A+b*B+c*C)*
                   #                                      (a*A+b*B-c*C)*
                   #                                      (a*A-b*B+c*C)*
                   #                                      (-a*A+b*B+c*C))
                   circumradius_sphere_radius = np.sqrt((a*A+b*B+c*C)*
                                                         (a*A+b*B-c*C)*
                                                         (a*A-b*B+c*C)*
                                                         (-a*A+b*B+c*C))/(24*volume)


                   AspectRatioBeta = circumradius_sphere_radius/(3.0 * inscribed_sphere_radius)
                   #print(volume)
                   #print(circumradius_sphere_radius)
                   #print(inscribed_sphere_radius)
                   #print(circumradius_sphere_radius)
                   #print(AspectRatioBeta)
                   mmax = max([mmax,AspectRatioBeta])
                   if AspectRatioBeta > 1000:
                       raise(Exception("Element " +str(el) + " has quality of " +str(AspectRatioBeta)) )
               #print("mmax")
               #print(mmax)

        else:
            raise

def ComputeMeshMinMaxLengthScale(mesh):
    (resMin,resMax) = (None,None)
    for name,data in mesh.elements.items():
        if data.GetNumberOfNodesPerElement() < 2: continue
        if data.GetNumberOfElements() == 0: continue
        posx = mesh.nodes[data.connectivity,0]
        posy = mesh.nodes[data.connectivity,1]
        meanx = np.sum(posx,axis=1)/data.GetNumberOfNodesPerElement()
        meany = np.sum(posy,axis=1)/data.GetNumberOfNodesPerElement()
        meanx.shape = (len(meanx),1)
        meany.shape = (len(meanx),1)


        distToBaricenter2 = (posx-meanx)**2 + (posy-meany)**2

        if mesh.nodes.shape[1] == 3:
            posz = mesh.nodes[data.connectivity,2]
            meanz = np.sum(posz,axis=1)/data.GetNumberOfNodesPerElement()
            meanz.shape = (len(meanx),1)
            distToBaricenter2 += (posz-meanz)**2

        distToBaricenter = np.sqrt(distToBaricenter2)
        mmin = np.min(distToBaricenter)
        mmax = np.max(distToBaricenter)
        if resMin is None:
            resMin = mmin
        else:
            resMin = min(resMin,mmin)

        if resMax is None:
            resMax = mmax
        else:
            resMax = min(resMax,mmax)
    return (2*resMin,2*resMax)

def PrintMeshInformation(mesh):
    from BasicTools.Helpers.TextFormatHelper import TFormat as TF

    def L25(text):
        return TF.Left(text ,fill=" ",width=25)
    def L15(text):
        return TF.Left(text ,fill=" ",width=15)

    mesh.VerifyIntegrity()

    print(TF.Center("Nodes information"))
    print(L25("nodes.shape:") , str(mesh.nodes.shape) )
    print(L25("nodes.dtype:") , str(mesh.nodes.dtype) )
    mesh.ComputeBoundingBox()
    print(L25("boundingMin:") + str(mesh.boundingMin) )
    print(L25("boundingMax:") + str(mesh.boundingMax) )
    print(L25("originalIDNodes.shape:") , str(mesh.originalIDNodes.shape) )
    print(L25("min(originalIDNodes):") , str(min(mesh.originalIDNodes) )  )
    print(L25("max(originalIDNodes):") , str(max(mesh.originalIDNodes) )  )
    print("---- nodes tags ---- " )
    print(mesh.nodesTags)
    for tag in mesh.nodesTags:
        res = L15(tag.name) + L15(" size:"+ str(len(tag)))
        if len(tag):
            res +=L15(" min:"+str(min(tag.GetIds()) ) )+ "  max:"+str(max(tag.GetIds()) )
        print(res)


    print(TF.Center("Elements information"))

    for name,data in mesh.elements.items():
        res = L25(str(type(data)).split("'")[1].split(".")[-1]+ " " + name+" ")
        res += L15("size:" + str(data.GetNumberOfElements()))
        res += L25("min(connectivity):" + str(min(data.connectivity.ravel())))
        res += L25("max(connectivity):" +  str(max(data.connectivity.ravel())))
        print(res)
        for tag in data.tags:
            res = L15(tag.name) + L15(" size:"+ str(len(tag)))
            if len(tag):
                res +=L15(" min:"+str(min(tag.GetIds()) ) )+ "  max:"+str(max(tag.GetIds()) )
            print("  Tag: " +res)

    print(TF.Center(""))
#------------------------- CheckIntegrity ------------------------
def CheckIntegrity_ExtractElementByTags(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.AddElementToTagUsingOriginalId(0,"first")
    res.AddElementToTagUsingOriginalId(1,"second")
    res.AddElementToTagUsingOriginalId(0,"all")
    res.AddElementToTagUsingOriginalId(1,"all")
    res.GetNodalTag("Point3").AddToTag(3)

    if ExtractElementByTags(res,["first"] ).GetNumberOfElements() != 1:
        raise# pragma: no cover
    if ExtractElementByTags(res,["Point3"] ).GetNumberOfElements() != 1:
        raise# pragma: no cover
    if ExtractElementByTags(res,["Point3"],allNodes=True ).GetNumberOfElements() != 0:
        raise# pragma: no cover

    if ExtractElementByTags(res,[],dimensionalityFilter=2 ).GetNumberOfElements() != 2:
        raise# pragma: no cover
    return "ok"

def CheckIntegrity_GetVolume(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf,CreateMeshFromConstantRectilinearMesh

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    vol = GetVolume(mesh)
    if vol != (1./6.):
        raise Exception('Error en the calculation of the volume')# pragma: no cover

    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh()
    myMesh.SetDimensions([3,3,3])
    myMesh.SetSpacing([0.5, 0.5,0.5])
    print(myMesh)
    vol  = GetVolume(CreateMeshFromConstantRectilinearMesh(myMesh))
    if abs(vol-1.) > 1e-8 :
        raise Exception('Error en the calculation of the volume')# pragma: no cover

    return "ok"


def CheckIntegrity_ExtractElementsByImplicitZone(GUI=False):
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshFromConstantRectilinearMesh
    myMesh = ConstantRectilinearMesh(dim=3)
    myMesh.SetDimensions([20,30,40])
    myMesh.SetOrigin([-1.0,-1.0,-1.0])
    myMesh.SetSpacing([2., 2.,2]/myMesh.GetDimensions())
    print(myMesh)
    res2 = CreateMeshFromConstantRectilinearMesh(myMesh,ofTetras=False)
    print(res2)

    class OPSphere(object):
        def __init__(self):
            self.center = np.array([0.0,0.0,0.0],dtype=float)
            self.radius = 0.5
        def __call__(self,pos):
            res = np.sqrt(np.sum((pos-self.center)**2,axis=1))-self.radius
            return res

    myOp = OPSphere()

    res = ExtractElementsByImplicitZone(res2,myOp)
    print(res)

    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    WriteMeshToXdmf(tempdir+"Test_ExtractElementsByImplicitZone.xdmf",res,PointFields=[res.originalIDNodes],PointFieldsNames=["originalIDNodes"] )
    print(tempdir)
    return "OK"


def CheckIntegrity_EnsureUniquenessElements(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    # this function must work
    EnsureUniquenessElements(mesh)

    tets = [[0,1,2,3],[3,1,0,2]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    # This function must not work
    try:
       EnsureUniquenessElements(mesh)
       return "No ok"
    except :
       pass
    return "ok"

def CheckIntegrity_GetElementsFractionInside(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles,CreateMeshOf

    meshTris = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    field = np.array([-1, -1, -1, 1])
    for name,elements in meshTris.elements.items():
        res = GetElementsFractionInside(field,meshTris.nodes,name,elements,range(elements.GetNumberOfElements()))
        print(res)

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,0] ]
    tets = [[0,1,2,3],]
    mesh3D = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)

    for name,elements in mesh3D.elements.items():
        res = GetElementsFractionInside(field,mesh3D.nodes,name,elements,range(elements.GetNumberOfElements()))
        print(res)

    return "ok"

def CheckIntegrity_GetDualGraph(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    dg, nused = GetDualGraph(res)

    return "ok"

def CheckIntegrity_ComputeMeshMinMaxLengthScale(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tets = [[0,1,2,3],[3,1,0,2]]
    mesh = CreateMeshOf(points,tets,ElementNames.Tetrahedron_4)
    print(ComputeMeshMinMaxLengthScale(mesh))
    PrintMeshInformation(mesh)
    return "ok"

def CheckIntegrity_ExtractElementsByMask(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("tri1").AddToTag(0)

    #tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([False,True]))
    #print(tri.connectivity)
    print(res.GetElementsOfType(ElementNames.Triangle_3).connectivity)

    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0],dtype=PBasicIndexType))
    tri = ExtractElementsByMask(res.GetElementsOfType(ElementNames.Triangle_3),np.array([0,1],dtype=bool))

    print(tri.connectivity)
    print(tri.originalIds)
    return "ok"
def CheckIntegrity(GUI=False):
    totest= [
    CheckIntegrity_EnsureUniquenessElements,
    CheckIntegrity_ExtractElementsByImplicitZone,
    CheckIntegrity_ExtractElementsByMask,
    CheckIntegrity_GetVolume,
    CheckIntegrity_GetDualGraph,
    CheckIntegrity_ComputeMeshMinMaxLengthScale,
    ]
    for f in totest:
        print("running test : " + str(f))
        res = f(GUI)
        if str(res).lower() != "ok":
            return "error in "+str(f) + " res"
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
