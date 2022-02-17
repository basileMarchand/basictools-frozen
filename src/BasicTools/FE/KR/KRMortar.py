# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import numpy as np

from scipy.sparse import coo_matrix, diags

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
from BasicTools.FE.IntegrationsRules import LagrangeIsoParam as LagrangeIsoParam
from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.Linalg.Transform import Transform
from BasicTools.FE.KR.KRBase import KRBaseVector
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Helpers.BaseOutputObject import froze_it

@froze_it
class KRMortar(KRBaseVector):
    def __init__(self):
        super(KRMortar,self).__init__()
        self.type = "Mortar"
#        self.inttype = "global" # ["PI","first","second","global","experimental"]

        self.useSurface = "first_surface" # [ "mean_surface", "first_surface", "second_surface","flat"]
        self.onII =[]
        self.ang = 30/180.*np.pi  # tolerance angle
        self.originSystem.keepOrthogonal = True
        self.originSystem.keepNormalised = True
        self.targetSystem.keepOrthogonal = True
        self.targetSystem.keepNormalised = True

        self._debug_IntegrationMesh = None

    def From(self,offset=None,first=None,second=None):
        if offset is not None:
            self.originSystem.SetOffset(offset)
        if first is not None:
            self.originSystem.SetFirst(first)
        if second is not None:
            self.originSystem.SetSecond(second)
        return self

    def To(self,offset=None,first=None,second=None):
        if offset is not None:
            self.targetSystem.SetOffset(offset)
        if first is not None:
            self.targetSystem.SetFirst(first)
        if second is not None:
            self.targetSystem.SetSecond(second)
        return self

    def SideI(self,zone):
        return self.On(zone)

    def SideII(self,zone):
        if type(zone) is list:
            self.onII.extend(zone)
        else:
            self.onII.append(zone)
        self.onII = list(set(self.onII))
        return self

    def GenerateEquations(self,meshI,fields,CH=None,meshII=None,fieldsII=None):

        CH = self._GetConstraintHolder(CH)

        if meshII is None:
            meshII = meshI
            fieldsII = fields

        fieldDicI = {f.name:f for f in fields }
        fieldDicII = {f.name:f for f in fields }

        offsetsI , self.__fieldOffsetsI, totalNumberOfDofsI  = self._ComputeOffsets(fields)
        offsetsII, self.__fieldOffsetsII, totalNumberOfDofsII = self._ComputeOffsets(fieldsII)

        def AddElementToViz(visuMesh, points, t,tag=None): # pragma: no cover
            if visuMesh is None:
                return
            n = visuMesh.GetNumberOfNodes()

            visuMesh.SetNodes(np.concatenate((visuMesh.nodes, points), axis=0),generateOriginalIDs=True)
            data = visuMesh.GetElementsOfType(t)
            cpt = data.AddNewElement(range(n,n+data.GetNumberOfNodesPerElement()),0)
            if tag is not None:
                data.GetTag(tag).AddToTag(cpt-1)

        bar_p, bar_w =  LagrangeIsoParam[ElementNames.Bar_2]
        bar_space_ipvalues = LagrangeSpaceGeo[ElementNames.Bar_2].SetIntegrationRule(bar_p,bar_w)

        tri_p,tri_w =  LagrangeIsoParam[ElementNames.Triangle_3]
        tri_space_ipvalues = LagrangeSpaceGeo[ElementNames.Triangle_3].SetIntegrationRule(tri_p,tri_w)

        tet_p,tet_w =  LagrangeIsoParam[ElementNames.Tetrahedron_4]
        tet_space_ipvalues = LagrangeSpaceGeo[ElementNames.Tetrahedron_4].SetIntegrationRule(tet_p,tet_w)

        meshI_IPoints = []
        meshII_IPoints = []
        weights_IPoints = []
        mesh_IElement = []
        mesh_IIElement = []

        # reference systems
        OS = self.originSystem.GetOrthoNormalBase()
        TS = self.targetSystem.GetOrthoNormalBase()

        # working memory
        # for bars
        integrationBar = np.zeros((2,3))
        iwBar = np.zeros_like(bar_w)
        ipBar = np.zeros((len(bar_w),3))
        # for triangles
        integrationTri = np.zeros((3,3))
        iwTri = np.zeros_like(tri_w)
        ipTri = np.zeros((len(tri_w),3))
        # for tetra
        integrationTet = np.zeros((4,3))
        iwTet = np.zeros_like(tet_w)
        ipTet = np.zeros((len(tet_w),3))

        # loop over every element in the two meshes
        from BasicTools.Helpers.ProgressBar import printProgressBar
        #computation of the integration points
        for name1,data1,ids1 in ElementFilter(meshI,tags=self.on):
            cpt = 0
            for elid1 in ids1:
                printProgressBar(cpt,len(ids1))
                cpt += 1
                #pos of node in the real spaces
                _nodes1 = meshI.nodes[data1.connectivity[elid1],:]
                #normal in the real space

                if ElementNames.dimension[name1] < 3:
                    _normal1,_baricentre1, incellvector1  = self.__ConputeNormal(meshI,name1,data1,elid1)
                    #baricenter in the projection space
                    #baricentre1 = OS.ApplyTransform(_baricentre1)
                    #normal in the projection space
                    normal1 = OS.ApplyTransformDirection(_normal1)
                #nodes in the projection space
                nodes1 = OS.ApplyTransform(_nodes1)
                #bounding box
                min1 = np.min(nodes1,axis=0)
                max1 = np.max(nodes1,axis=0)

                for name2,data2,ids2 in ElementFilter(meshII,tags=self.onII,dimensionality=ElementNames.dimension[name1]):
                    for elid2 in ids2:

                        _nodes2 = meshII.nodes[data2.connectivity[elid2],:]
                        nodes2 = TS.ApplyTransform(_nodes2)

                        # check if bounding box intersect
                        min2 = np.min(nodes2,axis=0)
                        max2 = np.max(nodes2,axis=0)

                        intersection = True
                        tol = 0.1*max( np.max(max1-min1),np.max(max2-min2))
                        for i in range(len(min1)):
                             intersection = intersection and ((min2[i] <= max1[i]+tol and min2[i] >= min1[i]-tol) or (
                                max2[i] <= max1[i]+tol and max2[i] >= min1[i]-tol) or (
                                min1[i] <= max2[i]+tol and min1[i] >= min2[i]-tol) or (
                                max1[i] <= max2[i]+tol and max1[i] >= min2[i]-tol))

                        if not intersection :
                            continue

                        if ElementNames.dimension[name2] < 3:
                            _normal2,_baricentre2, incellvector2 = self.__ConputeNormal(meshII,name2,data2,elid2)
                            normal2 = TS.ApplyTransformDirection(_normal2)
                            # only for 1d and 2d elements
                            # check normal are aligned using the self.ang
                            # the orientation is not important
                            angle = np.arccos(normal1.dot(normal2))
                            if abs(angle) > self.ang and abs(angle-np.pi) > self.ang:
                                continue

                        # we have intersection of 2 elements
                        if ElementNames.dimension[name2] == 1:
                            # element of dimension 1 (bars)
                            if self.useSurface == "mean_surface":
                                lineNormal = (normal1-normal2)/2
                                lineNormal /= np.linalg.norm(lineNormal)
                                base = np.array([-lineNormal[1], lineNormal[0], lineNormal[2]])
                                offset = (min1+min2+max1+max2)/4
                            # projection of the point in this line
                            elif self.useSurface == "first_surface":
                                base = np.array([-normal1[1],normal1[0],normal1[2]])
                                offset = (min1+max1)/2
                            elif self.useSurface == "second_surface":
                                base = np.array([-normal2[1],normal2[0],normal2[2]])
                                offset = (min2+max2)/2
                            else:
                                raise(Exception(f"Method {self.useSurface} not available!!!"))
                            second = base*1
                            for i,v in enumerate(second):
                                if v == 0:
                                    second[i] += 1
                                    break
                            else:
                                second[0] += 1

                            T = Transform(offset=offset,first=base, second=second )
                            proj1_ = T.ApplyTransform(nodes1)
                            proj2_ = T.ApplyTransform(nodes2)

                            proj1 = proj1_[:,0]
                            proj2 = proj2_[:,0]

                            #compute intersection
                            projmin = max(min(proj1),min(proj2))
                            projmax = min(max(proj1),max(proj2))
                            if projmax < projmin:
                                continue

                            integrationBar[0,0 ] = projmin
                            integrationBar[1,0 ] = projmax

                            # compute the coordinate of the integration points
                            for ip_nb in range(len(bar_w)):
                                Jack, Jdet, Jinv = bar_space_ipvalues.GetJackAndDetI(ip_nb,integrationBar)
                                iwBar[ip_nb] = Jdet*bar_w[ip_nb]
                                ipBar[ip_nb,:] = np.dot(bar_space_ipvalues.valN[ip_nb],integrationBar)

                            # if the integration bar is degenerated we skip it
                            if Jdet < 1e-8:
                                continue

                            # transfert the coordinate to the original meshes (1 and 2)
                            integration_coords_II = T.ApplyInvTransform(ipBar)

                            ## compute integration point in the meshI, meshII
                            original_int_points = OS.ApplyInvTransform(integration_coords_II)
                            target_int_points = TS.ApplyInvTransform(integration_coords_II)

                            if self._debug_IntegrationMesh is not None: # pragma: no cover
                                AddElementToViz(self._debug_IntegrationMesh,original_int_points,ElementNames.Bar_2,"eint_I")#,"I_E"+str(elid1)+"_E"+str(elid2) )
                                AddElementToViz(self._debug_IntegrationMesh,target_int_points,ElementNames.Bar_2,"eint_II")#,"I_E"+str(elid1)+"_E"+str(elid2) )

                            # Append points and weights
                            weights_IPoints.extend(iwBar)
                            mesh_IElement.extend(np.ones(len(bar_w))*elid1)
                            mesh_IIElement.extend(np.ones(len(bar_w))*elid2)
                            for ip_nb in range(len(bar_w)):
                                meshI_IPoints.append(original_int_points[ip_nb,:])
                                meshII_IPoints.append(target_int_points[ip_nb,:])

                                if self._debug_IntegrationMesh is not None:# pragma: no cover
                                    AddElementToViz(self._debug_IntegrationMesh,original_int_points[ip_nb:ip_nb+1,:],ElementNames.Point_1,"int_I")
                                    AddElementToViz(self._debug_IntegrationMesh,target_int_points[ip_nb:ip_nb+1,:],ElementNames.Point_1,"int_II")

                        elif ElementNames.dimension[name2] == 2:
                            # element of dimension 2 (tris)
                            if self.useSurface == "flat":
                                lineNormal = [0,0,1]
                            elif self.useSurface == "mean_surface":
                                lineNormal = (normal1-normal2)/2
                                lineNormal /= np.linalg.norm(lineNormal)
                            # projection of the point in this line
                            elif self.useSurface == "first_surface":
                                lineNormal  = normal1
                            elif self.useSurface == "second_surface":
                                lineNormal  = normal2
                            #T = Transform(first=[-lineNormal[1], lineNormal[0], lineNormal[2]],second=[-lineNormal[2], lineNormal[1], lineNormal[0]])
                            T = Transform()
                            T.SetOpUsingThird(lineNormal)

                            proj1 = T.ApplyTransform(nodes1)
                            proj2 = T.ApplyTransform(nodes2)

                            inter = Intersection(proj1*[1,1,0],range(data1.GetNumberOfNodesPerElement()),
                                                 proj2*[1,1,0],range(data2.GetNumberOfNodesPerElement()),tol/100000.)

                            if len(inter) < 4 :
                                # insufficent point to build as 2D domain
                                continue
                            if len(inter) == 4 :

                                # we have a triangle (first point is repetead at the end)
                                center = inter[0,:]
                                inter = inter[1:3,:]
                            else:
                                center = np.sum(inter[0:-1,:],axis=0)/(inter.shape[0]-1)

                            # we use the center as the first point for all the triangles
                            integrationTri[0,:] = center

                            for i in range(inter.shape[0]):
                                # copy the coord to generate a triangle
                                integrationTri[1:3,:] = inter[i:i+2,:]
                                # compute the coordinate of the integration points
                                for ip_nb in range(len(tri_w)):
                                    Jack, Jdet, Jinv = tri_space_ipvalues.GetJackAndDetI(ip_nb,integrationTri)
                                    iwTri[ip_nb] = Jdet*tri_w[ip_nb]
                                    ipTri[ip_nb,:] = np.dot(tri_space_ipvalues.valN[ip_nb],integrationTri)

                                # if the integration triangle is degenerated we skip it
                                if Jdet < 1e-8:
                                    continue

                                # transfert the coordinate to the original meshes (1 and 2)
                                integration_coords_II = T.ApplyInvTransform(ipTri)
                                if self._debug_IntegrationMesh is not None:# pragma: no cover
                                    AddElementToViz(self._debug_IntegrationMesh,OS.ApplyInvTransform(T.ApplyInvTransform(integrationTri)),ElementNames.Triangle_3,"eint_I")
                                    AddElementToViz(self._debug_IntegrationMesh,TS.ApplyInvTransform(T.ApplyInvTransform(integrationTri)),ElementNames.Triangle_3,"eint_II")

                                ## compute integration point in the meshI, meshII

                                original_int_points = OS.ApplyInvTransform(integration_coords_II)
                                target_int_points = TS.ApplyInvTransform(integration_coords_II)
                                if self._debug_IntegrationMesh is not None:# pragma: no cover
                                    AddElementToViz(self._debug_IntegrationMesh,original_int_points,ElementNames.Triangle_3,"eint_I")#,"I_E"+str(elid1)+"_E"+str(elid2) )
                                    AddElementToViz(self._debug_IntegrationMesh,target_int_points,ElementNames.Triangle_3,"eint_II")#,"I_E"+str(elid1)+"_E"+str(elid2) )

                                # Append points and weights
                                weights_IPoints.extend(iwTri)
                                mesh_IElement.extend(np.ones(len(tri_w))*elid1)
                                mesh_IIElement.extend(np.ones(len(tri_w))*elid2)
                                for ip_nb in range(len(tri_w)):
                                    meshI_IPoints.append(original_int_points[ip_nb,:])
                                    meshII_IPoints.append(target_int_points[ip_nb,:])
                                    AddElementToViz(self._debug_IntegrationMesh,original_int_points[ip_nb:ip_nb+1,:],ElementNames.Point_1,"int_I")
                                    AddElementToViz(self._debug_IntegrationMesh,target_int_points[ip_nb:ip_nb+1,:],ElementNames.Point_1,"int_II")

                        elif ElementNames.dimension[name2] == 3:
                            status, points, tets = IntersectionOf2CovexHulls(nodes1,nodes2)
                            if status == False:
                                continue

                            for i in range(tets.shape[0]):
                                # copy the coord to generate a tet
                                integrationTet[:,:] = points[tets[i,:],:]
                                for ip_nb in range(len(tet_w)):
                                    Jack, Jdet, Jinv = tet_space_ipvalues.GetJackAndDetI(ip_nb,integrationTet)
                                    iwTet[ip_nb] = abs(Jdet)*tet_w[ip_nb]
                                    ipTet[ip_nb,:] = np.dot(tet_space_ipvalues.valN[ip_nb],integrationTet)

                                # if the integration triangle is degenerated we skip it

                                if abs(Jdet) < 1e-8:
                                    continue

                                # transfert the coordinate to the original meshes (1 and 2)
                                if self._debug_IntegrationMesh is not None: # pragma: no cover
                                    AddElementToViz(self._debug_IntegrationMesh,OS.ApplyInvTransform(integrationTet),ElementNames.Tetrahedron_4,"eint_I_"+str(elid1)+"_"+str(elid2))
                                    #AddElementToViz(self._debug_IntegrationMesh,TS.ApplyInvTransform(integrationTet),ElementNames.Tetrahedron_4,"eint_II"+str(elid1)+"_"+str(elid2))
                                original_int_points = OS.ApplyInvTransform(ipTet)

                                target_int_points = TS.ApplyInvTransform(ipTet)
                                if self._debug_IntegrationMesh is not None: # pragma: no cover
                                    AddElementToViz(self._debug_IntegrationMesh,original_int_points,ElementNames.Point_1,"Left")#,"I_E"+str(elid1)+"_E"+str(elid2) )
                                    AddElementToViz(self._debug_IntegrationMesh,target_int_points,ElementNames.Point_1,"Right")#,"I_E"+str(elid1)+"_E"+str(elid2) )
                                # Append points and weights
                                weights_IPoints.extend(iwTet)

                                mesh_IElement.extend(np.ones(len(tet_w))*elid1)
                                mesh_IIElement.extend(np.ones(len(tet_w))*elid2)
                                for ip_nb in range(len(tet_w)):
                                    meshI_IPoints.append(original_int_points[ip_nb,:])
                                    meshII_IPoints.append(target_int_points[ip_nb,:])

                continue

            printProgressBar(len(ids1),len(ids1))

        if self._debug_IntegrationMesh is not None :
            self._debug_IntegrationMesh.PrepareForOutput()
            self._debug_IntegrationMesh.GenerateManufacturedOriginalIDs()

        weights_IPoints = np.array(weights_IPoints)
        meshI_IPoints = np.array(meshI_IPoints)
        meshII_IPoints= np.array(meshII_IPoints)
        mesh_IElement = np.array(mesh_IElement)
        mesh_IIElement = np.array(mesh_IIElement)

        totalNumberOfIP = weights_IPoints.shape[0]

        if len(meshI_IPoints) == 0:
            print("Warning! -> Zero elements in contact")
            return

        # need to code the transfert of the field to the integration points meshI
        from BasicTools.Containers.UnstructuredMeshFieldOperations import GetFieldTransferOp
        meshIOps = {}

        # we must apply the offsets to the operators  to build correctly the coupling terms.
        for f, offset  in zip(fields,offsetsI):
            op, status = GetFieldTransferOp(f,meshI_IPoints,method="Interp/Clamp",elementfilter=ElementFilter(meshI,tags=self.on) )
            # apply the offset to the op
            op.col += offset
            op.resize( (totalNumberOfIP,totalNumberOfDofsI) )
            meshIOps[f.name] = (op, status)

        # need to code the transfert of the field to the integration points meshII
        meshIIOps = {}
        for f,offset in zip(fieldsII,offsetsII):
            op, status = GetFieldTransferOp(f,meshII_IPoints,method="Interp/Clamp",elementfilter=ElementFilter(meshII,tags=self.onII) )
            op.col += offset
            op.resize( (totalNumberOfIP,totalNumberOfDofsII) )
            meshIIOps[f.name] = (op, status)

        diagWMatrix = diags(weights_IPoints)

        for arg in self.args:
            fieldI = fieldDicI.get(arg)
            fieldII = fieldDicII.get(arg)

            opI = meshIOps[arg][0]
            opII = meshIIOps[arg][0]

            A = diagWMatrix.dot(opI).T.dot(opI)
            B = diagWMatrix.dot(opI).T.dot(opII)

            ## differents types of behabiours
            if meshII is meshI:
                ## if 2 sides are in the same mesh (monolitic problem)
                ## the we add the constraint to the system
                data = (A-B)
                nums= np.where(np.sum(np.abs(data),axis=1)!= 0)[0]
                CH.AddEquationsFromOperatorAAndb((data[nums,:]).tocoo())
            else:
                pass
                ## in this case we are in "advance mode"
                ## the user is responsible of using the operators

        self.weights_IPoints = weights_IPoints
        self.meshIOps = meshIOps
        self.meshIIOps = meshIIOps
        return CH

    def __computeNormalSurface(self,nodes,conn):
        baricentre = np.sum(nodes[conn ,:],axis=0)/len(conn)
        edgeVector = nodes[conn[0],:] - nodes[conn[1],:]
        planeVector = baricentre - nodes[conn[1],:]
        normal = np.cross(edgeVector, planeVector)
        normal /= np.linalg.norm(normal)
        return normal, baricentre, planeVector/np.linalg.norm(planeVector)

    def __computeNormalEdge(self,nodes,conn):
        baricentre = np.sum(nodes[conn ,:],axis=0)/len(conn)
        planeVector = baricentre - nodes[conn[0],:]
        normal = np.array([planeVector[1], -planeVector[0],0])
        normal /= np.linalg.norm(normal)
        return normal, baricentre, planeVector/np.linalg.norm(planeVector)

    def __ConputeNormal(self,submesh,name,data,elid):
        if ElementNames.dimension[name] == 1:
            return self.__computeNormalEdge(submesh.nodes,data.connectivity[elid,:])
        elif ElementNames.dimension[name] == 2:
            return self.__computeNormalSurface(submesh.nodes,data.connectivity[elid,:])
        else: # pragma: no cover
            raise Exception(" Error ")

def AreCCW(p1,p2,p3):
    return bool((np.sign(np.cross(p2-p1,p3-p2)[2])+1)/2)

def Append(l,point,tol):
    p = np.array(point)
    if len(l):
        if np.linalg.norm(l[-1] - p ) > tol:
            l.append(p)
    else:
        l.append(point)

def CheckTriWinding(points, tri, allowReversed):
    trisq = np.ones((3,3))
    trisq[:,0:2] = np.array(points[tri,0:2])
    detTri = np.linalg.det(trisq)
    if detTri < 0.0:
        if allowReversed:
            #print(tri)
            a = tri[2]
            tri[2] = tri[1]
            tri[1] = a
            return tri
            a = trisq[2,:].copy()
            trisq[2,:] = trisq[1,:]
            trisq[1,:] = a

        else: # pragma: no cover
            raise ValueError("triangle has wrong winding direction")
    return tri


def Intersection(points1, _poly1,points2,_poly2, tol):
    ## cut poly1 by poly2

    poly1 = list(_poly1)
    poly1 = CheckTriWinding(points1,poly1,True)


    poly1 = [ points1[x,:] for x in poly1 ]

    poly2 = list(_poly2)
    poly2 = CheckTriWinding(points2,poly2,True)
    poly2.append(_poly2[0])

    res = []
    for cuts in range(len(poly2)-1):
        cp0 = points2[poly2[cuts]]
        cp1 = points2[poly2[cuts+1]]
        res = []

        Append(poly1,poly1[0],tol)
        for s in range(len(poly1)-1):
            sp0 = poly1[s]
            sp1 = poly1[s+1]
            if AreCCW(cp0,cp1,sp0):
                # point inside keep the point
                Append(res,sp0,tol)

            if AreCCW(cp0,cp1,sp0) != AreCCW(cp0,cp1,sp1):
                # segment must be cutted by cutter
                # keep the intersection
                x1 = cp0[0]
                y1 = cp0[1]
                x2 = cp1[0]
                y2 = cp1[1]

                x3 = sp0[0]
                y3 = sp0[1]
                x4 = sp1[0]
                y4 = sp1[1]
                if (( x1-x2 )*(y3-y4)-(y1-y2)*(x3-x4) ) == 0:
                    px = x4
                    py = y4
                else:
                    px = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4)) /(( x1-x2 )*(y3-y4)-(y1-y2)*(x3-x4) )
                    py = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)) /(( x1-x2 )*(y3-y4)-(y1-y2)*(x3-x4) )

                px = np.clip(px, min(x3,x4),max(x3,x4))
                py = np.clip(py, min(y3,y4),max(y3,y4))

                Append(res,[px,py,0],tol)
            if AreCCW(cp0,cp1,sp1):
                # the last point is treated by the next interation
                Append(res,sp1,tol)
        poly1 = res
        if len(poly1) < 3:
            break
    if len(res):
        Append(res,res[0],tol)
    return np.array(res)


def IntersectionOf2CovexHulls(pointsI,pointsII):
    from scipy.spatial import ConvexHull, HalfspaceIntersection
    from scipy.optimize import linprog

    # computation of the convex Hulls
    CHI = ConvexHull(points=pointsI)
    CHII = ConvexHull(points=pointsII)

    # computation of a interior_point
    # using method from
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.HalfspaceIntersection.html
    halfspaces = np.vstack((CHI.equations, CHII.equations))
    norm_vector = np.reshape(np.linalg.norm(halfspaces[:, :-1], axis=1), (halfspaces.shape[0], 1))
    c = np.zeros((halfspaces.shape[1],))
    c[-1] = -1
    A = np.hstack((halfspaces[:, :-1], norm_vector))
    b = - halfspaces[:, -1:]
    res = linprog(c, A_ub=A, b_ub=b, bounds=(None, None))

    # no intersection return false
    if res.success == False:
        return False, None, None

    x = res.x[:-1]

    # computation of the intersection
    try:
        HI = HalfspaceIntersection(halfspaces,interior_point = x, incremental=False)
    except:
        return False, None, None

    # we use all the intersection point to build the intersected convex hull
    CHIntersection = ConvexHull(points=HI.intersections)

    nbpoints = HI.intersections.shape[0]
    if nbpoints == 4:
        """ we know the only posibility in 3D is a tetra """
        points = HI.intersections
        tets = np.array([np.hstack((CHIntersection.simplices[0],[x  for x in range(4) if x not in CHIntersection.simplices[0]]))])
    else:
        nbpoints += 1
        points = np.vstack((HI.intersections,np.mean(HI.intersections,axis=0)))

        nbpoints = points.shape[0]
        tets = np.empty((len(CHIntersection.simplices),4),dtype=int)
        tets[:,0:3] = CHIntersection.simplices
        tets[:,3] = nbpoints -1

    ## verification of the order of the tets
    for i in range(tets.shape[0]):

        p0 = points[tets[i,0],:]
        p1 = points[tets[i,1],:]
        p2 = points[tets[i,2],:]
        p3 = points[tets[i,3],:]
        if np.dot(p3-(p0+p1+p2)/3, np.cross(p1-p0,p2-p1)) < 0 :
            # swap if the forth point is on the wrong side
            tets[i,1], tets[i,0] = tets[i,0], tets[i,1]

    return True, points, tets

def CheckIntegrity3D(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    nx,ny, nz = 1, 1, 1
    meshI = CreateCube(dimensions=[nx+1,ny+1,nz+1],origin=[0,0,0], spacing=[10./nx,10./ny,1/nz], ofTetras=True)
    meshII = CreateCube(dimensions=[nx+1,ny+1,nz+1],origin=[5,0,0], spacing=[1./nx,10./ny,10/nz], ofTetras=True)

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1
    from BasicTools.Actions.OpenInParaView import OpenInParaView

    numberingI = ComputeDofNumbering(meshI,LagrangeSpaceGeo,fromConnectivity=True)
    numberingII = ComputeDofNumbering(meshII,LagrangeSpaceGeo,fromConnectivity=True)

    tag1 = "3D"
    tag2 = "3D"

    space = LagrangeSpaceP1
    unknownNames = ["T"]
    unknownNamesII = ["T"]

    fieldsI = []
    for name in unknownNames:
        fieldsI.append(FEField(name,meshI,space,numberingI))

    fieldsII = []
    for name in unknownNamesII:
        fieldsII.append(FEField(name,meshII,space,numberingII))


    obj = KRMortar()
    obj.AddArg("T")
    obj.From()
    obj.To(offset=[0,0,0]  )
    obj.SideI(tag1)
    obj.SideII([tag2])
    obj._debug_IntegrationMesh = UnstructuredMesh()

    print(obj)
    print(obj.GenerateEquations(meshI,fieldsI,meshII=meshII,fieldsII=fieldsII))
    if GUI: # pragma: no cover
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=meshI)
        OpenInParaView(mesh=meshII)
        OpenInParaView(mesh=obj._debug_IntegrationMesh)

    return "ok"

def CheckIntegrityIntersectionConvexHull3D(GUI=False):

    from scipy.spatial import ConvexHull

    pointsI = np.array([[0.0,0.0,0.0],
                       [ 1.0,0.0,0.0],
                       [ 0.0,1.0,0.0],
                       [ 0.0,0.0,1.0]])

    pointsII = np.array([[0.5,0.0,0.1],
                       [1.5,0,0],
                       [1.5,1,0],
                       [1.5,0,1]])

    status, points, tets = IntersectionOf2CovexHulls(pointsI,pointsII)

    if status == False:# pragma: no cover
        raise Exception("Error in the detection of the intersection")

    mesh = UnstructuredMesh()
    mesh.nodes = np.vstack((pointsI,pointsII,points))
    elements = mesh.GetElementsOfType(ElementNames.Triangle_3)

    CHI = ConvexHull(points=pointsI)
    CHII = ConvexHull(points=pointsII)
    CHintersection = ConvexHull(points=points)

    cpt = 0
    for simplex in CHI.simplices:
        elements.AddNewElement(simplex,cpt)
        cpt += 1

    for simplex in CHII.simplices:
        elements.AddNewElement(simplex+pointsI.shape[0],cpt)
        cpt += 1


    if tets.shape[0] != 1:# pragma: no cover
        raise Exception("ERROR! Wrong number of tetras in the intersection ")

    elements = mesh.GetElementsOfType(ElementNames.Tetrahedron_4)
    elements.connectivity = tets+pointsI.shape[0]+pointsII.shape[0]
    elements.cpt = tets.shape[0]

    mesh.GenerateManufacturedOriginalIDs()
    mesh.PrepareForOutput()


    if GUI:# pragma: no cover
        from BasicTools.Containers.vtkBridge import PlotMesh
        PlotMesh(mesh)
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh=mesh)

    return "ok"

def CheckIntegrityIntersection(GUI=False):

    P = np.array([[0.,0,0],
                   [5.,0,0],
                   [0.,5,0],
                   [0,6,0]])

    data = Intersection(P, [0,1,3],P,[0,1,2],tol=0.0001)
    center = np.sum(data,axis=0)/(data.shape[0]-1)

    result = (center,data)

    P = np.array([[0.,0,0],
                   [4.,0,0],
                   [0.,3,0],
                   [2.,-1,0],
                   [3.,2,0],
                   [-1.,3,0]])


    data = Intersection(P, [0,1,2],P,[3,4,5],tol=0.0000001)
    center = np.sum(data[0:-1,:],axis=0)/(data.shape[0]-1)
    result = (center,data)
    return "ok"

def CheckIntegrity1DInterface2Meshes(GUI=False):
    """
    CheckIntegrity for
    Mesh 1

       y ^
         |
       1 |x       2 x
         ||
         ||       3 x
         ||
       0 -x-------4-x----> x

   y ^
         |
       1 |       2 x
         |         |
         |       3 x
         |         |
       0 --------4-x----> x


    """
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1

    ps = np.array([[0.,0.,0],
                    [0.,1.,0],
                    [1.,1,0],
                    [1.,0.5,0],
                    [1.,0.,0],
                    ])

    es = np.array([[0,1]])

    meshI = CreateMeshOf(ps,es, ElementNames.Bar_2)

    ps = np.array([ [1.,1,0],
                    [1.,0.5,0],
                    [1.,0.,0],
                    ])

    es = np.array([[0,1],[1,2]])
    meshII = CreateMeshOf(ps,es, ElementNames.Bar_2)

    tag1 = "Left"
    tag2 = "Right"
    meshI.GetElementsOfType(ElementNames.Bar_2).tags.CreateTag(tag1).SetIds(0)
    meshII.GetElementsOfType(ElementNames.Bar_2).tags.CreateTag(tag2).SetIds([0,1])
    print(meshI)
    print(meshII)

    numberingI = ComputeDofNumbering(meshI,LagrangeSpaceGeo)
    numberingII = ComputeDofNumbering(meshII,LagrangeSpaceGeo)
    print("ndofs I",numberingI.size)
    print("ndofs II",numberingII.size)

    print("------------------------")
    space = LagrangeSpaceP1
    unknownNames = ["T"]

    fieldsI = []
    fieldsII = []
    for name in unknownNames:
        fieldsI.append(FEField(name,meshI,space,numberingI))
        fieldsII.append(FEField(name,meshII,space,numberingII))

    obj = KRMortar()
    obj.AddArg("T")
    #obj.From()
    obj.To(offset=[1,0,0] )
    obj.SideI(tag1)
    obj.SideII(tag2)
    print(obj)
    obj.inttype = "PI"
    obj._debug_IntegrationMesh = UnstructuredMesh()
    CH = obj.GenerateEquations(meshI,fieldsI,CH=None,meshII=meshII,fieldsII=fieldsII)
    if GUI : # pragma: no cover
        from BasicTools.Containers.vtkBridge import PlotMesh
        PlotMesh(meshI)
        PlotMesh(meshII)
        PlotMesh(obj._debug_IntegrationMesh)

    return "ok"

def CheckIntegrity1DInterface(GUI=False):
    """
    CheckIntegrity for
    y ^
         |
       1 |x       2 x
         ||         |
         ||       3 x
         ||         |
       0 -x-------4-x----> x
    """
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1

    ps = np.array([[0.,0.,0],
                    [0.,1.,0],
                    [1.,1,0],
                    [1.,0.5,0],
                    [1.,0.,0],
                    ])

    es = np.array([[0,1],[2,3],[3,4]])

    mesh = CreateMeshOf(ps,es, ElementNames.Bar_2)

    tag1 = "Left"
    tag2 = "Right"
    mesh.GetElementsOfType(ElementNames.Bar_2).tags.CreateTag(tag1).SetIds(0)
    mesh.GetElementsOfType(ElementNames.Bar_2).tags.CreateTag(tag2).SetIds([1,2])
    print(mesh)

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)

    space = LagrangeSpaceP1
    unknownNames = ["T"]

    fields = []
    for name in unknownNames:
        fields.append(FEField(name,mesh,space,numbering))

    obj = KRMortar()
    obj.AddArg("T")
    obj.To(offset=[1,0,0] )
    obj.SideI(tag1)
    obj.SideII(tag2)
    obj._debug_IntegrationMesh = UnstructuredMesh()
    CH = obj.GenerateEquations(mesh,fields)
    CH.SetNumberOfDofs(5)
    print(CH.ToSparseFull().toarray())
    if GUI : # pragma: no cover
        from BasicTools.Containers.vtkBridge import PlotMesh
        PlotMesh(mesh)
        PlotMesh(obj._debug_IntegrationMesh)
    return "ok"

def CheckIntegrity2DScalar(GUI=False):
    from BasicTools.TestData import GetTestDataPath
    filename = GetTestDataPath()+"dent.msh"

    from BasicTools.IO.UniversalReader import ReadMesh
    import BasicTools.IO.GmshReader # to register the GmshReader
    mesh = ReadMesh(filename)
    #print(mesh.nodes.shape)
    nnodes = np.zeros((mesh.nodes.shape[0],3),dtype=float)
    nnodes[0:,:2] = mesh.nodes
    mesh.nodes = nnodes
    #print(mesh.nodes.shape)
    #raise

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)

    tag1 = "Left"
    tag2 = "Right"

    space = LagrangeSpaceP1
    unknownNames = ["T"]

    fields = []
    for name in unknownNames:
        fields.append(FEField(name,mesh,space,numbering))


    obj = KRMortar()
    obj.AddArg("T")
    obj.From(first=[-np.sin(np.pi/5),np.cos(np.pi/5),0 ], second=[-np.cos(np.pi/5),-np.sin(np.pi/5),0 ]  )
    angle = np.pi/5+0.01
    obj.To(first=[np.sin(angle),np.cos(angle),0 ], second=[-np.cos(angle),np.sin(angle),0 ]  )
    obj.SideI(tag1)
    obj.SideII(tag2)
    obj._debug_IntegrationMesh = UnstructuredMesh()
    print(obj)
    print(obj.GenerateEquations(mesh,fields))
    if GUI: # pragma: no cover
        from BasicTools.Containers.vtkBridge import PlotMesh
        PlotMesh(obj._debug_IntegrationMesh)

    return "ok"

def CheckIntegrity3DVector(GUI=False):
    from BasicTools.TestData import GetTestDataPath
    filename = GetTestDataPath()+"dent3D.msh"

    from BasicTools.IO.UniversalReader import ReadMesh
    import BasicTools.IO.GmshReader

    mesh = ReadMesh(filename)

    #print(mesh)
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP1

    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity=True)

    tag1 = "Left"
    tag2 = "Right"

    space = LagrangeSpaceP1
    unknownNames = ["u_0","u_1","u_2"]

    fields = []
    for name in unknownNames:
        fields.append(FEField(name,mesh,space,numbering))


    obj = KRMortar()
    obj.AddArg("u_0")
    obj.AddArg("u_1")
    obj.AddArg("u_2")
    obj.From(offset=[0,0,0], first=[-np.sin(np.pi/5),np.cos(np.pi/5),0 ], second=[-np.cos(np.pi/5),-np.sin(np.pi/5),0 ]  )
    angle = np.pi/5
    obj.To  (first=[np.sin(angle),np.cos(angle),0 ], second=[-np.cos(angle),np.sin(angle),0 ]  )
    obj.SideI(tag1)
    obj.SideII(tag2)
    obj._debug_IntegrationMesh = UnstructuredMesh()
    print(obj.GenerateEquations(mesh,fields))

    if GUI: # pragma: no cover
        from BasicTools.Containers.vtkBridge import PlotMesh
        PlotMesh(obj._debug_IntegrationMesh)

    return "ok"

def CheckIntegrity(GUI=False):
    func = [CheckIntegrity1DInterface, CheckIntegrity1DInterface2Meshes, CheckIntegrity2DScalar, CheckIntegrity3DVector, CheckIntegrityIntersection, CheckIntegrity3D, CheckIntegrityIntersectionConvexHull3D]
    for f in func:
        print("working on : ", f)
        res = f(GUI)
        print(res, f)
        if res != "ok": # pragma: no cover

            return res + "Error on " + str(f)
    return res

if __name__ == "__main__": # pragma: no cover
    CheckIntegrity(GUI=True)
