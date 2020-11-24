# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.Filters import ElementFilter, IntersectionElementFilter
import BasicTools.Containers.ElementNames as EN

def DistanceToSurface(mesh,surfMesh,out = None):
    from BasicTools.Containers.UnstructuredMeshFieldOperations import TransportPos
    from BasicTools.FE.FETools import PrepareFEComputation
    from BasicTools.FE.Fields.FEField import FEField

    space, numberings, offset, NGauss = PrepareFEComputation(surfMesh,numberOfComponents=1)


    tspace, tnumbering,_,_ = PrepareFEComputation(mesh,numberOfComponents=1)
    tnumbering = tnumbering[0]

    names = ["x","y","z"]
    InterfacePosFields = TransportPos(surfMesh,mesh,tspace,tnumbering)
    MeshPosFields = np.array([ FEField("pos_"+names[x],mesh, space=tspace, numbering=tnumbering,data=mesh.nodes[:,x]) for x in [0,1,2] ])

    if out is None:
        res = np.sqrt(np.sum((InterfacePosFields - MeshPosFields)**2)).data
        return res
    else:
        out[:] = np.sqrt(np.sum((InterfacePosFields - MeshPosFields)**2)).data
        return out

def Redistance(mesh,phi,out=None):

    IGTM = IGToMesh(mesh,phi )
#    IGTM.meshPartition = False
    interfaceMesh = IGTM.ComputeInterfaceMesh()
    sign = np.sign(phi)
    res = DistanceToSurface(mesh,interfaceMesh,out=out)
    res *= sign
    return res

class IGToMesh:
    def __init__(self,imesh=None,phi=None):
        self.inputMesh = imesh
        self.SetPhi(phi)
#        self.meshPartition = False

    def SetPhi(self,phi,snapValues=True, snapTol = 1e-5):
        if phi is None:
            self.phi = None
            return
        self.phi = np.copy(phi)
        phimin = min(self.phi)
        phimax = max(self.phi)
        self.phi[abs(self.phi) < abs(phimax-phimin)*snapTol] = 0.

    def ComputeInterfaceMesh(self):
        from sklearn.decomposition import PCA
        from scipy.spatial import Delaunay

        inputNodes = self.inputMesh.nodes
        self.volMesh = UnstructuredMesh()
        self.volMesh.CopyProperties(self.inputMesh)
        phi = self.phi
        pid = dict()
        elems = set()
        cpt = 0
        vcpt = 0
        xyz = []
        vxyz = []
        omesh = UnstructuredMesh()
        ef1 = ElementFilter(self.inputMesh,dimensionality=3,zones=[lambda x: phi<=0] )
        ef2 = ElementFilter(self.inputMesh,dimensionality=3,zones=[lambda x: -phi<=0] )
        ef1.zoneTreatment = "leastonenode"
        ef2.zoneTreatment = "leastonenode"
        ef = IntersectionElementFilter(self.inputMesh,[ef1,ef2])

        def AddElement(pos_id,neg_id,td,p):
                if pos_id is not None:
                    if np.cross(xyz[p[2]]-xyz[p[1]],xyz[p[0]]-xyz[p[1]]).dot(inputNodes[pos_id,:]-pcontrol_point ) > 0:
                        td.AddNewElement(p,0)
                    else:
                        td.AddNewElement(np.flip(p),0)
                elif neg_id is not None:
                    if np.cross(xyz[p[2]]-xyz[p[1]],xyz[p[0]]-xyz[p[1]]).dot(inputNodes[neg_id,:]-ncontrol_point ) > 0:
                        td.AddNewElement(np.flip(p),0)
                    else:
                        td.AddNewElement(p,0)
                else:
                    raise

        for name, data, ids in ef:
            for eid in ids:
                cutpoints = []
                neg_pointsO = []
                pos_pointsO = []
                neg_pointsG = []
                pos_pointsG = []

                lcoon = data.connectivity[eid,:]
                #if zero on the a point
                neg_id = None
                pos_id = None
                for n in lcoon:
                    if phi[n] == 0:
                        if n not in pid:
                            nid = cpt
                            xyz.append(inputNodes[n,:])
                            pid[n] = cpt
                            cpt +=1
                        else:
                            nid = pid[n]
                        cutpoints.append(nid)
                        neg_pointsO.append(n)
                        pos_pointsO.append(n)
                    elif phi[n] < 0:
                        neg_id = n
                        neg_pointsO.append(n)
                    elif phi[n] > 0:
                        pos_id = n
                        pos_pointsO.append(n)

                #if zero on a bars
                for n,d in EN.faces2[name]:
                    llcon=lcoon[d[0:2]]
                    lphi = phi[llcon]
                    pphi = lphi > 0
                    nphi = lphi < 0
                    if np.any(pphi) and np.any(nphi):
                        t = tuple(np.sort(llcon))
                        if t in pid:
                            nid,vnid = pid[t]
                        else:
                            a0 =  lphi[1]/(lphi[1] - lphi[0])
                            a1 = 1-a0
                            nid = cpt
                            vnid = vcpt
                            newpoint = inputNodes[llcon[1],:]*a1 +inputNodes[llcon[0],:]*a0
                            xyz.append(newpoint)
                            vxyz.append(newpoint)
                            pid[t] = (cpt,vcpt)
                            cpt +=1
                            vcpt +=1

                        cutpoints.append(nid)
                        neg_pointsG.append(vnid)
                        pos_pointsG.append(vnid)

                if len(cutpoints) < 3:
                    continue

                scp = tuple(np.sort(cutpoints))
                if scp in elems:
                    continue

                elems.update(scp)
                X = np.array([xyz[p] for p in cutpoints])

                pca = PCA(n_components=2)
                XX = pca.fit_transform(X)
                # use this to check the orientation

                if pos_id is not None:
                    pcontrol_point =  pca.inverse_transform(pca.transform( inputNodes[[pos_id],:] ))[0]
                elif neg_id is not None:
                    ncontrol_point =  pca.inverse_transform(pca.transform( inputNodes[[neg_id],:] ))[0]

                if len(cutpoints)==3 or len(cutpoints)==4:
                    ang = [np.arctan2(x,y) for x,y in XX]
                    s = np.argsort(ang)
                    cutpoints2 = [cutpoints[x] for x in s]
                    if len(cutpoints)==3:
                        td = omesh.elements.GetElementsOfType(EN.Triangle_3)
                    else:
                        td = omesh.elements.GetElementsOfType(EN.Quadrangle_4)
                    AddElement(pos_id,neg_id,td,cutpoints2)
                else:
                    tri = Delaunay(XX)
                    td = omesh.elements.GetElementsOfType(EN.Triangle_3)
                    for simple in tri.simplices:
                        p = [cutpoints[x] for x in simple]
                        AddElement(pos_id,neg_id,td,p)


#                this algorithm is local and generate incompatible meshes of tets
#                if self.meshPartition:
#                    P = np.vstack( (self.inputMesh.nodes[neg_pointsO],
#                                    [vxyz[x] for x in neg_pointsG]))
#                    index = np.hstack( (neg_pointsO,
#                                    [self.inputMesh.GetNumberOfNodes()+x for x in neg_pointsG]))
#                    tets = Delaunay(P)
#                    td = self.volMesh.elements.GetElementsOfType(EN.Tetrahedron_4)
#                    for simple in tets.simplices:
#                        p = [index[x] for x in simple]
#                        n= td.AddNewElement(p,0)
#                        td.tags.CreateTag("Inside",False).AddToTag(n-1)
#
#                    P = np.vstack( (self.inputMesh.nodes[pos_pointsO],
#                                    [vxyz[x] for x in pos_pointsG]))
#                    index = np.hstack( (pos_pointsO,
#                                    [self.inputMesh.GetNumberOfNodes()+x for x in pos_pointsG]))
#                    tets = Delaunay(P)
#                    td = self.volMesh.elements.GetElementsOfType(EN.Tetrahedron_4)
#                    for simple in tets.simplices:
#                        p = [index[x] for x in simple]
#                        n= td.AddNewElement(p,0)
#                        td.tags.CreateTag("Outside",False).AddToTag(n-1)
#
#        self.volMesh.PrepareForOutput()

#        if self.meshPartition:
#            vmesh = UnstructuredMesh()
#            #print(self.inputMesh.nodes)
#            #print( np.array(vxyz) )
#            vmesh.nodes = np.vstack( (self.inputMesh.nodes, np.array(vxyz) ) )
#            vmesh.MergeElements(self.inputMesh,force=True)
#            vmesh.MergeElements(self.volMesh,force=True)
#            vmesh.PrepareForOutput()
#            self.volMesh = vmesh

        omesh.nodes = np.array(xyz)
        omesh.GenerateManufacturedOriginalIDs()
        omesh.PrepareForOutput()

        return  omesh

def CheckIntegrity(GUI=False):

    import BasicTools.Containers.UnstructuredMeshCreationTools as UMCT
    from BasicTools.ImplicitGeometry.ImplicitGeometryObjects import ImplicitGeometrySphere
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.IO.XdmfWriter import  WriteMeshToXdmf

    print("Creating mesh")

    n = 15
    nn = n-1
    myMesh = UMCT.CreateCube(dimensions=[n,n,n], origin=[0,0,0], spacing=[1./nn,1./nn,1./nn], ofTetras=True)

    #print(myMesh)

    print("Creating phi")
    sB0 = ImplicitGeometrySphere(radius=0.5,center=[0,0,0])
    phi = sB0(myMesh)

    oo = IGToMesh(myMesh,phi )

    omesh =  oo.ComputeInterfaceMesh()
    print(omesh)

    tempdir = TestTempDir.GetTempPath()

    WriteMeshToXdmf(tempdir+'SphereIsoZero.xdmf',omesh)
    print("Iso zero mesh in " +tempdir+'SphereIsoZero.xdmf')
    phi2 = phi*2
    NewPhi = Redistance(myMesh,phi2)

    #print("NewPhi")
    #print(NewPhi)

    WriteMeshToXdmf(tempdir+'SpherePhiNew.xdmf',myMesh,PointFieldsNames=["Phi","PhiX2","NewPhi"],PointFields=[phi,phi2,NewPhi])
    print("phis on original mesh in " +tempdir+'SpherePhiNew.xdmf')

    error = max(np.abs(phi-NewPhi))
    if error > 1e-2:
        return "KO"

    return  'OK'

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=True))
