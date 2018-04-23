# -*- coding: utf-8 -*-
import numpy as np
import BasicTools.FE.ElementNames as EN
from  BasicTools.FE.UnstructuredMesh import AllElements

def Hash(space,j,lconn,name=None,i=None):
    T,n,h = space.dofAttachments[j]
    if h is not None :
        h = Hash(h,lconn[n])
        raise

    if T == "P":
        T = "P"
        n = lconn[n]
    elif T == "C":
        T = name
        n = i
    elif T == "E":
        raise
        edge = EN.Faces[name][n]
        T = edge[0]
        n = np.sort(lconn[edge[1]])
    T += ("_"+ str(i)) if space.classification.discontinuous else ""
    #print (T,n,h),
    return (T,n,h)


def ComputeDofNumbering(mesh,Space,dofs=None,tag=AllElements,sign=1,fromConnectivity=False):
    if fromConnectivity:
        if dofs is not None or tag != AllElements or sign != 1:
            raise(Exception("cant take dofs tag or sign different from the default values"))
        dofs = {}
        dofs["size"] = mesh.GetNumberOfNodes()
        dofs["dirichelet"] = -1
        dofs["almanac"] = {}
        dofs["doftopointLeft"] = slice(mesh.GetNumberOfNodes())
        dofs["doftopointRight"] = dofs["doftopointLeft"]
        for name,data in mesh.elements.items():
            dofs[name] = data.connectivity

        return dofs

    if dofs is None:
        dofs = {}
        dofs["size"] = 0
        dofs["dirichelet"] = -1
        dofs["almanac"] = {}
        dofs["Hash"] = Hash
        if sign == 1:
            cpt = 0
        else:
            cpt = -1

    else:
        if sign == 1:
          cpt= dofs["size"]
        else:
          cpt = dofs["dirichelet"]

    almanac = dofs["almanac"]

    for name,data in mesh.elements.items():
        sp = Space[name]

        if not name in dofs:
            dof = np.empty((data.GetNumberOfElements(),sp.GetNumberOfShapeFunctions()), dtype=int)
        else:
            dof = dofs[name]

        if tag is AllElements:
            fil = range(data.GetNumberOfElements())
        elif tag in data.tags:
            fil = data.GetTag(tag).GetIds()
        else:
            continue

        numberOfShapeFunctions = sp.GetNumberOfShapeFunctions()
        for i in fil:
            conn = data.connectivity[i,:]
            for j in range(numberOfShapeFunctions):
                h = Hash(sp,j,conn,name,i)

                if h in almanac:
                    d = almanac[h]
                else:
                    d = cpt
                    almanac[h] = d
                    cpt += 1*sign
                dof[i,j] = d

        dofs[name] = dof


    if sign ==1:
        #print("size : " + str(cpt) )
        dofs["size"] = cpt
        # two point to take into account:
        #     the numbering can be applied only to  a subdomain
        #     the space used can be bigger than the mesh (p2 for example)
        # so we need a two side extractor in the form of
        #
        # PointSizeVector[doftopointLeft] = solution[doftopointRight]
        #
        # this way we can have a extractor on P2 or P3 solution into a P1 mesh
        # the user is responsible to allocate and initialise (of zeros) the
        # PointSizeVector of the right size
        extractorLeftSide = np.empty(cpt,dtype=np.int)
        extractorRightSide = np.empty(cpt,dtype=np.int)

        tmpcpt = 0
        # if k[0] is 'P' then k[1] is the node number
        for k,v in almanac.items():
            if k[0] == 'P':
                extractorLeftSide[tmpcpt] = k[1]
                extractorRightSide[tmpcpt] = v
                tmpcpt += 1

        extractorLeftSide.resize(tmpcpt)
        extractorRightSide.resize(tmpcpt)
        dofs["doftopointLeft"] =   extractorLeftSide
        dofs["doftopointRight"] =  extractorRightSide


    else:
        dofs["dirichelet"] = cpt


    return dofs

def NodesPermutation(mesh,per):
    nnodes = mesh.nodes[per,:]

    perII =np.argsort(per)

    mesh.nodes = nnodes
    for tag in mesh.nodesTags:
        ids = tag.GetIds()
        nids = perII[ids]
        tag.SetIds(nids)

    for name,data in mesh.elements.items():
        data.connectivity = perII[data.connectivity]



def CheckIntegrity(GUI=False):
    from BasicTools.FE.UnstructuredMeshTools import CreateCube

    #res2 = CreateCube([100.,100.,100.],[-1.0,-1.0,-1.0],[2./46, 2./46,2./46])
    res2 = CreateCube([2.,2.,2.],[-1.0,-1.0,-1.0],[2./46, 2./46,2./46])
    #NodesPermutation(res2,np.array([0, 4, 6,2 ,1 ,5 ,7, 3]))

    print(res2)

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

    print("Numbering ")
    import time
    st = time.time()
    #numbering = ComputeDofNumbering(res2,LagrangeSpaceGeo,tag="X0",sign=-1)
    numbering = ComputeDofNumbering(res2,LagrangeSpaceGeo)
    numbering = ComputeDofNumbering(res2,LagrangeSpaceGeo,fromConnectivity=True)
    print(time.time()-st)

    print(numbering)
    for name,data in res2.elements.items():
        if name in numbering:
            print(name)
            print(data.connectivity)
            print(numbering[name])


    import BasicTools.FE.ElementNames as EN
    print(res2.GetElementsOfType(EN.Quadrangle_4).tags["X1"].GetIds() )
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
