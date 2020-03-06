# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers import Filters

def Hash(dofAttachment,lconn,name=None,i=None,discontinuous=False):
    """
    Function to create a unique tuple for a given shape function (dof) this
    tuple must be independent of the element that use this shape function ( an
    exception is the discontinues galerkin approximation).

    """
    T,n,h = dofAttachment
    if T == "P":
        """P is for point"""
        T = "P"
        n = lconn[n]
    elif T == "C":
        """C is for cell"""
        T = name
        if h is not None:
            T2, n2, h2 = h
            h2 = Hash(h,lconn=lconn)
            return (T,tuple(np.sort(lconn)),h2)
        return (T,tuple(np.sort(lconn)),n)

    elif T == "F" :
        """is for face  (face for a 3D element, edge for a 2D element """
        edge = EN.faces[name][n]
        T = edge[0]
        n = tuple(np.sort(lconn[edge[1]]))
        h = 0

    elif T == "F2":
        """is for face second level (edge for a 3D element, point for a 2D element """
        edge = EN.faces2[name][n]
        T = edge[0]
        n = tuple(np.sort(lconn[edge[1]]))
        h = 0

    elif T == "G":
        """G is for global """
        return (T,0,h)

    elif T == "IP":
        return (T,tuple(lconn),i)

    else:
        raise(Exception(" type of dof unknown : " + str(T) ) )

    T += ("_"+ str(i)) if discontinuous else ""
    return (T,n,h)


def ComputeDofNumbering(mesh,Space,dofs=None,sign=1,fromConnectivity=False,elementFilter=None):
    """
    Function to compute a unique numbering of dofs. The user must provide:

    * mesh : (UnstructuredMesh), discretisation
    * space : (from FESpaces.py for example) approximation space
    * dofs : () previous dofNumbering computation
    * tag : (str) tag name to treat
    * sign : (default 1) user can give (-1) to use negative numbering
    * fromConnectivity : (default False) for fast computation of isoparametric
      numbering for all element (space, dofs,tag,sign are ignored). All Elements
      are treated.
    """

    if fromConnectivity:
        if dofs is not None or sign != 1 or elementFilter is not None:
            raise(Exception("cant take dofs tag sign or elementFilter different from the default values"))
        dofs = {}
        dofs["size"] = mesh.GetNumberOfNodes()
        dofs["dirichelet"] = -1
        dofs["doftopointLeft"] = slice(mesh.GetNumberOfNodes())
        dofs["doftopointRight"] = dofs["doftopointLeft"]
        for name,data in mesh.elements.items():
            dofs[name] = data.connectivity

        cpt = dofs["size"]
        almanac = {}
        for i in range(cpt):
            almanac[('P', i, None)] = i
        dofs["almanac"] = almanac
    else:

        if elementFilter is None:
            elementFilter = Filters.ElementFilter(mesh)

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

        for name,data,fil in elementFilter:
            if len(fil) == 0:
                continue

            sp = Space[name]

            if name in dofs:
                dof = dofs[name]
            else:
                dof = np.zeros((data.GetNumberOfElements(),sp.GetNumberOfShapeFunctions()), dtype=np.int_) -1


            numberOfShapeFunctions = sp.GetNumberOfShapeFunctions()
            for i in fil:
                conn = data.connectivity[i,:]
                for j in range(numberOfShapeFunctions):
                    h = Hash(sp.dofAttachments[j],lconn=conn,name=name,i=i,discontinuous = sp.classification.discontinuous)

                    if h in almanac:
                        d = almanac[h]
                    else:
                        d = cpt
                        almanac[h] = d
                        cpt += 1*sign
                    dof[i,j] = d

            dofs[name] = dof

    if sign ==1:
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

        extractorLeftSide = np.resize(extractorLeftSide, (tmpcpt,))
        extractorRightSide = np.resize(extractorRightSide, (tmpcpt,))
        dofs["doftopointLeft"] =   extractorLeftSide
        dofs["doftopointRight"] =  extractorRightSide

        extractorLeftSide = np.empty(cpt,dtype=np.int)
        extractorRightSide = np.empty(cpt,dtype=np.int)

        tmpcpt = 0
        # if k[0] is the elementname then k[1] is the connecivity
        # we generate the same almanac with the number of each element
        elemDic = {}
        for name,data in mesh.elements.items():
            #num = elemDic.get(name,np.zeros(data.GetNumberOfElements()) -1)
            elemDic[name] = {}
            elemDic2 = elemDic[name]

            for i in range(data.GetNumberOfElements()):
                elemDic2[tuple(np.sort(data.connectivity[i,:]))] = i

        for k,v in almanac.items():
            #if not k[0] in {'P',"F","F2","G"} :
            #we need the global number of the element (not the local to the element container)
            if k[0] in elemDic.keys():
                localdic = elemDic[k[0]]
                if k[1] in localdic.keys():
                    extractorLeftSide[tmpcpt] = mesh.elements[k[0]].globaloffset + localdic[k[1]]
                    extractorRightSide[tmpcpt] = v
                    tmpcpt += 1

        extractorLeftSide = np.resize(extractorLeftSide, (tmpcpt,))
        extractorRightSide = np.resize(extractorRightSide, (tmpcpt,))

        dofs["doftocellLeft"] =   extractorLeftSide
        dofs["doftocellRight"] =  extractorRightSide

    else:
        dofs["dirichelet"] = cpt

    return dofs

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube

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


    import BasicTools.Containers.ElementNames as EN
    print(res2.GetElementsOfType(EN.Quadrangle_4).tags["X1"].GetIds() )

    from BasicTools.FE.Spaces.FESpaces import ConstantSpaceGlobal

    t_numbering = ComputeDofNumbering(res2,ConstantSpaceGlobal)
    print(t_numbering)

    from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceP2

    t_numbering = ComputeDofNumbering(res2,LagrangeSpaceP2)
    print(t_numbering)
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
