# -*- coding: utf-8 -*-
import numpy as np

from scipy.sparse import coo_matrix

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Linalg.MatOperations import deleterowcol
from BasicTools.FE.Fields.ConstantField import ConstantField
from BasicTools.FE.Fields.FEField import FEField

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.DofNumbering import ComputeDofNumbering
from BasicTools.FE.Fields.IntegrationPointField import IntegrationPointField
from BasicTools.FE.IntegrationsRules import LagrangeP1


import sys


class UMatFemProblem(BaseOutputObject):
    """
    Class to solve a finit element problem using the local behaviour law subroutine
    Zmat. For the moment only 3D solid problem are treated.

    mesh    : the discretisation of the domain
    dtime   : time variation of the current step
    temp    : current temperature
    dtemp   : variation of the temperature
    dstrain : variation of the straint field in the currect step
    tag     : the name of the domain to be used in the calculation
    u       : displacement field
    du      : variation of the displacement field
    pyumat  : the wraper to the behaviour law

    """

    cpt = 0

    def __init__(self):
        super(UMatFemProblem,self).__init__()
        self.mesh = None

        self.dtime = 1

        self.temp = ConstantField("temp",20)
        self.dtemp = ConstantField("dtemp",0)

        self.dstrain = None

        self.tag = "3D"

        self.u = [ FEField() for i in range(3)]
        self.du = [ FEField() for i in range(3)]

        if sys.version_info[0] == 2:
            import pyumat
            #return __import__('pyumat')
        else:
            import py3umat as pyumat
            #return __import__('py3umat')
        self.pyumat = pyumat

    def SetMesh(self,mesh):
        """
        Function to set the discratisation.

        This function will:
            - define the spaces (LagrangeSpaceGeo)
            - compute the numbering (il uses the point ids for the dof numbering)
            - allocate the fields for the strain
            - allocate the fields for the displacement the the displacment increment
        """
        self.mesh = mesh
        dim = mesh.GetDimensionality()
        self.spaces = [LagrangeSpaceGeo]*dim

        numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo,fromConnectivity =True)
        self.numberings = [numbering]*dim

        dstrain = [IntegrationPointField("dstrain"+i) for i in ["11", "22", "33", "12", "31", "23"]  ]
        for i in range(6):
            dstrain[i].Allocate(self.mesh,LagrangeP1,tag=self.tag)
        self.dstrain = dstrain

        for i in range(3):
            self.u[i].mesh = self.mesh
            self.u[i].space = self.spaces[i]
            self.u[i].numbering = self.numberings[i]
            self.u[i].Allocate(0)
            self.du[i].mesh = self.mesh
            self.du[i].space = self.spaces[i]
            self.du[i].numbering = self.numberings[i]
            self.du[i].Allocate(0)


    def SetMaterial(self,filename=None,string=None):
        """
        Function to initialise the behavoiur law, the user can give a filename,
        or a string containing the behavoiur law parameters.
        """

        if filename is not None:
            self.matfile = filename
        elif string is not None:
            from BasicTools.Helpers.Tests import TestTempDir
            import os
            fullfilename = TestTempDir.GetTempPath() + str(UMatFemProblem.cpt) + ".m"


            f = open(fullfilename,"w")
            f.write(string)

            print("material file in " + str(fullfilename) )
            f.close()

            self.matfile = fullfilename.split(os.sep)[-1]
        else:
            raise(Exception("need a file or a string"))

        code="""
import sys
if sys.version_info[0] == 2:
    import pyumat
else:
    import py3umat as pyumat

import numpy as np

cmname = '"""+self.matfile+"""'

nstatv = 83
ndi = 3
nshr = 3
ntens = ndi+nshr

stress = np.zeros(6)

statev = np.zeros(nstatv)

ddsdde = np.zeros((6,6),dtype=np.float)

sse = 0.
spd = 0.
scd = 0.

rpl = 0.
ddsddt = np.zeros(ntens)
drplde = np.zeros(ntens)
drpldt = 0.

stran = np.zeros(6)
dstran = np.zeros(6)
r2 = 2


timesim = np.array([0.,0.1])
dtime = 0.1# self.data['dt']
temperature = 20.
dtemp = 0.
predef = np.zeros(1)
dpred = np.zeros(1)


nprops = 1
props = np.zeros(nprops)
coords = np.zeros(3, dtype= float)
drot  =  np.zeros((3,3), dtype= float)
pnewdt = 1.
celent = 1.
dfgrd0 = np.zeros((3,3), dtype = float)
dfgrd1 = np.zeros((3,3), dtype = float)
noel = -1
npt = -1
kslay = 1
kspt = 1
kstep = np.array([1,1,0,0], dtype=int)
kinc  = 1

ddsddeNew = pyumat.umat(stress=stress,statev=statev,ddsdde=ddsdde,sse=sse,spd=spd,scd=scd,rpl=rpl,ddsddt=ddsddt,drplde=drplde,drpldt=drpldt,stran=stran,dstran=dstran,time=timesim,dtime=dtime,
                        temp=temperature,dtemp=dtemp,predef=predef,dpred=dpred,cmname=cmname,ndi=ndi,nshr=nshr,ntens=ntens,nstatv=nstatv,props=props,nprops=nprops,coords=coords,drot=drot,pnewdt=pnewdt,celent=celent,dfgrd0=dfgrd0,
        dfgrd1=dfgrd1,noel=noel,npt=npt,kslay=kslay,kspt=kspt,kstep=kstep,kinc=kinc)
"""
        f = open("materialtest.py","w")
        f.write(code)
        f.close()

        import sys
        print('Running zset to get info ')
        if sys.version_info[0] == 2:
            # hack for python2
            from subprocess import check_output
            out = check_output([sys.executable, "materialtest.py"])
        else:
            import subprocess
            out = subprocess.run([sys.executable, "materialtest.py"], stdout=subprocess.PIPE).stdout.decode("utf-8")

        print('Running zset to get info DONE')
        outlines = out.split(u"\n")


        names = ['Flux', 'Grad','var_int','var_aux','Extra Zmat']

        def parser(fname, obj):
            cont = 1
            line = 0
            while line < len(outlines):
                if fname in outlines[line]:
                    line +=1
                    while any(word in outlines[line].strip() for word in names) == False and len(outlines[line].strip())>0  and outlines[line][0] != "=":
                        if outlines[line][0] != " ":
                            break
                        fnames = outlines[line].strip().split()
                        for i in range(len(fnames)):
                          if '--' in fnames[i]:
                            cont = 0
                          if cont == 1:
                            obj.append(fnames[i].split("(")[0])
                        line +=1
                else:
                    line +=1

        self.flux = []
        self.grad = []
        self.var_int = []
        self.var_aux = []
        self.var_extra = []
        parser("Flux", self.flux)
        parser("Grad", self.grad)
        parser("var_int", self.var_int)
        parser("var_aux", self.var_aux)
        parser("Extra Zmat", self.var_extra)
        print("Flux :\n",self.flux)
        print("Grad :\n",self.grad)
        print("var_int :\n",self.var_int)
        print("var_aux :\n",self.var_aux)
        print("Extra Zmat :\n",self.var_extra)



    def AllocateMaterialData(self):
        """
        This function must be called to allocate the fields for the integration
        points data

        Attention, the order of the data is not the same as in the umat !!!!
        """

        def create(obj):
          from BasicTools.FE.Fields.IntegrationPointField import IntegrationPointField as IPF
          for i in range(len(obj)):
            name = obj[i]
            obj[i] = IPF(name)
            obj[i].Allocate(self.mesh,LagrangeP1,tag=self.tag)


        #stran = np.array([eto11[count], eto22[count], eto33[count], eto12[count], eto31[count], eto23[count]])
        create(self.flux)
        #reorder flux
        self.flux = [ self.flux[i] for i in [0,1,2,3,5,4]]

        create(self.grad)
        #reorder grad
        self.grad = [ self.grad[i] for i in [0,1,2,3,5,4]]

        create(self.var_int)
        create(self.var_aux)
        create(self.var_extra)



    def GetDofsForElementTag(self,tagname,dof=None):
        """
        Get the dofs numbers of all the nodes present in a tag

        tagname : (string) name of the tag
        dof : [0,1,2,None] Ux,Uy,Uz or all the dofs for the current points

        """
        nbdofs = self.mesh.GetNumberOfNodes()
        dofs = []
        for name,data in self.mesh.elements.items():

            if tagname in data.tags:
                ids = data.connectivity[data.tags[tagname].GetIds(),:].flatten()
                if dof == 0:
                    dofs.extend(ids)
                elif dof == 1:
                    dofs.extend(ids+nbdofs)
                elif dof == 2:
                    dofs.extend(ids+2*nbdofs)
                elif dof is None:
                    dofs.extend(np.concatenate((ids,ids+nbdofs,ids+2*nbdofs)))
                else:
                    raise
        return np.unique(dofs)

    def CallUmatOnIntegrationPointI(self,elemtype,el,ip,updateInternals=True):
        """
        Function to call the umat routine on a specifique integration point of
        one element:

        elemtype : (string) the element type
        et       :  (int)   the local id of the element
        ip       : (int)    the id of the integration point
        updateInternals : (bool) if we update the internal variable with the output
            of the routine
        """
        cmname = self.matfile

        nstatv = len(self.var_int) + len(self.var_aux) +  len(self.var_extra)

        ndi = self.mesh.GetDimensionality()# 3
        if ndi == 3:
            nshr = 3
        else:
            nshr = 1

        ntens = ndi+nshr

        stress = np.zeros(6)

        statev = np.zeros(nstatv)  # <-----
        cpt =0
        for var in range(len(self.var_int)):
            statev[cpt] = self.var_int[var].GetValueAtIP(elemtype,el,ip)
            cpt += 1

        for var in range(len(self.var_aux)):
            statev[cpt] = self.var_aux[var].GetValueAtIP(elemtype,el,ip)
            cpt += 1

        for var in range(len(self.var_extra)):
            statev[cpt] = self.var_extra[var].GetValueAtIP(elemtype,el,ip)
            cpt += 1

        ddsdde = np.zeros((6,6),dtype=np.float)

        sse = 0.
        spd = 0.
        scd = 0.

        rpl = 0.
        ddsddt = np.zeros(ntens)
        drplde = np.zeros(ntens)
        drpldt = 0.

        stran = np.zeros(6)  # <-----
        dstran = np.zeros(6) # <-----
        temp = self.dstrain
        for i in range(6):
            stran[i]  = self.grad[i].GetValueAtIP(elemtype,el,ip)
            dstran[i] = temp[i].GetValueAtIP(elemtype,el,ip)

        r2 = 2

        timesim = np.array([0.,0.1])
        dtime = self.dtime
        temperature = self.temp.GetValueAtIP(elemtype,el,ip)
        dtemp = self.dtemp.GetValueAtIP(elemtype,el,ip)
        predef = np.zeros(1)
        dpred = np.zeros(1)


        nprops = 1
        props = np.zeros(nprops)
        coords = np.zeros(3, dtype= float)
        drot  =  np.zeros((3,3), dtype= float)
        pnewdt = 1.
        celent = 1.
        dfgrd0 = np.zeros((3,3), dtype = float)
        dfgrd1 = np.zeros((3,3), dtype = float)
        noel = el
        npt = ip
        kslay = 1
        kspt = 1
        kstep = np.array([1,1,0,0], dtype=int)
        kinc  = 1

        #print(stran)
        ddsddeNew = self.pyumat.umat(stress=stress,statev=statev,ddsdde=ddsdde,sse=sse,spd=spd,scd=scd,rpl=rpl,ddsddt=ddsddt,drplde=drplde,drpldt=drpldt,stran=stran,dstran=dstran,time=timesim,dtime=dtime,
        temp=temperature,dtemp=dtemp,predef=predef,dpred=dpred,cmname=cmname,ndi=ndi,nshr=nshr,ntens=ntens,nstatv=nstatv,props=props,nprops=nprops,coords=coords,drot=drot,pnewdt=pnewdt,celent=celent,dfgrd0=dfgrd0,
        dfgrd1=dfgrd1,noel=noel,npt=npt,kslay=kslay,kspt=kspt,kstep=kstep,kinc=kinc)


        if updateInternals:

            for var in range(len(self.flux)):
                self.flux[var].SetValueAtIP(elemtype,el,ip,stress[var] )
                #print(stress[var])

            cpt =0
            for var in range(len(self.var_int)):
                self.var_int[var].SetValueAtIP(elemtype,el,ip,statev[cpt] )
                cpt +=1

            for var in range(len(self.var_aux)):
                self.var_aux[var].SetValueAtIP(elemtype,el,ip,statev[cpt] )
                cpt +=1

            for var in range(len(self.var_extra)):
                self.var_extra[var].SetValueAtIP(elemtype,el,ip,statev[cpt] )
                cpt +=1
        return ddsddeNew, stress



    def ViewInParaView(self):
        """
        Just a helper function to view the mesh
        """
        from BasicTools.Actions.OpenInParaView import OpenInParaView
#        F = f.view()
#        F.shape = (3,mesh.GetNumberOfNodes())
#        F = F.T
#        mesh.nodeFields["normalflux"] =  F[numbering["permutation"],:]
#        fixedValues[dofs] = res
#        fixedValues.shape = (3,mesh.GetNumberOfNodes())
#        fixedValues= fixedValues.T
#        mesh.nodeFields["sol"] =  fixedValues[numbering["permutation"],:]
        OpenInParaView(self.mesh)

    def GenerateKF(self,computeK=True,computeF=True,updateInternals=True):
        """
        General function to compute lhs, rhs and to update the internal variables
        of the material

        computeK : (bool) to compute the lhs
        computeF : (bool) to compute the rhs
        updateInternals : (bool) to update the internal variables of the material


        """
        #total number of dofs and offset calculation
        offset = []
        totaldofs = 0
        for n in self.numberings:
            offset.append(totaldofs)
            totaldofs += n["size"]

        if computeF :
            F = np.zeros(totaldofs,dtype=np.float)
        else:
            F = None

        #matrix allocation
        if computeK:
            numberOfVIJ = 0
            for name,data in self.mesh.elements.items():
                tag = self.tag
                if tag in data.tags:
                    numberOfUsedElements = len(data.tags[tag])
                    ss = 0
                    for space in self.spaces:
                        ss +=space[name].GetNumberOfShapeFunctions()
                    numberOfVIJ += numberOfUsedElements*(ss**2)

            vK = np.empty(numberOfVIJ,dtype=np.float)
            iK = np.empty(numberOfVIJ,dtype=int)
            jK = np.empty(numberOfVIJ,dtype=int)

            vij = [vK,iK,jK,0]
        else:
            vij = None

        for name, data in self.mesh.elements.items():
            elspace = [b[name] for b in self.spaces]
            elnumbering = [b[name] for b in self.numberings]
            if self.tag in data.tags:
               idstotreat = data.tags[tag].GetIds()
               elF = self.ElementsIntegral(data,vij,elspace,elnumbering,offset,totaldofs,idstotreat,computeK=computeK,F=F,updateInternals=updateInternals )
            if elF is not None:
                F += F
        if computeK:
            K = coo_matrix((vij[0][0:vij[3]], (vij[1][0:vij[3]],vij[2][0:vij[3]])), shape=(totaldofs, totaldofs))
        else :
            K = None

        return (K,F)

    def ElementsIntegral(self,domain,vij,elspace,elnumbering,offset,totaldofs,idstotreat,computeK,F,updateInternals ):
        """
        Internal function to compute the tanget matrix using the umat

        domain : (ElementsContainer) a storage of element of the same type
        vij    : (list(np.array(float)),np.array(int)),np.array(int)))
                 list containing the 3 vectors (values,Is,Js), for the construction
                 of the sparse matrix
        elspace : (SpaceBase) approximation space for the current type of elements
        elnumbering : list[np.arrays(ints)] numbering of the dofs
        offset : list(ints) the offset of the numbering with respect to the precedent field
        totaldofs : (int) total number of dofs (the size of the tanget matrix)
        idstotreat : list(ints) the element  to be treated
        computeK : (bool) if true the filling of vij is done
        F : [None, np.array(float)] if not none then F is fillet with the rhs constribution
        updateInternals : (bool) to update the internal variable of the class with
            the calculated data from the umat

        """
        from BasicTools.FE.IntegrationsRules import LagrangeP1
        import BasicTools.Containers.ElementNames as EN

        p,w =  LagrangeP1[EN.geoSupport[domain.elementType]]
        elspace[0].SetIntegrationRule(p,w)
        space = elspace[0]
        nbsf = space.GetNumberOfShapeFunctions()
        B = np.zeros((6,nbsf*3), dtype=np.float)
#
#                        [[dN_x[0], 0      , dN_x[1], 0      , dN_x[2], 0      ],
#                      [0      , dN_y[0], 0      , dN_y[1], 0      , dN_y[2]],
#                      [dN_y[0], dN_x[0], dN_y[1], dN_x[1], dN_y[2], dN_x[2]]]);
        for n in idstotreat:
            ev = []
            ei = []
            ej = []
            xcoor = self.mesh.nodes[domain.connectivity[n],:]

            leftNumbering = np.concatenate((elnumbering[0][n,:]+ offset[0],
                                           elnumbering[1][n,:]+ offset[1],
                                           elnumbering[2][n,:]+ offset[2]))

            rightNumbering = leftNumbering
            #print(list(offset))
            #print(list(rightNumbering))

            for ip in range(len(w)):
                 Jack, Jdet, Jinv = space.GetJackAndDetI(ip,xcoor)
                 BxByBzI = Jinv(space.valdphidxi[ip])
                 for i in range(nbsf):
                     B[0,i] =   BxByBzI[0,i]
                     B[1,i+nbsf] = BxByBzI[1,i]
                     B[2,i+2*nbsf] = BxByBzI[2,i]

                     B[3,i] =   BxByBzI[1,i]
                     B[3,i+nbsf] =   BxByBzI[0,i]

                     B[4,i] =   BxByBzI[2,i]
                     B[4,i+2*nbsf] =   BxByBzI[0,i]

                     B[5,i+nbsf] =   BxByBzI[2,i]
                     B[5,i+2*nbsf] =   BxByBzI[1,i]

                 u = np.concatenate((self.u[0].data[elnumbering[0][n,:]] ,
                                     self.u[1].data[elnumbering[1][n,:]],
                                     self.u[2].data[elnumbering[2][n,:]]))

                 epsilon = np.dot(B,u)
                 for i in range(6):
                     self.grad[i].SetValueAtIP(domain.elementType,n,ip,epsilon[i])

                 ddsddeNew,sigma = self.CallUmatOnIntegrationPointI(domain.elementType,n,ip,updateInternals=updateInternals)

                 if computeK:
                     vals = w[ip]*Jdet*B.T.dot(ddsddeNew.dot(B))
                     ev.extend(vals.ravel())
                     ones = np.ones( len(rightNumbering) ,dtype=int)
                     for i in leftNumbering:
                         ej.extend( i * ones)
                         ei.extend(rightNumbering.ravel())

                 if F is not None:
                     vals = w[ip]*Jdet*B.T.dot(sigma)
                     F[leftNumbering] = vals

            if len (ev) and computeK:

                data = coo_matrix((ev, (ei,ej)), shape=( totaldofs,totaldofs)).tocsr().tocoo()
                start = vij[3]
                stop = start+len(data.data)

                vij[0][start:stop] = data.data
                vij[1][start:stop] = data.row
                vij[2][start:stop] = data.col
                vij[3] += len(data.data)
        return

    def SolveUsing(self,k,f,fixedValues=None,freeMask=None):
        """
        Simple function to solve the problem ku=f qith the fixed values.

        k : (matrix(floats)) linear operator
        f : (vector(floats)) rhs
        fixedValues : (vector(floats)) values of Dirichlet boundary conditions
        freeMask : (vector(bool)) true to impose Dirichlet value to the courrent dof

        """
        if  fixedValues is not None:
            [K, rhsfixed] = deleterowcol(k.tocsr(), np.logical_not(freeMask), np.logical_not(freeMask), fixedValues)
            rhs = f[freeMask]-rhsfixed[freeMask]

            from BasicTools.Linalg.LinearSolver import LinearProblem

            prob = LinearProblem()
            prob.SetAlgo("Direct")
            prob.SetOp(K.tocsc())
            res = prob.Solve(rhs)
            fixedValues[freeMask] = res


def CheckIntegrity(GUI=False):

    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("ZSET_NO_FAIL"): return "ok"

    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.IO.PathControler import TemporalChdir
    with TemporalChdir(TestTempDir.GetTempPath()):

        from BasicTools.Containers.UnstructuredMeshTools import CreateCube
        nx= 10
        ny= 10
        nz= 10
        mesh = CreateCube([nx,ny,nz],[0,0,0],[1./(nx-1) ,1./(ny-1) ,1./(nz-1)])

        from BasicTools.FE.Behaviours.ZMatLinearMaterials import GetMaterial

        props = {"YOUNG":4.e5,"POISSON":0.3}
        matstring = GetMaterial('Elastic',props)
        matstring = """
***behavior gen_evp
 **elasticity isotropic
   young 2.1e5
   poisson 0.3
# **potential gen_evp ep
#  *criterion mises
#  *flow plasticity
#  *isotropic nonlinear
#   R0 180.
#   Q  360.
#   b   30.
***return"""
        prob = UMatFemProblem()
        print("init mesh")
        prob.SetMesh(mesh)
        print("init material")
        prob.SetMaterial(string=matstring)
        prob.AllocateMaterialData()



        import BasicTools.Containers.ElementNames as EN
        prob.CallUmatOnIntegrationPointI(EN.Hexaedron_8,0,0)

        nx = 10
        stran = [None]*nx
        stress = [None]*nx
        for i in range(nx):
            dx = nx*i/100000000.
            prob.dstrain[0].SetValueAtIP(EN.Hexaedron_8,0,0,dx)
            prob.dstrain[1].SetValueAtIP(EN.Hexaedron_8,0,0,-dx/3*0)
            prob.dstrain[2].SetValueAtIP(EN.Hexaedron_8,0,0,-dx/3*0)
            prob.CallUmatOnIntegrationPointI(EN.Hexaedron_8,0,0)
            stran[i] = dx
            stress[i] = prob.flux[0].GetValueAtIP(EN.Hexaedron_8,0,0)
        import time
        stime = time.time()
        k, f = prob.GenerateKF()
        print(k)
        print(f)
        print("time" + str(stime-time.time()) )

        fixedValues = np.zeros(f.shape, dtype=np.float)
        freeMask = np.ones(f.shape, dtype=np.bool)

        block0 = prob.GetDofsForElementTag("X0")
        print("---")
        print(block0 )
        freeMask[block0] = False
        block1 = prob.GetDofsForElementTag("X1",dof=0)
        freeMask[block1] = False
        fixedValues[block1] = 0.1

        prob.SolveUsing(k,f, fixedValues,freeMask)

        print("stran  =", stran)
        print("stress =", stress)
        if GUI :
            from BasicTools.Actions.OpenInParaView import OpenInParaView
            #F = f.view()
            #F.shape = (3,mesh.GetNumberOfNodes())
            #F = F.T
            #mesh.nodeFields["normalflux"] =  F[numbering["permutation"],:]
            #fixedValues[dofs] = res
            mesh.nodeFields["solx"] = fixedValues[0:mesh.GetNumberOfNodes()]
            fixedValues.shape = (3,mesh.GetNumberOfNodes())
            fixedValues= fixedValues.T
            mesh.nodeFields["sol"] =  fixedValues;

            OpenInParaView(mesh)

            import matplotlib.pyplot as plt
            plt.plot(stran,stress)
            plt.xlabel('strain11')
            plt.ylabel('stress11')
            plt.show()

        return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))#pragma: no cover
