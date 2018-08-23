# -*- coding: utf-8 -*-
import numpy as np

from scipy.sparse import coo_matrix
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.FE.ElementNames as EN

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.Fields.NodalField import NodalField

#UseCpp = False
UseCpp = True

def Integrate( mesh, wform, constants, fields, dofs,spaces,numbering, tag="3D",integrationRuleName=None,onlyEvaluation=False):
    #conversion of the old api to the new one

    UnkownFields = []
    for i in range(len(dofs)):
        field = NodalField()
        field.numbering = numbering[i]
        field.name = dofs[i]
        field.data = None
        field.mesh = mesh
        field.space = spaces[i]
        UnkownFields.append(field)

    return IntegrateGeneral( mesh, wform, constants, list(fields.values()), unkownFields=UnkownFields,testFields=None, tag=tag,integrationRuleName=integrationRuleName,onlyEvaluation=onlyEvaluation)


def IntegrateGeneral( mesh, wform, constants, fields, unkownFields,testFields=None, tag="3D",integrationRuleName=None,onlyEvaluation=False):

    import time
    st = time.time()
    BaseOutputObject().PrintDebug("Integration ")

    try:
        from BasicTools.FE.WeakFormNumerical import PyWeakForm
        typeCpp = (type(wform) == PyWeakForm)
    except :
        typeCpp = False

    if UseCpp and typeCpp:
        import BasicTools.FE.NativeIntegration as NI
        integrator = NI.PyMonoElementsIntegralCpp()
    else:
        import BasicTools.FE.PythonIntegration as PI
        integrator = PI.MonoElementsIntegral()

    integrator.SetOnlyEvaluation(onlyEvaluation)
    integrator.SetUnkownFields(unkownFields)
    integrator.SetTestFields(testFields)
    integrator.SetExtraFields(fields )
    integrator.SetConstants(constants)
    integrator.SetIntegrationRule(integrationRuleName)

    numberOfVIJ = integrator.ComputeNumberOfVIJ(mesh,tag)

    vK = np.empty(numberOfVIJ,dtype=np.float)
    iK = np.empty(numberOfVIJ,dtype=np.int)
    jK = np.empty(numberOfVIJ,dtype=np.int)

    F = np.zeros(integrator.GetTotalTestDofs(),dtype=np.float)

    integrator.PrepareFastIntegration(mesh,wform,vK,iK,jK,0,F)

    for name, data in mesh.elements.items():

        if data.GetNumberOfElements() == 0:
            continue

        if tag in data.tags:
            idstotreat = data.tags[tag].GetIds()
            if len(idstotreat) == 0:
                    continue
        elif tag == "ALL":
                idstotreat = np.arange(data.GetNumberOfElements(),dtype=int)
        else:
            continue

        integrator.ActivateElementType(data)
        integrator.Integrate(wform,idstotreat)

    numberOfUsedvij = integrator.GetNumberOfUsedIvij()
    data = (vK[0:numberOfUsedvij], (iK[0:numberOfUsedvij],jK[0:numberOfUsedvij]))
    K = coo_matrix(data, shape=(integrator.GetTotalTestDofs(), integrator.GetTotalUnkownDofs())) .tocsr()
    BaseOutputObject().PrintDebug("Integration Done        " +str(time.time()-st))
    return K,F


def CheckIntegrityNormalFlux(GUI=False):
    from BasicTools.FE.UnstructuredMeshTools import CreateMeshOf

    #points = [[0,0,0],[6,-8,5],[6,2,3],[0,5,2] ]
    #mesh = CreateMeshOf(points,[[0,1,2,3]],EN.Quadrangle_4)
    #mesh.GetElementsOfType(EN.Quadrangle_4).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Quadrangle_4).GetNumberOfElements() ) )

    points = [[0,0,0],[1,0,1],[0,1,1] ]
    mesh = CreateMeshOf(points,[[0,1,2]],EN.Triangle_3)
    mesh.GetElementsOfType(EN.Triangle_3).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Triangle_3).GetNumberOfElements() ) )

    sdim = 3

    space = LagrangeSpaceGeo
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    numbering = ComputeDofNumbering(mesh,space)
    dofs= ["u_" + str(i) for i in range(sdim)]
    spaces = [space]*sdim
    numberings = [numbering]*sdim

    from BasicTools.FE.WeakForm import GetField
    from BasicTools.FE.WeakForm import GetTestField
    from BasicTools.FE.WeakForm import GetNormal
    from BasicTools.FE.WeakForm import SymWeakToNumWeak


    u = GetField("u",sdim)
    ut = GetTestField("u",sdim)
    p = GetField("p",1)

    Normal = GetNormal(3)

    wformflux = p*Normal.T*ut

    print(numberings)
    print(wformflux)
    constants = {"alpha":1.0}
    from BasicTools.FE.Fields.NodalField import NodalField as Field
    pf = Field("p",mesh,space,numbering)
    pf.Allocate(1)

    fields  = {"p":pf}
    wf = SymWeakToNumWeak(wformflux)
    print(wf)

    K,F = Integrate(mesh=mesh,wform=wf, tag="ALL", constants=constants, fields=fields, dofs=dofs,spaces=spaces,numbering=numberings)
    print(list(F))

    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        F.shape = (3,3)
        F = F.T
        mesh.nodeFields["normalflux"] =  F
        OpenInParaView(mesh)
    print(np.sum(np.abs(F) ) )
    return "ok"

def CheckIntegrityKF(case, sdim,nmesh=1):
    from math import sqrt
    import numpy as np
    from BasicTools.FE.UnstructuredMeshTools import CreateMeshOfTriangles
    from BasicTools.FE.UnstructuredMeshTools import CreateMeshOf
    import BasicTools.FE.ElementNames as EN
    from BasicTools.FE.IntegrationsRules import LagrangeP1

    from BasicTools.FE.MaterialHelp import HookeIso


    if case == 1:


        nu = 0.25
        #https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch20.d/IFEM.Ch20.pdf
        if sdim == 1:
            E = 1
            A = 1
            points = np.array([[0],[1],])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            K = np.matrix([[E],])
        elif sdim == 2:
            E = 1000.
            A = 5.
            points = np.array([[0,0],[30,40], ])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            K = (E)*np.matrix([[1,0,0],[0,1,0],[0,0,0]])
            K = A*E*np.matrix([[1 ,1,0],[1,1 ,0],[0 ,0,0]])
#            print(K)

        elif sdim == 3:
            points = np.array([[0,0,0],[2.,3.,6.], ])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            E = 343.
            A = 10.
            K = A*E*np.matrix([[1,1,1,0,0,0],
                               [1,1,1,0,0,0],
                               [1,1,1,0,0,0],
                               [0,0,0,0,0,0],
                               [0,0,0,0,0,0],
                               [0,0,0,0,0,0],])


        #HookeIso(E,nu,dim=sdim)

        mesh = CreateMeshOf(points,[[0,1],],EN.Bar_2 )
        mesh.GetElementsOfType(EN.Bar_2).tags.CreateTag("1D").SetIds([0])
    elif case == 2:
        if sdim == 2:
        #    points = [[0,0],[3,1],[2,2]]
        #    mesh = CreateMeshOfTriangles(points,[[0,1,2],])
        #elif sdim == 3:,3
            if nmesh == 1:
                #http://www.unm.edu/~bgreen/ME360/2D%20Triangular%20Elements.pdf
                points = [[3,0],[3,2],[0,0],[0,2] ]
                #,[3,2,1][0,1,2]
                mesh = CreateMeshOfTriangles(points,[[0,1,2],[3,2,1]])
                E = 3.0/2
                nu = 0.25
                K = HookeIso(E,nu,dim=2)



            elif nmesh == 2:
                points = [[0,0],[3,1],[2,2]]
                mesh = CreateMeshOfTriangles(points,[[0,1,2],])
                K = np.matrix([[8,2,0],[2,8,0],[0,0,3]])*8

        elif sdim ==3:
            if nmesh == 1:
                points = [[0,0,0],[6,-8,5],[6,2,3],[0,5,2] ]
                mesh = CreateMeshOf(points,[[0,1,2,3]],EN.Quadrangle_4)
                E = 3.0/2
                nu = 0.25
                K = HookeIso(E,nu,dim=sdim)
                mesh.GetElementsOfType(EN.Quadrangle_4).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Quadrangle_4).GetNumberOfElements() ) )
            elif nmesh == 2:
                points = [[3,0,0],[3,2,0],[0,0,1],[0,2,1] ]
                #,[3,2,1][0,1,2]
                mesh = CreateMeshOfTriangles(points,[[0,1,2],[3,2,1]])
                E = 3.0/2
                nu = 0.25
                K = HookeIso(E,nu,dim=2)
        mesh.GetElementsOfType(EN.Triangle_3).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Triangle_3).GetNumberOfElements() ) )

    elif case ==3:
        #from BasicTools.FE.UnstructuredMeshTools import CreateCube
        #mesh = CreateCube(origin=[0,0,0])
        data=u"""
***geometry
**node
8 3
1 0.0000000000000000e+00          0.0000000000000000e+00          0.0000000000000000e+00
2 1.0000000000000000e+00          0.0000000000000000e+00          0.0000000000000000e+00
3 0.0000000000000000e+00          1.0000000000000000e+00          0.0000000000000000e+00
4 1.0000000000000000e+00          1.0000000000000000e+00          0.0000000000000000e+00
5 0.0000000000000000e+00          0.0000000000000000e+00          1.0000000000000000e+00
6 1.0000000000000000e+00          0.0000000000000000e+00          1.0000000000000000e+00
7 0.0000000000000000e+00          1.0000000000000000e+00          1.0000000000000000e+00
8 1.0000000000000000e+00          1.0000000000000000e+00          1.0000000000000000e+00
**element
1
1 c3d8 1 2 4 3 5 6 8 7
**elset 3D
1
***return

"""
        from BasicTools.IO.GeofReader import ReadGeof
        mesh = ReadGeof(string=data)

#        print(mesh)
#        print(mesh.nodes)
        E = 2.E9
        nu = 0.3
#        print(mesh.GetElementsOfType(EN.Hexaedron_8).connectivity)
        K = HookeIso(E,nu,dim=sdim)

    else :
        raise()


    #CompureVolume(mesh)
#    print("----------------------- K -------------------------")
#    print(K)
    space = LagrangeSpaceGeo
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    numbering = ComputeDofNumbering(mesh,space)

    if sdim == 1:
        dofs = ["u"]
    else:
        dofs= ["u_" + str(i) for i in range(sdim)]

    spaces = [space]*sdim
    numberings = [numbering]*sdim

    from BasicTools.FE.WeakForm import GetField
    from BasicTools.FE.WeakForm import GetTestField

    from BasicTools.FE.WeakForm import Gradient
    from BasicTools.FE.WeakForm import Strain
    from BasicTools.FE.WeakForm import ToVoigtEpsilon
    from BasicTools.FE.WeakForm import SymWeakToNumWeak

    from sympy import pprint
    from sympy import Symbol
    #init_session()
#    print(space)

    u = GetField("u",sdim)
    ut = GetTestField("u",sdim)

#    print(u)
#    print(u.shape)
#
#    print("-----------------")
#    pprint(u)
#    pprint(Gradient(u,sdim=sdim))
#    pprint(Strain(u,sdim=sdim))
#    pprint(ToVoigtEpsilon(Strain(u,sdim)))



#    print(K)

    #weak = 3.*ut[1]
    #weak = 3*u[0].diff(Symbol("x"))* ut[1].diff(Symbol("x"))
    weak = ToVoigtEpsilon(Strain(u,sdim)).T*K*ToVoigtEpsilon(Strain(ut,sdim))
    #from sympy.matrices import Matrix
    #+ Matrix([alpha*ut[0]])#f.T*ut*alpha

    #ener = Matrix([alpha*Tt[0]])#f.T*ut*alpha
    #+f.diff(PhysicalSpace[0]).T*ut
    #ener = f[0]*u.T*ut
    pprint(weak)
    print(weak)



    wf = SymWeakToNumWeak(weak)
    print(wf)

    rwf = wf.GetRightPart(dofs)
#    print([str(r) for r in rwf.form])

    lwf = wf.GetLeftPart(dofs)
#    print([str(r) for r in lwf.form])
#
#    print("lwf")
#    print(lwf)
#    print("rwf")
#    print(rwf)

    import numpy as np

    #alpha = GetConstant("alpha")

    constants = {}
    fields  = {}
    #f0 = Field("f_0",mesh,space,numbering)
    #f0.Allocate(1)
    #fields["f0"] = f0

#    print(mesh)
    import time
    startt = time.time()
    print("lwf")
    print(lwf)
    print(mesh)
    Kinteg,F = Integrate(mesh=mesh,wform=lwf, tag=(str(case)+"D"), constants=constants, fields=fields, dofs=dofs,spaces=spaces,numbering=numberings)
    stopt = time.time() - startt
    print(F)

#    print("!-!-!-!-!-!-!-!-!-!-!")


#    if sdim ==2 and case ==2:
#        name = EN.Triangle_3
#        domain = mesh.GetElementsOfType(name)
#        geoSpace = LagrangeSpaceGeo[name]
#        #construction of the Bmatrix for the triangle
#        p,w =  LagrangeP1[EN.geoSupport[name]]
#        geoSpace.SetIntegrationRule(p,w)
#        xcoor = mesh.nodes[domain.connectivity[0],:]
#        Jinv, Jdet = geoSpace.GetJackAndDetI(0,xcoor)
#        Nfder =  np.array([geoSpace.valdNdxi[0].T, geoSpace.valdNdeta[0].T])
#        #print(Nfder)
#        BxByBz = np.dot(Jinv,Nfder)/Jdet
#        #print(BxByBz)
#        B = np.matrix([[BxByBz[0][0], 0 ,BxByBz[0][1],0,BxByBz[0][2],0 ],
#               [0, BxByBz[1][0], 0 ,BxByBz[1][1],0,BxByBz[1][2]],
#               #[BxByBz[1][0]/2,BxByBz[0][0]/2, BxByBz[1][1] /2,BxByBz[0][1]/2,BxByBz[1][2]/2,BxByBz[0][2]/2 ],
#               [BxByBz[1][0],BxByBz[0][0], BxByBz[1][1],BxByBz[0][1],BxByBz[1][2],BxByBz[0][2]],
#               ])
#        t = 1
#        A = 3.
#        print("+++++****************++++++++++")
#        print(t)
#        print("A")
#        print(A)
#        print(B)
#        print(K)
#        print("B.T*K*B")
#        print(B.T*K*B)
#        Kref = t*A*B.T*K*B
        #print(Kref)
        #error = np.sum(abs(Kinteg.todense()[:,[0,4,1,5,2,6,3,7]][[0,4,1,5,2,6,3,7],:] -Kref))
        #print(error)
        #if error  >1e-15 :
        #    pass
        #    raise ("error")




#    print("---------------------------------time : "),
#    print(stopt)
    #print(F)
    #print("KintegOriginal")
    #print(Kinteg.todense())

    if sdim ==1 and case ==1:
        KK = Kinteg.todense()
        KValidation = (E*A/L)*np.matrix([[1, -1],[-1, 1]])
    elif sdim ==2 and case ==1:
        KK = Kinteg.todense()[:,[0,2,1,3]][[0,2,1,3],:]
        KValidation = np.matrix([[36,48,-36,-48],
                                 [48,64,-48,-64],
                                 [-36,-48,36,48],
                                 [-48,-64,48,64]])
    elif sdim ==3 and case ==1:
        KK = Kinteg.todense()[:,[0,2,4,1,3,5]][[0,2,4,1,3,5],:]
        KValidation = np.matrix([[40,60,120,-40,-60,-120],
                                 [60,90,180,-60,-90,-180],
                                 [120,180,360,-120,-180,-360],
                                 [-40,-60,-120,40,60,120],
                                 [-60,-90,-180,60,90,180],
                                 [-120,-180,-360,120,180,360]])

    elif sdim ==1 and case ==2:
        KK = Kinteg.todense()
        KValidation = np.matrix([[0.5, -0.5],[-0.5, 0.5]])
    elif sdim ==3 and case == 2 and nmesh == 1 :
        KK = Kinteg.todense()
        KValidation = np.matrix([[0.5, -0.5],[-0.5, 0.5]])

    elif sdim ==2 and case == 2 and nmesh == 1 :
        if mesh.GetNumberOfElements() ==1:
            ind =[0,3,1,4,2,5]
            name = EN.Triangle_3
            domain = mesh.GetElementsOfType(name)
            geoSpace = LagrangeSpaceGeo[name]
            #construction of the Bmatrix for the triangle
            p,w =  LagrangeP1[EN.geoSupport[name]]
            geoSpace.SetIntegrationRule(p,w)
            xcoor = mesh.nodes[domain.connectivity[0],:]
            Jack, Jdet, Jinv = geoSpace.GetJackAndDetI(0,xcoor)
            Nfder =  np.array([geoSpace.valdNdxi[0].T, geoSpace.valdNdeta[0].T])
            BxByBz = np.dot(Jinv,Nfder)/Jdet

            B = np.matrix([[BxByBz[0][0], 0 ,BxByBz[0][1],0,BxByBz[0][2],0 ],
               [0, BxByBz[1][0], 0 ,BxByBz[1][1],0,BxByBz[1][2]],
               #[BxByBz[1][0]/2,BxByBz[0][0]/2, BxByBz[1][1] /2,BxByBz[0][1]/2,BxByBz[1][2]/2,BxByBz[0][2]/2 ],
               [BxByBz[1][0],BxByBz[0][0], BxByBz[1][1],BxByBz[0][1],BxByBz[1][2],BxByBz[0][2]],
               ])
            t = 1
            A = 3.
            KValidation = t*A*B.T*K*B
#            print("Jinv")
#            print(Jinv)
#            print("Jdet")
#            print(Jdet)
#            print("B")
#            print(B)*6
        else:
            #ind =[0,4,2,6,3,7,1,5]
            ind =[0,4,1,5,3,7,2,6]

            KValidation = np.matrix([[0.9833333333333333333,-0.5, -.45, .2, 0 ,0,-0.5333333333333333333, 0.3],
                                [-0.5,1.4,.3,-1.2,0,0,.2,-.2],
                                [-0.45,.3,0.9833333333333333333,0,-0.5333333333333333333,0.2,0,-0.5],
                                [.2,-1.2,0,1.4,.3,-.2,-0.5,0],
                                [0,0,-0.5333333333333333333,.3,0.983333333333333333,-0.5,-0.45,.2],
                                [0,0,0.2,-0.2,-0.5,1.4,.3,-1.2],
                                [-0.5333333333333333333,0.2,0.0,-0.5,-0.45,0.3,0.9833333333333333333,0],
                                [0.3,-0.2,-0.5,0,0.2,-1.2,0,1.4]       ])
        #                             ^^^
        #attention il y a un erreur dans le pdf  (ce termet est possitive dans le pdf)
        KK = Kinteg.todense()[:,ind][ind,:]
        #print(KK)
        np.set_printoptions(threshold=np.nan)


    elif sdim ==2 and case == 2 and nmesh == 2 :
        KK = Kinteg.todense()[:,[0,3,1,4,2,5]][[0,3,1,4,2,5],:]
        KValidation = np.matrix([[11,5,-10,-2,-1,-3],
                                 [5,11,2,10,-7,-21],
                                 [-10,2,44,-20,-34,18],
                                 [-2,10,-20,44,22,-54],
                                 [-1,-7,-34,22,35,-15],
                                 [-3,-21,18,-54,-15,75]])

    elif sdim == 3 and case == 2 and nmesh == 1 :
        KK = Kinteg.todense()[:,[0,3,6,1,4,7,2,5,8]][[0,3,6,1,4,7,2,5,8],:]

    elif sdim == 3 and case == 3  :
        libtozsetnum = [None]*(3*8)
        libtozsetnum[slice(0,None,3)] =  np.arange(0,8)
        libtozsetnum[slice(1,None,3)] =  np.arange(0,8)+8
        libtozsetnum[slice(2,None,3)] =  np.arange(0,8)+2*8

        KK = Kinteg.todense()[:,libtozsetnum][libtozsetnum,:]
        from scipy.io import mmread
        print(KK)
        KValidation = mmread("/scratch/fcasenave/partage/For_Felipe/0_Zmatrix_1_0").todense()
        #np.set_printoptions(threshold=np.nan)
    else:
        raise

#
    print("K Validation")
    print(KValidation)
    print("K Calcul")
    print(KK)
    print((case,sdim,nmesh))
    error = np.sum(abs(KK-KValidation))/(np.sum(abs(KK))  )


    print("ERROR : "),
    print(error)

    if error > 1e-14:
          print(KK-KValidation)
          print(numbering)
          raise
    raise
    return "ok"




def CompureVolume(mesh):

    import BasicTools.FE.Elements.SymElement as SE
    import BasicTools.FE.ElementNames as EN

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo)

    dofs= ["T"]
    spaces = [LagrangeSpaceGeo]
    numberings = [numbering]


    from BasicTools.FE.WeakForm import GetField
    from BasicTools.FE.WeakForm import GetTestField
    from BasicTools.FE.WeakForm import SymWeakToNumWeak

    from sympy import pprint
    #init_session()
#    print(LagrangeSpaceGeo)

    T = GetField("T",1)
    F = GetField("F",1)
    Tt = GetTestField("T",1)

    ener = T.T*Tt + F.T*Tt


    wf = SymWeakToNumWeak(ener)

    constants = {}
    fields  = {}
    from BasicTools.FE.Fields.NodalField import NodalField as Field
    f = Field("f",mesh,LagrangeSpaceGeo,numbering)
    f.Allocate(1)
    fields["F"] = f

    import time
    startt = time.time()
    K,F = Integrate(mesh=mesh,wform=wf, tag="ALL", constants=constants, fields=fields, dofs=dofs,spaces=spaces,numbering=numberings)
    stopt = time.time() - startt
    volk = np.sum(K)
    print("volume (k): " + str(volk))
    volf = np.sum(F)
    print("volume (f): " + str(volf))
    if abs(volk - volf) > 1e-15 :
        print(volk-volf)
        raise


def CheckIntegrity(GUI=False):
    if CheckIntegrityNormalFlux(GUI).lower() != "ok":
        return "Not ok in the Normal Calculation"

    print("Normal Calculation OK")
    problems = [#(3,3,1), # hexa in 3D
                #(2,3,1), # quad in 3D
                #(2,3,2), # tri in 3D
                (2,2,2), # triangles in 2D http://www.unm.edu/~bgreen/ME360/2D%20Triangular%20Elements.pdf
                (2,2,1), # triangles in 2D
                (1,3,1), # bar , espace 3D
                (1,2,1), # bar , espace 2D
                (1,1,1), # bar , espace 1d
                ]
    import time
    startt = time.time()

    for c,d,m in problems:
        print("*-"*80)
        print("*-"*80)
        print((c,d,m))
        if CheckIntegrityKF(case=c,sdim = d,nmesh=m).lower() !=  "ok":
            return "not ok"

    stopt = time.time() - startt
    print("Total time : ")
    print(stopt)
    print("ALL ok")

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(False))# pragma: no cover
