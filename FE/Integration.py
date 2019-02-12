# -*- coding: utf-8 -*-
import numpy as np

from scipy.sparse import coo_matrix
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.Containers.ElementNames as EN

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.Fields.FEField import FEField

#UseCpp = False
UseCpp = True

def Integrate( mesh, wform, constants, fields, dofs,spaces,numbering, tag="3D",integrationRuleName=None,onlyEvaluation=False):
    #conversion of the old api to the new one

    UnkownFields = []
    for i in range(len(dofs)):
        field = FEField()
        field.numbering = numbering[i]
        field.name = dofs[i]
        field.data = None
        field.mesh = mesh
        field.space = spaces[i]
        UnkownFields.append(field)

    return IntegrateGeneral( mesh, wform, constants, list(fields.values()), unkownFields=UnkownFields,testFields=None, tag=tag,integrationRuleName=integrationRuleName,onlyEvaluation=onlyEvaluation)


def IntegrateGeneral( mesh, wform, constants, fields, unkownFields,testFields=None, tag="3D",integrationRuleName=None,onlyEvaluation=False):

    import BasicTools.FE.WeakForm as WeakForm
    if wform is None :
        return

    if not isinstance(wform, WeakForm.PyWeakForm):
        from BasicTools.FE.WeakForm import SymWeakToNumWeak
        wform = SymWeakToNumWeak(wform)


    import time
    st = time.time()
    BaseOutputObject().PrintDebug("Integration ")

    try:
        from BasicTools.FE.WeakFormNumerical import PyWeakForm
        typeCpp = (type(wform) == PyWeakForm)
        import BasicTools.FE.NativeIntegration
    except :
        typeCpp = False
        print("Warning : Using Python Integration (very slow)")

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

    if numberOfVIJ == 0:
        print("Warning!!! System with zero dofs")

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
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOf

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

    constants = {"alpha":1.0}

    pf = FEField("p",mesh,space,numbering)
    pf.Allocate(1)

    fields  = {"p":pf}
    wf = SymWeakToNumWeak(wformflux)

    K,F = Integrate(mesh=mesh,wform=wformflux, tag="ALL", constants=constants, fields=fields, dofs=dofs,spaces=spaces,numbering=numberings)

    if GUI :
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        F.shape = (3,3)
        F = F.T
        mesh.nodeFields["normalflux"] =  F
        OpenInParaView(mesh)

    return "ok"

def CheckIntegrityKF(edim, sdim,testCase):
    from math import sqrt
    import numpy as np
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOfTriangles
    from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOf
    import BasicTools.Containers.ElementNames as EN
    from BasicTools.FE.IntegrationsRules import LagrangeP1

    from BasicTools.FE.MaterialHelp import HookeIso


    if edim == 1:
        nu = 0.25
        #https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch20.d/IFEM.Ch20.pdf
        if sdim == 1:
            E = 1
            A = 1
            points = np.array([[0,],[2,],])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            K = np.matrix([[E],])

            KValidation = (E*A/L)*np.matrix([[1, -1],[-1, 1]])
            permut = None

        elif sdim == 2:
            E = 1000.
            A = 5.
            points = np.array([[0,0],[30,40], ])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            K = A*E*np.matrix([[1 ,1,0],[1,1 ,0],[0 ,0,0]])
            KValidation = np.matrix([[36,48,-36,-48],
                                 [48,64,-48,-64],
                                 [-36,-48,36,48],
                                 [-48,-64,48,64]])
            permut = [0,2,1,3]

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

            permut = [0,2,4,1,3,5]
            KValidation = np.matrix([[40,60,120,-40,-60,-120],
                                 [60,90,180,-60,-90,-180],
                                 [120,180,360,-120,-180,-360],
                                 [-40,-60,-120,40,60,120],
                                 [-60,-90,-180,60,90,180],
                                 [-120,-180,-360,120,180,360]])


        mesh = CreateMeshOf(points,[[0,1],],EN.Bar_2 )
        mesh.GetElementsOfType(EN.Bar_2).tags.CreateTag("1D").SetIds([0])
    elif edim == 2:
        if sdim == 2:
            if testCase[0] =="A" :

                #http://www.unm.edu/~bgreen/ME360/2D%20Triangular%20Elements.pdf

                points = [[3,0],[3,2],[0,0],[0,2] ]
                #,[3,2,1][0,1,2]
                mesh = CreateMeshOfTriangles(points,[[0,1,2],[3,2,1]])
                mesh.GetElementsOfType(EN.Triangle_3).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Triangle_3).GetNumberOfElements() ) )

                E = 3.0/2
                nu = 0.25
                K = HookeIso(E,nu,dim=2)


                permut = [0,4,1,5,3,7,2,6]
#
                KValidation = np.matrix([[0.9833333333333333333,-0.5, -.45, .2, 0 ,0,-0.5333333333333333333, 0.3],
                                [-0.5,1.4,.3,-1.2,0,0,.2,-.2],
                                [-0.45,.3,0.9833333333333333333,0,-0.5333333333333333333,0.2,0,-0.5],
                                [.2,-1.2,0,1.4,.3,-.2,-0.5,0],
                                [0,0,-0.5333333333333333333,.3,0.983333333333333333,-0.5,-0.45,.2],
                                [0,0,0.2,-0.2,-0.5,1.4,.3,-1.2],
                                [-0.5333333333333333333,0.2,0.0,-0.5,-0.45,0.3,0.9833333333333333333,0],
                                [0.3,-0.2,-0.5,0,0.2,-1.2,0,1.4]       ])
                #                     ^^^
                #attention il y a un erreur dans le pdf  (ce termet est possitive dans le pdf)

            elif  testCase[0] == "B" :
                #https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch15.d/IFEM.Ch15.pdf
                #pages 15-11
                points = [[0,0],[3,1],[2,2]]
                mesh = CreateMeshOfTriangles(points,[[0,1,2],])
                mesh.GetElementsOfType(EN.Triangle_3).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Triangle_3).GetNumberOfElements() ) )

                K = np.matrix([[8,2,0],[2,8,0],[0,0,3]])*8
                permut = [0,3,1,4,2,5]

                KValidation = np.matrix([[11,5,-10,-2,-1,-3],
                                 [5,11,2,10,-7,-21],
                                 [-10,2,44,-20,-34,18],
                                 [-2,10,-20,44,22,-54],
                                 [-1,-7,-34,22,35,-15],
                                 [-3,-21,18,-54,-15,75]])
            else:
                raise


        elif sdim == 3:
            if testCase[0] == "A" :
                points = [[0,0,0],[6,-8,5],[6,2,3],[0,5,2] ]
                mesh = CreateMeshOf(points,[[0,1,2,3]],EN.Quadrangle_4)
                E = 3.0/2
                nu = 0.25
                K = HookeIso(E,nu,dim=sdim)
                mesh.GetElementsOfType(EN.Quadrangle_4).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Quadrangle_4).GetNumberOfElements() ) )
                permut = [0,3,6,1,4,7,2,5,8]

            else:
                raise

    elif edim == 3:
        ## https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
        ## page  9-17
        points = [[2.,3.,4],[6,3,2],[2,5,1],[4,3,6]]
        mesh = CreateMeshOf(points,[[0,1,2,3]],EN.Tetrahedron_4)
        E = 480;
        nu = 1./3.

        K = HookeIso(E,nu,dim=sdim)
        mesh.GetElementsOfType(EN.Tetrahedron_4).tags.CreateTag("3D").SetIds(np.arange(mesh.GetElementsOfType(EN.Tetrahedron_4).GetNumberOfElements() ) )

        permut = [0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11]
        KValidation = np.matrix( [[745, 540, 120, -5, 30, 60, -270, -240, 0, -470, -330, -180],
                                 [540, 1720, 270, -120, 520, 210, -120, -1080, -60, -300, -1160, -420],
                                 [120, 270, 565, 0, 150, 175, 0, -120, -270, -120, -300, -470],
                                 [-5, -120, 0, 145, -90, -60, -90, 120, 0, -50, 90, 60],
                                 [30, 520, 150, -90, 220, 90, 60, -360, -60, 0, -380, -180],
                                 [60, 210, 175, -60, 90, 145, 0, -120, -90, 0, -180, -230],
                                 [-270, -120, 0, -90, 60, 0, 180, 0, 0, 180, 60, 0],
                                 [-240, -1080, -120, 120, -360, -120, 0, 720, 0, 120, 720, 240],
                                 [0, -60, -270, 0, -60, -90, 0, 0, 180, 0, 120, 180],
                                 [-470, -300, -120, -50, 0, 0, 180, 120, 0, 340, 180, 120],
                                 [-330, -1160, -300, 90, -380, -180, 60, 720, 120, 180, 820, 360],
                                 [-180, -420, -470, 60, -180, -230, 0, 240, 180, 120, 360, 520]])


    else :
        raise


    #CompureVolume(mesh)
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

    u = GetField("u",sdim)
    ut = GetTestField("u",sdim)

    weak = ToVoigtEpsilon(Strain(u,sdim)).T*K*ToVoigtEpsilon(Strain(ut,sdim))

    wf = SymWeakToNumWeak(weak)

    rwf = wf.GetRightPart(dofs)

    lwf = wf.GetLeftPart(dofs)

    import numpy as np


    constants = {}
    fields  = {}

    import time
    startt = time.time()


    Kinteg,F = Integrate(mesh=mesh,wform=lwf, tag=(str(edim)+"D"), constants=constants, fields=fields, dofs=dofs,spaces=spaces,numbering=numberings)
    stopt = time.time() - startt

    KK = Kinteg.todense()
    if permut is not None:
        KK = KK[:,permut][permut,:]

    error = np.sum(abs(KK-KValidation))/(np.sum(abs(KK))  )

    if error > 1e-14 or error is np.nan:

          print((edim,sdim,testCase))
          print("K Validation")
          print(KValidation)
          print("K Calcul")
          print(KK)

          print("ERROR : "),
          print(error)

          print(KK-KValidation)
          print(numbering)

          raise

    return "ok"

def CompureVolume(mesh):

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo)

    dofs= ["T"]
    spaces = [LagrangeSpaceGeo]
    numberings = [numbering]


    from BasicTools.FE.WeakForm import GetField
    from BasicTools.FE.WeakForm import GetTestField
    from BasicTools.FE.WeakForm import SymWeakToNumWeak


    T = GetField("T",1)
    F = GetField("F",1)
    Tt = GetTestField("T",1)

    ener = T.T*Tt + F.T*Tt


    wf = SymWeakToNumWeak(ener)

    constants = {}
    fields  = {}
    from BasicTools.FE.Fields.FEField import FEField
    f = FEField("f",mesh,LagrangeSpaceGeo,numbering)
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
    global UseCpp
    saveOldeStateOfUseCpp = UseCpp
    UseCpp = False
    if CheckIntegrityNormalFlux(GUI).lower() != "ok":
        return "Not ok in the Normal Calculation"
    UseCpp = True
    if CheckIntegrityNormalFlux(GUI).lower() != "ok":
        return "Not ok in the Normal Calculation"

    print("Normal Calculation OK")

    problems = [ (1,1,"A bars 1D"),
                 (1,2,"A bars 2D"),
                 (1,3,"A bars 3D"),
                 (2,2,"A Triangles in 2D"),
                 (2,2,"B Triangle in 2D"),
                 #(2,3,"A Triangle in 3D"),
                 (3,3,"A tet in 3D"),

                ]
    import time
    startt = time.time()

    for ed,sd,m in problems:
        print("*-"*80)
        print("*-"*80)
        print((ed,sd,m))

        UseCpp = False
        print(" --- python  integration --")
        if CheckIntegrityKF(edim=ed,sdim = sd,testCase=m).lower() !=  "ok":
            return "not ok python "
        else :
            print(" --- python  integration -- OK ")
        print(" --- cpp integration --")
        UseCpp = True
        if CheckIntegrityKF(edim=ed,sdim = sd,testCase=m).lower() !=  "ok":
            return "not ok cpp"
        else:
            print(" --- cpp integration -- OK")

    stopt = time.time() - startt
    print("Total time : ")
    print(stopt)
    print("ALL ok")

    UseCpp = saveOldeStateOfUseCpp
    return "ok"

if __name__ == '__main__':
    print("Start")# pragma: no cover
    print(CheckIntegrity(True))# pragma: no cover
    print("Stop")# pragma: no cover
