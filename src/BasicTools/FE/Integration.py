# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from scipy.sparse import coo_matrix
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers.Filters import ElementFilter, NodeFilter

from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo
from BasicTools.FE.Fields.FEField import FEField

#UseCpp = False
UseCpp = True

def Integrate( mesh, wform, constants, fields, dofs,spaces,numbering, integrationRuleName=None,onlyEvaluation=False,elementFilter=None):
    #conversion of the old api to the new one
    print("BasicTools.FE.Integration:Integrate(...) Warging!! This function will be removed in future versions")

    UnkownFields = []
    for i in range(len(dofs)):
        field = FEField()
        field.numbering = numbering[i]
        field.name = dofs[i]
        field.data = None
        field.mesh = mesh
        field.space = spaces[i]
        UnkownFields.append(field)


    return IntegrateGeneral( mesh,
                            wform,
                            constants,
                            list(fields.values()),
                            unkownFields=UnkownFields,
                            testFields=None,
                            integrationRuleName=integrationRuleName,
                            onlyEvaluation=onlyEvaluation,
                            elementFilter=elementFilter)


def IntegrateGeneral( mesh, wform, constants, fields, unkownFields, testFields=None, integrationRuleName=None,onlyEvaluation=False, elementFilter=None,userIntegrator=None, integrationRule=None):

    if elementFilter is None:
        # if no filter the the integral is over the bulk element
        # 3D elements if the mesh is 3D
        # 2D elements if the mesh is 2D
        elementFilter = ElementFilter()
        elementFilter.SetDimensionality(mesh.GetDimensionality())

    elementFilter.mesh = mesh
    import BasicTools.FE.WeakForms.NumericalWeakForm as WeakForm
    if wform is None :
        return

    if userIntegrator is None:
        ttest = [WeakForm.PyWeakForm]
        try:
            import BasicTools.FE.WearForms.NativeNumericalWeakForm as NativeNumericalWeakForm
            ttest.append(NativeNumericalWeakForm.PyWeakForm)
        except :
            pass

        if not isinstance(wform, tuple(ttest) ):
            from BasicTools.FE.WeakForms.NumericalWeakForm import SymWeakToNumWeak
            wform = SymWeakToNumWeak(wform)

        import time
        st = time.time()
        BaseOutputObject().PrintDebug("Integration ")

        try:
            from BasicTools.FE.WeakForms.NativeNumericalWeakForm import PyWeakForm
            typeCpp = (type(wform) == PyWeakForm)
            import BasicTools.FE.Integrators.NativeIntegration
        except Exception as err :
            typeCpp = False
            print("Error loading c++ integration")
            print(str(err))
            print("Warning : Using Python Integration (very slow)")

        if UseCpp and typeCpp:
            import BasicTools.FE.Integrators.NativeIntegration as NI
            integrator = NI.PyMonoElementsIntegralCpp()
        else:
            import BasicTools.FE.Integrators.PythonIntegration as PI
            integrator = PI.MonoElementsIntegral()
    else:
        integrator = userIntegrator

    if integrationRuleName is None:
       if integrationRule is None:
           from BasicTools.FE.IntegrationsRules import LagrangeIsoParam
           integrationRule = LagrangeIsoParam
    else:
       if integrationRule is None:
           from BasicTools.FE.IntegrationsRules import IntegrationRulesAlmanac as IntegrationRulesAlmanac
           integrationRule = IntegrationRulesAlmanac[integrationRuleName]
       else:
           raise(Exception("must give integrationRuleName or integrationRule not both"))

    ## verification of the integration rule for the ip fields:
    from BasicTools.FE.Fields.IPField import IPField
    for f in fields:
        if isinstance(f,IPField):
            if f.rule != integrationRule:
                print("f.rule")
                print(f.rule)
                print("integrationRule")
                print(integrationRule)
                if integrationRuleName is not None:
                    print("integrationRuleName")
                    print(integrationRuleName)
                raise(Exception("Integration rule of field {f.GetName()} not compatible with the integration" ))

    integrator.SetOnlyEvaluation(onlyEvaluation)
    integrator.SetUnkownFields(unkownFields)
    integrator.SetTestFields(testFields)
    integrator.SetExtraFields(fields )
    integrator.SetConstants(constants)
    integrator.SetIntegrationRule(integrationRule)

    numberOfVIJ = integrator.ComputeNumberOfVIJ(mesh,elementFilter)
    if numberOfVIJ == 0:
        nodalFilter = NodeFilter(mesh=elementFilter.mesh,
                                  tags=elementFilter.tags,
                                  zones=elementFilter.zones)
        numberOfVIJ =  len(nodalFilter.GetIdsToTreat())

    if numberOfVIJ == 0:
        print("Warning!!! System with zero dofs")
        raise(Exception("Error!!! System with zero dofs"))

    vK = np.zeros(numberOfVIJ,dtype=np.float)
    iK = np.empty(numberOfVIJ,dtype=np.int)
    jK = np.empty(numberOfVIJ,dtype=np.int)

    F = np.zeros(integrator.GetTotalTestDofs(),dtype=np.float)

    from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
    if not isinstance(mesh, UnstructuredMesh):
        mesh.GetPosOfNodes()

    integrator.PrepareFastIntegration(mesh,wform,vK,iK,jK,0,F)

    tagFound = False
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearElementContainer

    for name,data,idstotreat in elementFilter:

        if data.GetNumberOfElements() == 0:
            continue
        tagFound = True
        if isinstance(data, ConstantRectilinearElementContainer):
            #TODO : Generate the elementary matrix
            #
            nodesPerElement = data.nodesPerElement
            dofpernode = len(unkownFields)

            nt = (nodesPerElement*dofpernode)**2*len(idstotreat)

            offset = 0
            offsets = []
            for uf in unkownFields:
                offsets.append(offset)
                offset += uf.numbering["size"]

            from BasicTools.FE.FETools import GetElementaryMatrixForFormulation
            elementaryMatrix, rhs = GetElementaryMatrixForFormulation(name,wform, unknownNames = [f.name for f in unkownFields], geoFactor= mesh.GetSpacing())
            #
            edofMat = np.concatenate( [f.numbering[name][idstotreat,:]+offset for f,offset in zip(unkownFields,offsets)] ,axis=1)
            from BasicTools.Containers.UnstructuredMeshInspectionTools import GetElementsFractionInside
            if elementFilter.zonesField is None:
                Eeff = np.ones(len(idstotreat))
            else:
                Eeff = GetElementsFractionInside(elementFilter.zoneField,mesh.nodes,name,data,idstotreat)

            minthreshold = 0.9e-3
            bool_Eeff = (Eeff>=minthreshold)
            nEeff = Eeff[bool_Eeff]
            #TODO check numbering "f" or "c" ???
            sM = ((elementaryMatrix.flatten()[np.newaxis]).T * nEeff.ravel()).flatten(order='F')
            sF = ((rhs.flatten()[np.newaxis]).T * nEeff.ravel()).flatten(order='F')
            nnt = nodesPerElement*dofpernode
            one = np.ones((nnt, 1), dtype=np.int_)
            local_iK = np.kron(edofMat[bool_Eeff,:], one).flatten()
            one.shape = (1,nnt)
            local_jK = np.kron(edofMat[bool_Eeff,:], one).flatten()

            totalvijcpt = integrator.GetNumberOfUsedIvij()
            vK[totalvijcpt:totalvijcpt+nt] = sM
            iK[totalvijcpt:totalvijcpt+nt] = local_iK
            jK[totalvijcpt:totalvijcpt+nt] = local_jK
            np.add.at(F, edofMat.flatten(), sF)
            integrator.AddToNumbefOfUsedIvij(nt)

        else:

            ids = np.array(idstotreat,dtype=int,copy=False)
            if len(ids ) == 0:
                continue
            integrator.ActivateElementType(data)
            integrator.Integrate(wform,ids)

    # in case we have zero element treated, we try to do integrate an integration
    # by points
    # for now only work on rhs terms with no external field
    if tagFound == False:
        ids = nodalFilter.GetIdsToTreat()
        if len(ids) == 0:
            raise(Exception("Tag not found to do the integration"))

        from BasicTools.FE.SymWeakForm import testcharacter

        if testFields is None:
           testFields = []
           for f in unkownFields:
              testFields.append(FEField(name=f.name+testcharacter,mesh=f.mesh,space=f.space,numbering=f.numbering,data=f.data) )

        for monom in wform:
            factor = monom.prefactor
            for term in monom:
                if term.internalType == term.EnumNormal :
                    raise(Exception("no normal"))
                elif  term.internalType == term.EnumConstant :
                    raise(Exception("no constant numerical"))
                elif  term.internalType == term.EnumUnknownField :
                    raise(Exception("no right unknown"))
                elif  term.internalType == term.EnumTestField :
                    if term.derDegree == 1:
                        raise(Exception("No derivative"))
                    offset = 0
                    for tf in testFields:
                        if tf.name == term.fieldName:
                            break
                        offset += tf.numbering["size"]
                    #idx = testFields.find(lambda x:x.name == term.fieldName)
                    for x in ids:
                        leftNumbering = tf.numbering["almanac"][("P",x,None)] + offset
                    continue
                else:
                    print("type " + str(term.internalType) + "Not coded yet")
                    raise(Exception("not coded yet"))

            F[leftNumbering] += factor


    if userIntegrator is None:
        numberOfUsedvij = integrator.GetNumberOfUsedIvij()
        data = (vK[0:numberOfUsedvij], (iK[0:numberOfUsedvij],jK[0:numberOfUsedvij]))
        K = coo_matrix(data, shape=(integrator.GetTotalTestDofs(), integrator.GetTotalUnkownDofs())) .tocsr()
        BaseOutputObject().PrintDebug("Integration Done        " +str(time.time()-st))
        return K,F
    else:
        return

def CheckIntegrityNormalFlux(GUI=False):
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf

    points = [[0,0,0],[1,0,1],[0,1,1] ]
    mesh = CreateMeshOf(points,[[0,1,2]],EN.Triangle_3)

    sdim = 3

    space = LagrangeSpaceGeo
    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    numbering = ComputeDofNumbering(mesh,space)
    dofs= ["u_" + str(i) for i in range(sdim)]
    spaces = [space]*sdim
    numberings = [numbering]*sdim

    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    from BasicTools.FE.SymWeakForm import GetNormal
    from BasicTools.FE.WeakForms.NumericalWeakForm import SymWeakToNumWeak


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
    print(mesh)
    ff = ElementFilter(mesh=mesh) # all elements
    K,F = Integrate(mesh=mesh,wform=wformflux, constants=constants, fields=fields,
                    dofs=dofs,spaces=spaces,numbering=numberings,
                    elementFilter=ff )
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
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles, CreateMeshOf
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
            K = np.array([[E],])

            KValidation = (E*A/L)*np.array([[1, -1],[-1, 1]])
            permut = None

        elif sdim == 2:
            E = 1000.
            A = 5.
            points = np.array([[0,0],[30,40], ])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            K = A*E*np.array([[1 ,1,0],[1,1 ,0],[0 ,0,0]])
            KValidation = np.array([[36,48,-36,-48],
                                 [48,64,-48,-64],
                                 [-36,-48,36,48],
                                 [-48,-64,48,64]])
            permut = [0,2,1,3]

        elif sdim == 3:
            points = np.array([[0,0,0],[2.,3.,6.], ])
            L = sqrt(np.sum((points[1,:] - points[0,:])**2))
            E = 343.
            A = 10.
            K = A*E*np.array([[1,1,1,0,0,0],
                               [1,1,1,0,0,0],
                               [1,1,1,0,0,0],
                               [0,0,0,0,0,0],
                               [0,0,0,0,0,0],
                               [0,0,0,0,0,0],])

            permut = [0,2,4,1,3,5]
            KValidation = np.array([[40,60,120,-40,-60,-120],
                                 [60,90,180,-60,-90,-180],
                                 [120,180,360,-120,-180,-360],
                                 [-40,-60,-120,40,60,120],
                                 [-60,-90,-180,60,90,180],
                                 [-120,-180,-360,120,180,360]])


        mesh = CreateMeshOf(points,[[0,1],],EN.Bar_2 )

    elif edim == 2:
        if sdim == 2:
            if testCase[0] =="A" :

                #http://www.unm.edu/~bgreen/ME360/2D%20Triangular%20Elements.pdf

                points = [[3,0],[3,2],[0,0],[0,2] ]
                #,[3,2,1][0,1,2]
                mesh = CreateMeshOfTriangles(points,[[0,1,2],[3,2,1]])


                E = 3.0/2
                nu = 0.25
                K = HookeIso(E,nu,dim=2)


                permut = [0,4,1,5,3,7,2,6]
#
                KValidation = np.array([[0.9833333333333333333,-0.5, -.45, .2, 0 ,0,-0.5333333333333333333, 0.3],
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
                #mesh.GetElementsOfType(EN.Triangle_3).tags.CreateTag("2D").SetIds(np.arange(mesh.GetElementsOfType(EN.Triangle_3).GetNumberOfElements() ) )

                K = np.array([[8,2,0],[2,8,0],[0,0,3]])*8
                permut = [0,3,1,4,2,5]

                KValidation = np.array([[11,5,-10,-2,-1,-3],
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
        #mesh.GetElementsOfType(EN.Tetrahedron_4).tags.CreateTag("3D").SetIds(np.arange(mesh.GetElementsOfType(EN.Tetrahedron_4).GetNumberOfElements() ) )

        permut = [0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11]
        KValidation = np.array( [[745, 540, 120, -5, 30, 60, -270, -240, 0, -470, -330, -180],
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

    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField

    from BasicTools.FE.SymWeakForm import Gradient
    from BasicTools.FE.SymWeakForm import Strain
    from BasicTools.FE.SymWeakForm import ToVoigtEpsilon
    from BasicTools.FE.WeakForms.NumericalWeakForm import SymWeakToNumWeak

    from sympy import pprint
    from sympy import Symbol

    u = GetField("u",sdim)
    ut = GetTestField("u",sdim)

    weak = ToVoigtEpsilon(Strain(u,sdim)).T@K@ToVoigtEpsilon(Strain(ut,sdim))

    wf = SymWeakToNumWeak(weak)

    rwf = wf.GetRightPart(dofs)

    lwf = wf.GetLeftPart(dofs)

    import numpy as np


    constants = {}
    fields  = {}

    import time
    startt = time.time()

    ff = ElementFilter(mesh,tag=(str(edim)+"D"))

    Kinteg,F = Integrate(mesh=mesh,wform=lwf, constants=constants, fields=fields,
                    dofs=dofs,spaces=spaces,numbering=numberings,
                    elementFilter=ff)
    stopt = time.time() - startt

    KK = Kinteg.todense()
    if permut is not None:
        KK = KK[:,permut][permut,:]

    error = np.sum(abs(KK-KValidation))/(np.sum(abs(KValidation))  )

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

def ComputeVolume(mesh):

    from BasicTools.FE.DofNumbering import ComputeDofNumbering
    numbering = ComputeDofNumbering(mesh,LagrangeSpaceGeo)

    dofs= ["T"]
    spaces = [LagrangeSpaceGeo]
    numberings = [numbering]


    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    from BasicTools.FE.WeakForms.NumericalWeakForm import SymWeakToNumWeak


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

def CheckIntegrityIntegrationWithIntegrationPointField(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1
    #from BasicTools.FE.Spaces.FESpaces import LagrangeSpaceGeo

    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube
    mesh = CreateCube([2.,3.,4.],[-1.0,-1.0,-1.0],[2./10, 2./10,2./10])

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, offset, NGauss = PrepareFEComputation(mesh,numberOfComponents=1)
    from BasicTools.FE.Fields.IPField import IPField

    rho_field = IPField("rho",mesh=mesh,rule=LagrangeP1)
    factor = .1
    rho_field.Allocate(factor)

    from BasicTools.FE.Fields.FEField import FEField



    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    T = GetField("T",1)
    rho = GetField("rho",1)
    Tt = GetTestField("T",1)

    wf = T.T*Tt + rho.T*Tt


    constants = {}
    fields  = {}
    fields = [ rho_field]
    unkownFields = [FEField("T",mesh=mesh,space=space,numbering=numberings[0]) ]
    import time
    startt = time.time()
    ff = ElementFilter(mesh,tag="3D")
    K,F = IntegrateGeneral(mesh=mesh,
                    wform=wf,
                    constants=constants,
                    fields=fields,
                    unkownFields=unkownFields,
                    integrationRuleName="LagrangeP1",
                    elementFilter=ff)

    stopt = time.time() - startt
    volk = np.sum(K)
    print("volume (k): " + str(volk))
    volf = np.sum(F)
    print("volume (f): " + str(volf))
    if volk*factor - volf <   volf/100000000:
        return "ok"

    return "KO"



def CheckIntegrity(GUI=False):
    global UseCpp
    saveOldeStateOfUseCpp = UseCpp
    UseCpp = False
    if CheckIntegrityNormalFlux(GUI).lower() != "ok":
        return "Not ok in the Normal Calculation"
    if CheckIntegrityIntegrationWithIntegrationPointField() != "ok":
        return "Not ok in the integration with IPField"

    UseCpp = True
    if CheckIntegrityNormalFlux(GUI).lower() != "ok":
        return "Not ok in the Normal Calculation"
    if CheckIntegrityIntegrationWithIntegrationPointField() != "ok":
        return "Not ok in the integration with IPField"

    print("Normal Calculation OK")


    print("Integration with IPField OK")

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
        print(" --- cpp integration --",(ed,sd,m ))
        UseCpp = True
        if CheckIntegrityKF(edim=ed,sdim = sd,testCase=m).lower() !=  "ok":
            return "not ok cpp"
        else:
            print(" --- cpp integration -- OK")

    stopt = time.time() - startt
    print("Total time : ")
    print(stopt)
    print("ALL ok")
    CheckIntegrityConstantRectilinearIntegration()
    if CheckIntegrityConstantRectilinearIntegration().lower() !=  "ok":
        return "CheckIntegrityConstantRectilinearIntegration not ok"
    else:
        print(" --- CheckIntegrityConstantRectilinearIntegration -- OK")
    UseCpp = saveOldeStateOfUseCpp

    return "ok"

def CheckIntegrityConstantRectilinearIntegration(GUI=False):
    print("Test integrator on Constant rectilinear mesh")
    from BasicTools.Containers.ConstantRectilinearMesh import ConstantRectilinearMesh
    myCRMesh = ConstantRectilinearMesh(3)
    nx = 3; ny = 4; nz =5;

    myCRMesh.SetDimensions([nx,ny,nz]);
    myCRMesh.SetSpacing([1./(nx-1), 1./(ny-1), 1./(nz-1)]);
    myCRMesh.SetOrigin([0, 0, 0]);

    from BasicTools.FE.FETools import PrepareFEComputation
    space, numberings, offset, NGauss = PrepareFEComputation(myCRMesh,numberOfComponents=1)

    from BasicTools.FE.Fields.FEField import FEField
    from BasicTools.FE.SymWeakForm import GetField
    from BasicTools.FE.SymWeakForm import GetTestField
    T = GetField("T",1)
    Tt = GetTestField("T",1)

    wf = T.T*Tt + 1*Tt

    constants = {}
    fields = [ ]
    unkownFields = [FEField("T",mesh=myCRMesh,space=space,numbering=numberings[0]) ]
    import time
    startt = time.time()
    ff = ElementFilter(myCRMesh)
    print(startt)
    K,F = IntegrateGeneral(mesh=myCRMesh,
                    wform=wf,
                    constants=constants,
                    fields=fields,
                    unkownFields=unkownFields,
                    integrationRuleName="LagrangeP1",
                    elementFilter=ff)

    stopt = time.time() - startt
    volk = np.sum(K)
    print("volume (k): " + str(volk))
    volf = np.sum(F)
    print("volume (f): " + str(volf))
    print(stopt)
    if volk*1 - volf <   volf/100000000:
        return "ok"
    print("volk: "+str(volk))
    print("volf: "+str(volf))

    return "Not Ok"
if __name__ == '__main__':
    print("Start")# pragma: no cover
    print(CheckIntegrity(True))# pragma: no cover
    print("Stop")# pragma: no cover
