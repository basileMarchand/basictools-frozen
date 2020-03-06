# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
import BasicTools.Containers.ElementNames as EN


def GenerateSpaceForIntegrationPointInterpolation(integrationRule):
    res={}
    names = [EN.Point_1,
             EN.Bar_2,
             EN.Bar_3,
             EN.Triangle_3,
             EN.Triangle_6,
             EN.Tetrahedron_4,
             EN.Tetrahedron_10,
             EN.Quadrangle_4,
             EN.Hexaedron_8]

    from BasicTools.FE.Spaces.SymSpace import SymSpaceBase
    class IntegrationPointSpace(SymSpaceBase):
        def __init__(self,integrationRules, geoSupport):
            super(IntegrationPointSpace,self).__init__()

            self.geoSupport = geoSupport
            from sympy.matrices import Matrix
            from sympy import DiracDelta
            integrationPoints  = integrationRules[geoSupport][0]
            coord  = [self.xi, self.eta, self.phi]

            self.symN = Matrix([ np.prod([DiracDelta(c-x) for x,c in zip(y,coord)]) for y in integrationPoints])
            self.posN = np.array(integrationPoints)
            self.dofAttachments = [("IP",i,None) for i in range(len(integrationPoints)) ]
            self.Create()

    for name in names:
        if EN.geoSupport[name] in integrationRule:
            res[name] = IntegrationPointSpace(integrationRule,EN.geoSupport[name])
        else:
            continue

    return res


def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1

    a =GenerateSpaceForIntegrationPointInterpolation(LagrangeP1)
    #print(a)
    a[EN.Triangle_3].SetIntegrationRule(LagrangeP1[EN.GeoTri][0],LagrangeP1[EN.GeoTri][1])
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
