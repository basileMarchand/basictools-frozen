# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np
import BasicTools.Containers.ElementNames as EN


def GenerateSpaceForIntegrationPointInterpolation(integrationRule):
    res={}

    from BasicTools.FE.Spaces.SymSpace import SymSpaceBase
    class IntegrationPointSpace(SymSpaceBase):
        def __init__(self,integrationRules, name):
            super(IntegrationPointSpace,self).__init__()

            self.geoSupport = EN.geoSupport[name]
            from sympy.matrices import Matrix
            from sympy import DiracDelta
            integrationPoints  = integrationRules[name][0]
            coord  = [self.xi, self.eta, self.phi]

            self.symN = Matrix([ np.prod([DiracDelta(c-x) for x,c in zip(y,coord)]) for y in integrationPoints])
            self.posN = np.array(integrationPoints)
            self.dofAttachments = [("IP",i,None) for i in range(len(integrationPoints)) ]

    for name in integrationRule:
            res[name] = IntegrationPointSpace(integrationRule,name)
    return res


def CheckIntegrity(GUI=False):
    from BasicTools.FE.IntegrationsRules import LagrangeP1

    a = GenerateSpaceForIntegrationPointInterpolation(LagrangeP1)
    #print(a)
    a[EN.Triangle_3].SetIntegrationRule(LagrangeP1[EN.Triangle_3][0],LagrangeP1[EN.Triangle_3][1])
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
