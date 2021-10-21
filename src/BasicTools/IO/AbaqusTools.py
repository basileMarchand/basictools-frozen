# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import BasicTools.Containers.ElementNames as EN


InpNameToBasicTools = {}

InpNameToBasicTools["S3"] = EN.Triangle_3
InpNameToBasicTools['CONN3D2'] = EN.Bar_2
InpNameToBasicTools['T2D2'] = EN.Bar_2
InpNameToBasicTools['CPS4R'] = EN.Quadrangle_4
InpNameToBasicTools["C3D4"] = EN.Tetrahedron_4
InpNameToBasicTools["C3D8"] = EN.Hexaedron_8
InpNameToBasicTools["C3D8R"] = EN.Hexaedron_8
InpNameToBasicTools["C3D8I"] = EN.Hexaedron_8
InpNameToBasicTools["C3D10"] = EN.Tetrahedron_10
InpNameToBasicTools["C3D10M"] = EN.Tetrahedron_10
InpNameToBasicTools["C3D20"] = EN.Hexaedron_20
InpNameToBasicTools["C3D6"] = EN.Wedge_6
InpNameToBasicTools["T3D2"] = EN.Bar_2
InpNameToBasicTools["CPS3"] = EN.Triangle_3


permutation = {}
#permutation[ EN.Tetrahedron_4] = [0, 1, 3, 2]

BasicToolsToInpName = {}

for k,v in InpNameToBasicTools.items():
    if k in BasicToolsToInpName:
        continue
    BasicToolsToInpName[v] = k


