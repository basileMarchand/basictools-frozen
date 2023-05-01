# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Gmsh tools
"""

import BasicTools.Containers.ElementNames as EN

gmshName = {}
gmshName[EN.Bar_2]          = '1'
gmshName[EN.Triangle_3]     = '2'
gmshName[EN.Quadrangle_4]   = '3'
gmshName[EN.Tetrahedron_4]  = '4'
gmshName[EN.Hexaedron_8]    = '5'
gmshName[EN.Wedge_6]        = '6'
gmshName[EN.Pyramid_5]      = '7'
gmshName[EN.Bar_3]          = '8'
gmshName[EN.Triangle_6]     = '9'
gmshName[EN.Quadrangle_9]   = '10'
gmshName[EN.Tetrahedron_10] = '11'
gmshName[EN.Point_1]        = '15'
gmshName[EN.Quadrangle_8]   = '16'
gmshName[EN.Hexaedron_20]   = '17'
gmshName[EN.Wedge_15]       = '18'
gmshName[EN.Pyramid_13]     = '19'

gmshNumber = {}
gmshNumber['1'] = EN.Bar_2
gmshNumber['2'] = EN.Triangle_3
gmshNumber['3'] = EN.Quadrangle_4
gmshNumber['4'] = EN.Tetrahedron_4
gmshNumber['5'] = EN.Hexaedron_8
gmshNumber['6'] = EN.Wedge_6
gmshNumber['7'] = EN.Pyramid_5
gmshNumber['8'] = EN.Bar_3
gmshNumber['9'] = EN.Triangle_6
gmshNumber['10'] = EN.Quadrangle_9
gmshNumber['11'] = EN.Tetrahedron_10
gmshNumber['15'] = EN.Point_1
gmshNumber['16'] = EN.Quadrangle_8
gmshNumber['17'] = EN.Hexaedron_20
gmshNumber['18'] = EN.Wedge_15
gmshNumber['19'] = EN.Pyramid_13

PermutationGmshToBasicTools = {}
PermutationGmshToBasicTools[EN.Tetrahedron_10] = [2, 0, 1, 3, 6, 4, 5, 8, 7, 9]
PermutationGmshToBasicTools[EN.Hexaedron_20]   = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15]
PermutationGmshToBasicTools[EN.Wedge_15]       = [0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11]
PermutationGmshToBasicTools[EN.Pyramid_13]     = [0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12]

def CheckIntegrity():
    return "OK"
