# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Ansys tools
"""

import BasicTools.Containers.ElementNames as EN

PermutationAnsysToBasicTools = {}
PermutationAnsysToBasicTools[ EN.Tetrahedron_10 ] = [2, 0, 1, 3, 6, 4, 5, 9, 7, 8]

def CheckIntegrity():
    return "OK"
