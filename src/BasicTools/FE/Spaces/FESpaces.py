# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import BasicTools.Containers.ElementNames as EN
import BasicTools.FE.Spaces.PointSpaces as PointSpaces
import BasicTools.FE.Spaces.BarSpaces as BarSpaces
import BasicTools.FE.Spaces.TriSpaces as TriSpaces
import BasicTools.FE.Spaces.TetSpaces as TetSpaces
import BasicTools.FE.Spaces.QuadSpaces as QuadSpaces
import BasicTools.FE.Spaces.HexaSpaces as HexaSpaces
import BasicTools.FE.Spaces.WedgeSpaces as WedgeSpaces

LagrangeSpaceGeo = {}
LagrangeSpaceGeo[EN.Point_1] = PointSpaces.Point_P0_Lagrange()
LagrangeSpaceGeo[EN.Bar_2] = BarSpaces.Bar_P1_Lagrange()
LagrangeSpaceGeo[EN.Bar_3] = BarSpaces.Bar_P2_Lagrange()
LagrangeSpaceGeo[EN.Triangle_3] = TriSpaces.Tri_P1_Lagrange()
LagrangeSpaceGeo[EN.Triangle_6] = TriSpaces.Tri_P2_Lagrange()
LagrangeSpaceGeo[EN.Tetrahedron_4]  = TetSpaces.Tet_P1_Lagrange()
LagrangeSpaceGeo[EN.Tetrahedron_10] = TetSpaces.Tet_P2_Lagrange()
LagrangeSpaceGeo[EN.Quadrangle_4] = QuadSpaces.Quad_P1_Lagrange()
LagrangeSpaceGeo[EN.Quadrangle_8] = QuadSpaces.Quad8_P2_Lagrange()
LagrangeSpaceGeo[EN.Quadrangle_9] = QuadSpaces.Quad_P2_Lagrange()
LagrangeSpaceGeo[EN.Hexaedron_8]  = HexaSpaces.Hexa_P1_Lagrange()
LagrangeSpaceGeo[EN.Hexaedron_20] = HexaSpaces.Hexa20_P2_Lagrange()
LagrangeSpaceGeo[EN.Hexaedron_27] = HexaSpaces.Hexa_P2_Lagrange()
LagrangeSpaceGeo[EN.Wedge_6] = WedgeSpaces.Wedge_P1_Lagrange()


ConstantSpaceGlobal = {}
ConstantSpaceGlobal[EN.Point_1] = PointSpaces.Point_P0_Global()
ConstantSpaceGlobal[EN.Bar_2] = BarSpaces.Bar_P0_Global()
ConstantSpaceGlobal[EN.Bar_3] = BarSpaces.Bar_P0_Global()
ConstantSpaceGlobal[EN.Triangle_3] = TriSpaces.Tri_P0_Global()
ConstantSpaceGlobal[EN.Triangle_6] = TriSpaces.Tri_P0_Global()
ConstantSpaceGlobal[EN.Tetrahedron_4]  = TetSpaces.Tet_P0_Global()
ConstantSpaceGlobal[EN.Tetrahedron_10] = TetSpaces.Tet_P0_Global()
ConstantSpaceGlobal[EN.Quadrangle_4] = QuadSpaces.Quad_P0_Global()
ConstantSpaceGlobal[EN.Quadrangle_8] = QuadSpaces.Quad_P0_Global()
ConstantSpaceGlobal[EN.Quadrangle_9] = QuadSpaces.Quad_P0_Global()
ConstantSpaceGlobal[EN.Hexaedron_8]  = HexaSpaces.Hexa_P0_Global()
ConstantSpaceGlobal[EN.Hexaedron_20] = HexaSpaces.Hexa_P0_Global()
ConstantSpaceGlobal[EN.Hexaedron_27] = HexaSpaces.Hexa_P0_Global()
ConstantSpaceGlobal[EN.Wedge_6] = WedgeSpaces.Wedge_P0_Global()

LagrangeSpaceP0 = {}
LagrangeSpaceP0[EN.Point_1] = PointSpaces.Point_P0_Lagrange()
LagrangeSpaceP0[EN.Bar_2] = BarSpaces.Bar_P0_Lagrange()
LagrangeSpaceP0[EN.Bar_3] = BarSpaces.Bar_P0_Lagrange()
LagrangeSpaceP0[EN.Triangle_3] = TriSpaces.Tri_P0_Lagrange()
LagrangeSpaceP0[EN.Triangle_6] = TriSpaces.Tri_P0_Lagrange()
LagrangeSpaceP0[EN.Tetrahedron_4]  = TetSpaces.Tet_P0_Lagrange()
LagrangeSpaceP0[EN.Tetrahedron_10] = TetSpaces.Tet_P0_Lagrange()
LagrangeSpaceP0[EN.Quadrangle_4] = QuadSpaces.Quad_P0_Lagrange()
LagrangeSpaceP0[EN.Quadrangle_8] = QuadSpaces.Quad_P0_Lagrange()
LagrangeSpaceP0[EN.Quadrangle_9] = QuadSpaces.Quad_P0_Lagrange()
LagrangeSpaceP0[EN.Hexaedron_8]  = HexaSpaces.Hexa_P0_Lagrange()
LagrangeSpaceP0[EN.Hexaedron_20] = HexaSpaces.Hexa_P0_Lagrange()
LagrangeSpaceP0[EN.Hexaedron_27] = HexaSpaces.Hexa_P0_Lagrange()
LagrangeSpaceP0[EN.Wedge_6] = WedgeSpaces.Wedge_P0_Lagrange()

LagrangeSpaceP1 = {}
LagrangeSpaceP1[EN.Point_1] = PointSpaces.Point_P0_Lagrange()
LagrangeSpaceP1[EN.Bar_2] = BarSpaces.Bar_P1_Lagrange()
LagrangeSpaceP1[EN.Bar_3] = BarSpaces.Bar_P1_Lagrange()
LagrangeSpaceP1[EN.Triangle_3] = TriSpaces.Tri_P1_Lagrange()
LagrangeSpaceP1[EN.Triangle_6] = TriSpaces.Tri_P1_Lagrange()
LagrangeSpaceP1[EN.Tetrahedron_4]  = TetSpaces.Tet_P1_Lagrange()
LagrangeSpaceP1[EN.Tetrahedron_10] = TetSpaces.Tet_P1_Lagrange()
LagrangeSpaceP1[EN.Quadrangle_4] = QuadSpaces.Quad_P1_Lagrange()
LagrangeSpaceP1[EN.Quadrangle_8] = QuadSpaces.Quad_P1_Lagrange()
LagrangeSpaceP1[EN.Quadrangle_9] = QuadSpaces.Quad_P1_Lagrange()
LagrangeSpaceP1[EN.Hexaedron_8]  = HexaSpaces.Hexa_P1_Lagrange()
LagrangeSpaceP1[EN.Hexaedron_20] = HexaSpaces.Hexa_P1_Lagrange()
LagrangeSpaceP1[EN.Hexaedron_27] = HexaSpaces.Hexa_P1_Lagrange()
LagrangeSpaceP1[EN.Wedge_6]      = WedgeSpaces.Wedge_P1_Lagrange()

LagrangeSpaceP2 = {}
LagrangeSpaceP2[EN.Point_1] = PointSpaces.Point_P0_Lagrange()
LagrangeSpaceP2[EN.Bar_2] = BarSpaces.Bar_P2_Lagrange()
LagrangeSpaceP2[EN.Bar_3] = BarSpaces.Bar_P2_Lagrange()
LagrangeSpaceP2[EN.Triangle_3] = TriSpaces.Tri_P2_Lagrange()
LagrangeSpaceP2[EN.Triangle_6] = TriSpaces.Tri_P2_Lagrange()
LagrangeSpaceP2[EN.Tetrahedron_4]  = TetSpaces.Tet_P2_Lagrange()
LagrangeSpaceP2[EN.Tetrahedron_10] = TetSpaces.Tet_P2_Lagrange()
LagrangeSpaceP2[EN.Quadrangle_4] = QuadSpaces.Quad_P2_Lagrange()
LagrangeSpaceP2[EN.Quadrangle_8] = QuadSpaces.Quad_P2_Lagrange()
LagrangeSpaceP2[EN.Quadrangle_9] = QuadSpaces.Quad_P2_Lagrange()
LagrangeSpaceP2[EN.Hexaedron_8]  = HexaSpaces.Hexa_P2_Lagrange()
LagrangeSpaceP2[EN.Hexaedron_20] = HexaSpaces.Hexa_P2_Lagrange()
LagrangeSpaceP2[EN.Hexaedron_27] = HexaSpaces.Hexa_P2_Lagrange()
LagrangeSpaceP2[EN.Wedge_6] = WedgeSpaces.Wedge_P2_Lagrange()

def CheckIntegrity(GUI=False):
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
