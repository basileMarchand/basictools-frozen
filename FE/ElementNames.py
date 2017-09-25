# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"


numberOfNodes = {}

# permutation of index to make a valid element again after a mirror operation
mirrorPermutation = {}
dimension = {}
linear = {}
faces = {}

#0d
Point_1  = 'point1'
numberOfNodes[Point_1] = 1;
mirrorPermutation[Point_1] = [0]
dimension[Point_1] = 0
linear[Point_1] = True
#1d
#linear
Bar_2 = 'bar2'
numberOfNodes[Bar_2] = 2;
mirrorPermutation[Bar_2] = [1,0]
dimension[Bar_2] = 1
linear[Bar_2] = True
#quadratic
Bar_3 = 'bar3'
numberOfNodes[Bar_3] = 3;
mirrorPermutation[Bar_3] = [1,0,2]
dimension[Bar_3] = 1
linear[Bar_3] = False
#2d
#linear
Triangle_3 = 'tri3'
numberOfNodes[Triangle_3] = 3;
mirrorPermutation[Triangle_3] = [0,2,1]
dimension[Triangle_3] = 2
linear[Triangle_3] = True

Quadrangle_4  = 'quad4'
numberOfNodes[Quadrangle_4] = 4;
dimension[Quadrangle_4] = 2
linear[Quadrangle_4] = True

#quadratic
Triangle_6 = 'tri6'
numberOfNodes[Triangle_6] = 6;
dimension[Triangle_6] = 2
linear[Triangle_6] = True


Quadrangle_8  = 'quad8'
numberOfNodes[Quadrangle_8] = 8;
dimension[Quadrangle_8] = 2
Quadrangle_9  = 'quad9'
numberOfNodes[Quadrangle_9] = 9;
dimension[Quadrangle_9] = 2
#3d
#linear
Tetrahedron_4 = 'tet4'
numberOfNodes[Tetrahedron_4] = 4;
mirrorPermutation[Tetrahedron_4] = [0,2,1,3]
dimension[Tetrahedron_4] = 3
linear[Tetrahedron_4] = True
faces[Tetrahedron_4] = [(Triangle_3,[0, 2, 1]),
                        (Triangle_3,[0, 1, 3]),
                        (Triangle_3,[1, 2, 3]),
                        (Triangle_3,[2, 0, 3]),
     ]

Pyramid_5  = 'pyr5'
numberOfNodes[Pyramid_5] = 5;
dimension[Pyramid_5] = 3
linear[Pyramid_5] = True

Wedge_6 = 'wed6'
numberOfNodes[Wedge_6] = 6;
dimension[Wedge_6] = 3
linear[Wedge_6] = True

Hexaedron_8 = 'hex8'
numberOfNodes[Hexaedron_8] = 8;
dimension[Hexaedron_8] = 3
linear[Hexaedron_8] = True

#quadratic
Tetrahedron_10 = 'tet10'
numberOfNodes[Tetrahedron_10] = 10;
dimension[Tetrahedron_10] = 3
linear[Tetrahedron_10] = False

Pyramid_13  = 'pyr13'
numberOfNodes[Pyramid_13] = 13;
dimension[Pyramid_13] = 3
linear[Pyramid_13] = False

Wedge_15 = 'hex15'
numberOfNodes[Wedge_15] = 15;
dimension[Wedge_15] = 3
linear[Wedge_15] = False

Wedge_18 = 'hex18'
numberOfNodes[Wedge_18] = 18;
dimension[Wedge_18] = 3
linear[Wedge_18] = False


Hexaedron_20 = 'hex20'
numberOfNodes[Hexaedron_20] = 20;
dimension[Hexaedron_20] = 3
linear[Hexaedron_20] = False

Hexaedron_27 = 'hex27'
numberOfNodes[Hexaedron_27] = 27;
dimension[Hexaedron_27] = 3
linear[Hexaedron_27] = False



def CheckIntegrity():
    return "ok"
