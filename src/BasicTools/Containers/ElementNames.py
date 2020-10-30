# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


class GeoSupport(object):
    def __init__(self,data):
        super(GeoSupport,self).__init__()
        self.name = data[0]
        self.dimensionality = data[1]
    def __rep__(self):
        res = "GeoSuport( " + self.name + ")"
        return res
    def __str__(self):
        return self.__rep__()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, GeoSupport):
            return self.name == other.name
        return False

    def __hash__(self):
        return id(self.name)

#GeoEntities = [
#        GeoSupport(("point",1)),   #0
#        GeoSupport(("bar"  ,1)),   #1
#        GeoSupport(("tri"  ,2)),   #2
#        GeoSupport(("quad" ,2)),   #3
#        GeoSupport(("tet"  ,3)),   #4
#        GeoSupport(("pyr"  ,3)),   #5
#        GeoSupport(("wed"  ,3)),   #6
#        GeoSupport(("hex"  ,3 ))]  #7


GeoPoint = GeoSupport(("point",0))   #0
GeoBar   = GeoSupport(("bar"  ,1))   #1
GeoTri   = GeoSupport(("tri"  ,2))   #2
GeoQuad  = GeoSupport(("quad" ,2))   #3
GeoTet   = GeoSupport(("tet"  ,3))   #4
GeoPyr   = GeoSupport(("pyr"  ,3))   #5
GeoWed   = GeoSupport(("wed"  ,3))   #6
GeoHex   = GeoSupport(("hex"  ,3 ))  #

numberOfNodes = {}

# permutation of index to make a valid element again after a mirror operation
mirrorPermutation = {}
dimension = {}
linear = {}
faces = {}
faces2 = {}
geoSupport = {}

#0d
Point_1  = 'point1'
geoSupport[Point_1] = GeoPoint
numberOfNodes[Point_1] = 1
mirrorPermutation[Point_1] = [0]
dimension[Point_1] = 0
linear[Point_1] = True
#1d
#linear
Bar_2 = 'bar2'
geoSupport[Bar_2] = GeoBar
numberOfNodes[Bar_2] = 2
mirrorPermutation[Bar_2] = [1,0]
dimension[Bar_2] = 1
linear[Bar_2] = True
faces[Bar_2] = [(Point_1,[0]), (Point_1,[1])]

#quadratic
Bar_3 = 'bar3'
geoSupport[Bar_3] = GeoBar
numberOfNodes[Bar_3] = 3
mirrorPermutation[Bar_3] = [1,0,2]
dimension[Bar_3] = 1
linear[Bar_3] = False
faces[Bar_3] = [(Point_1,[0]), (Point_1,[1])]

#2d
#linear
Triangle_3 = 'tri3'
geoSupport[Triangle_3] = GeoTri
numberOfNodes[Triangle_3] = 3
mirrorPermutation[Triangle_3] = [0,2,1]
dimension[Triangle_3] = 2
linear[Triangle_3] = True
faces[Triangle_3] = [(Bar_2,[0, 1]),
                     (Bar_2,[1, 2]),
                     (Bar_2,[2, 0])]

Quadrangle_4  = 'quad4'
geoSupport[Quadrangle_4] = GeoQuad
numberOfNodes[Quadrangle_4] = 4
mirrorPermutation[Quadrangle_4] = [1,0,3,2]
dimension[Quadrangle_4] = 2
linear[Quadrangle_4] = True
faces[Quadrangle_4] = [(Bar_2,[0, 1]),
                       (Bar_2,[1, 2]),
                       (Bar_2,[2, 3]),
                       (Bar_2,[3, 0])]

#quadratic
Triangle_6 = 'tri6'
geoSupport[Triangle_6] = GeoTri
numberOfNodes[Triangle_6] = 6
mirrorPermutation[Triangle_6] = [0,2,1,5,4,3]
dimension[Triangle_6] = 2
linear[Triangle_6] = False
faces[Triangle_6] = [(Bar_3,[0, 1,3]),
                     (Bar_3,[1, 2,4]),
                     (Bar_3,[2, 0,5])]

Quadrangle_8  = 'quad8'
geoSupport[Quadrangle_8] = GeoQuad
numberOfNodes[Quadrangle_8] = 8
mirrorPermutation[Quadrangle_8] = [0,2,1,5,4,3]
dimension[Quadrangle_8] = 2
linear[Quadrangle_8] = False
faces[Quadrangle_8] = [(Bar_3,[0, 1,4]),
                       (Bar_3,[1, 2,5]),
                       (Bar_3,[2, 3,6]),
                       (Bar_3,[3, 0,7])]

Quadrangle_9  = 'quad9'
geoSupport[Quadrangle_9] = GeoQuad
numberOfNodes[Quadrangle_9] = 9
dimension[Quadrangle_9] = 2
linear[Quadrangle_9] = False
faces[Quadrangle_9] = [(Bar_3,[0, 1,4]),
                       (Bar_3,[1, 2,5]),
                       (Bar_3,[2, 3,6]),
                       (Bar_3,[3, 0,7])]
#3d
#linear
Tetrahedron_4 = 'tet4'
geoSupport[Tetrahedron_4] = GeoTet
numberOfNodes[Tetrahedron_4] = 4
mirrorPermutation[Tetrahedron_4] = [0,2,1,3]
dimension[Tetrahedron_4] = 3
linear[Tetrahedron_4] = True
faces[Tetrahedron_4] = [(Triangle_3,[0, 2, 1]),
                        (Triangle_3,[0, 1, 3]),
                        (Triangle_3,[1, 2, 3]),
                        (Triangle_3,[2, 0, 3]),
     ]
faces2[Tetrahedron_4] = [(Bar_2,[0, 1]),
                         (Bar_2,[1, 2]),
                         (Bar_2,[2, 0]),
                         (Bar_2,[0, 3]),
                         (Bar_2,[1, 3]),
                         (Bar_2,[2, 3])
     ]

Pyramid_5  = 'pyr5'
geoSupport[Pyramid_5] = GeoPyr
numberOfNodes[Pyramid_5] = 5
dimension[Pyramid_5] = 3
linear[Pyramid_5] = False
faces[Pyramid_5] = [(Quadrangle_4,[0, 1, 2,3]),
                        (Triangle_3,[0, 1, 4]),
                        (Triangle_3,[1, 2, 4]),
                        (Triangle_3,[2, 3, 4]),
                        (Triangle_3,[3, 0, 4]),
     ]
Wedge_6 = 'wed6'
geoSupport[Wedge_6] = GeoWed
numberOfNodes[Wedge_6] = 6
dimension[Wedge_6] = 3
linear[Wedge_6] = True

Hexaedron_8 = 'hex8'
geoSupport[Hexaedron_8] = GeoHex
numberOfNodes[Hexaedron_8] = 8
dimension[Hexaedron_8] = 3
linear[Hexaedron_8] = True
faces[Hexaedron_8] = [(Quadrangle_4,[3, 0, 4, 7]),
                      (Quadrangle_4,[1, 2, 6, 5]),
                      (Quadrangle_4,[0, 1, 5, 4]),
                      (Quadrangle_4,[2, 3, 7, 6]),
                      (Quadrangle_4,[0, 3, 2, 1]),
                      (Quadrangle_4,[4, 5, 6, 7]),
     ]
faces2[Hexaedron_8] = [(Bar_2,[0,1]),
                       (Bar_2,[1,2]),
                       (Bar_2,[2,3]),
                       (Bar_2,[3,0]),

                       (Bar_2,[4,5]),
                       (Bar_2,[5,6]),
                       (Bar_2,[6,7]),
                       (Bar_2,[7,4]),

                       (Bar_2,[0,4]),
                       (Bar_2,[1,5]),
                       (Bar_2,[2,6]),
                       (Bar_2,[3,7]),
                           ]

#quadratic
Tetrahedron_10 = 'tet10'
geoSupport[Tetrahedron_10] = GeoTet
numberOfNodes[Tetrahedron_10] = 10
dimension[Tetrahedron_10] = 3
linear[Tetrahedron_10] = False
faces[Tetrahedron_10] = [(Triangle_6,[0, 2, 1, 6, 5, 4]),
                         (Triangle_6,[0, 1, 3, 4, 8, 7]),
                         (Triangle_6,[1, 2, 3, 5, 9, 8]),
                         (Triangle_6,[2, 0, 3, 6, 7, 9])]

Pyramid_13  = 'pyr13'
geoSupport[Pyramid_13] = GeoPyr
numberOfNodes[Pyramid_13] = 13
dimension[Pyramid_13] = 3
linear[Pyramid_13] = False

Wedge_15 = 'wed15'
geoSupport[Wedge_15] = GeoWed
numberOfNodes[Wedge_15] = 15
dimension[Wedge_15] = 3
linear[Wedge_15] = False

Wedge_18 = 'wed18'
geoSupport[Wedge_18] = GeoWed
numberOfNodes[Wedge_18] = 18
dimension[Wedge_18] = 3
linear[Wedge_18] = False


Hexaedron_20 = 'hex20'
geoSupport[Hexaedron_20] = GeoHex
numberOfNodes[Hexaedron_20] = 20
dimension[Hexaedron_20] = 3
linear[Hexaedron_20] = False
faces[Hexaedron_20] = [(Quadrangle_8,[3, 0, 4, 7,11,16,15,19]),
                       (Quadrangle_8,[1, 2, 6, 5, 9,18,13,17]),
                       (Quadrangle_8,[0, 1, 5, 4, 8,17,12,16]),
                       (Quadrangle_8,[2, 3, 7, 6,10,19,14,18]),
                       (Quadrangle_8,[0, 3, 2, 1,11,10, 9, 8]),
                       (Quadrangle_8,[4, 5, 6, 7,12,13,14,15])]


Hexaedron_27 = 'hex27'
geoSupport[Hexaedron_27] = GeoHex
numberOfNodes[Hexaedron_27] = 27
dimension[Hexaedron_27] = 3
linear[Hexaedron_27] = False
faces[Hexaedron_27] = [(Quadrangle_9,[3, 0, 4, 7,11,16,15,19,20]),
                       (Quadrangle_9,[1, 2, 6, 5, 9,18,13,17,21]),
                       (Quadrangle_9,[0, 1, 5, 4, 8,17,12,16,22]),
                       (Quadrangle_9,[2, 3, 7, 6,10,19,14,18,23]),
                       (Quadrangle_9,[0, 3, 2, 1,11,10, 9, 8,24]),
                       (Quadrangle_9,[4, 5, 6, 7,12,13,14,15,25])]


def CheckIntegrity(GUI=False):
    print(GeoPoint)
    print(GeoPoint==GeoBar)
    print(GeoPoint==1)
    print(GeoPoint!=GeoBar)
    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=True))
