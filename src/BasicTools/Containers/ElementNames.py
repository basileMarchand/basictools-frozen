# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

# If this file is modified a compilation must executed to export this data
# to the cpp portion of BasicTools
import numpy as np
from typing import Dict, Tuple, List

class GeoSupport():
    """Class to store basic information about the geometrical support.
    This class in not intender for end user.
    """
    def __init__(self,data: Tuple[str,int])-> None:
        super().__init__()
        self.name = data[0]
        self.dimensionality = data[1]

    def __rep__(self) -> str:
        res = "GeoSuport( " + self.name + ")"
        return res

    def __str__(self) -> str:
        return self.__rep__()

    def __eq__(self, other:object)-> bool:
        """Overrides the default implementation"""
        if isinstance(other, GeoSupport):
            return self.name == other.name
        return False

    def __hash__(self) -> int:
        return id(self.name)

GeoPoint = GeoSupport(("point",0))   #0
"""Point geometrical support : name point, dimensionality = 0"""
GeoBar   = GeoSupport(("bar"  ,1))   #1
"""Bar geometrical support : name bar, dimensionality = 1"""
GeoTri   = GeoSupport(("tri"  ,2))   #2
"""Triangle geometrical support : name tri, dimensionality = 2"""
GeoQuad  = GeoSupport(("quad" ,2))   #3
"""Quadrangle geometrical support : name quad, dimensionality = 2"""
GeoTet   = GeoSupport(("tet"  ,3))   #4
"""Tetrahedral geometrical support : name tet, dimensionality = 3"""
GeoPyr   = GeoSupport(("pyr"  ,3))   #5
"""Pyramidal (square base) geometrical support : name pyr, dimensionality = 3"""
GeoWed   = GeoSupport(("wed"  ,3))   #6
"""Wedge geometrical support : name wed, dimensionality = 3"""
GeoHex   = GeoSupport(("hex"  ,3 ))  #7
"""Hexahedral geometrical support : name hex, dimensionality = 3"""

class ElementInformation():
    """Class to store information about the different type of elements.
    """
    def __init__(self, name, geoSupport):
        self.name = name
        """Element Name """

        self.geoSupport = geoSupport
        """ the geometrical support of this element """
        self.dimension = geoSupport.dimensionality
        """Dimensionality 0,1,2,3 """

        self.numberOfNodes = 0
        """number of nodes of this elements"""
        self.linear = True
        """ true if the the jacobian matrix is constant in the element"""
        self.degree = 1
        """ the higher degree of the polynomial """
        self.faces:List[Tuple[str, np.ndarray]]  = []
        """First level faces. (triangles for a tetra for example)
        A list containing a tuples of the elements name and the local connectivity """
        self.faces2 = []
        """Second level faces. (bars for a tetra for example)
        A list containing a tuples of the elements name and the local connectivity """
        self.faces3 = []
        """third level faces. (points for a tetra for example)
        A list containing a tuples of the elements name and the local connectivity """
        self.mirrorPermutation = []
        """Permutation of index to make a valid element again after a mirror operation"""


numberOfNodes = {}     # type: Dict[str, int]
mirrorPermutation = {} # type: Dict[str, List[int]]
dimension = {}         # type: Dict[str, int]
linear = {}            # type: Dict[str, bool]
degree = {}            # type: Dict[str, int]
faces  = {}            # type: Dict[str, List[Tuple[str, List[int] ]] ]
faces2 = {}            # type: Dict[str, List[Tuple[str, List[int] ]] ]
faces3 = {}            # type: Dict[str, List[Tuple[str, List[int] ]] ]
geoSupport = {}        # type: Dict[str, GeoSupport]

#0D
Point_1  = 'point1'
geoSupport[Point_1] = GeoPoint
numberOfNodes[Point_1] = 1
mirrorPermutation[Point_1] = [0]
dimension[Point_1] = geoSupport[Point_1].dimensionality
linear[Point_1] = True
degree[Point_1] = 0
faces[Point_1] = []
faces2[Point_1] = []
faces3[Point_1] = []

#1D
#linear
Bar_2 = 'bar2'
geoSupport[Bar_2] = GeoBar
numberOfNodes[Bar_2] = 2
mirrorPermutation[Bar_2] = [1,0]
dimension[Bar_2] = 1
linear[Bar_2] = True
degree[Bar_2] = 1
faces[Bar_2] = [(Point_1,[0]),
                (Point_1,[1])]
faces2[Bar_2] = []
faces3[Bar_2] = []

#quadratic
Bar_3 = 'bar3'
geoSupport[Bar_3] = GeoBar
numberOfNodes[Bar_3] = 3
mirrorPermutation[Bar_3] = [1,0,2]
dimension[Bar_3] = 1
linear[Bar_3] = False
degree[Bar_3] = 2
faces[Bar_3] = [(Point_1,[0]),
                (Point_1,[1])]
faces2[Bar_3] = []
faces3[Bar_3] = []

#2D
#linear
Triangle_3 = 'tri3'
geoSupport[Triangle_3] = GeoTri
numberOfNodes[Triangle_3] = 3
mirrorPermutation[Triangle_3] = [0,2,1]
dimension[Triangle_3] = 2
linear[Triangle_3] = True
degree[Triangle_3] = 1
faces[Triangle_3] = [(Bar_2,[0, 1]),
                     (Bar_2,[1, 2]),
                     (Bar_2,[2, 0])]
faces2[Triangle_3] = [(Point_1,[0]),
                      (Point_1,[1]),
                     (Point_1,[2])]
faces3[Triangle_3] = []

#non linear
Quadrangle_4  = 'quad4'
geoSupport[Quadrangle_4] = GeoQuad
numberOfNodes[Quadrangle_4] = 4
mirrorPermutation[Quadrangle_4] = [1,0,3,2]
dimension[Quadrangle_4] = 2
linear[Quadrangle_4] = False
degree[Quadrangle_4] = 1
faces[Quadrangle_4] = [(Bar_2,[0, 1]),
                       (Bar_2,[1, 2]),
                       (Bar_2,[2, 3]),
                       (Bar_2,[3, 0])]
faces2[Quadrangle_4] = [(Point_1,[0]),
                       (Point_1,[1]),
                       (Point_1,[2]),
                       (Point_1,[3])]
faces3[Quadrangle_4] = []

#quadratic
Triangle_6 = 'tri6'
geoSupport[Triangle_6] = GeoTri
numberOfNodes[Triangle_6] = 6
mirrorPermutation[Triangle_6] = [0,2,1,5,4,3]
dimension[Triangle_6] = 2
linear[Triangle_6] = False
degree[Triangle_6] = 2
faces[Triangle_6] = [(Bar_3,[0, 1,3]),
                     (Bar_3,[1, 2,4]),
                     (Bar_3,[2, 0,5])]
faces2[Triangle_6] = [(Point_1,[0]),
                      (Point_1,[1]),
                      (Point_1,[2])]
faces3[Triangle_6] = []

Quadrangle_8  = 'quad8'
geoSupport[Quadrangle_8] = GeoQuad
numberOfNodes[Quadrangle_8] = 8
mirrorPermutation[Quadrangle_8] = [0,3,2,1,7,6,5,4]
dimension[Quadrangle_8] = 2
linear[Quadrangle_8] = False
degree[Quadrangle_8] = 2
faces[Quadrangle_8] = [(Bar_3,[0, 1,4]),
                       (Bar_3,[1, 2,5]),
                       (Bar_3,[2, 3,6]),
                       (Bar_3,[3, 0,7])]
faces2[Quadrangle_8] = [(Point_1,[0]),
                        (Point_1,[1]),
                        (Point_1,[2]),
                        (Point_1,[3])]
faces3[Quadrangle_8] = []

Quadrangle_9  = 'quad9'
geoSupport[Quadrangle_9] = GeoQuad
numberOfNodes[Quadrangle_9] = 9
mirrorPermutation[Quadrangle_9] = [0,3,2,1,7,6,5,4,8]
dimension[Quadrangle_9] = 2
linear[Quadrangle_9] = False
degree[Quadrangle_9] = 2
faces[Quadrangle_9] = [(Bar_3,[0, 1,4]),
                       (Bar_3,[1, 2,5]),
                       (Bar_3,[2, 3,6]),
                       (Bar_3,[3, 0,7])]
faces2[Quadrangle_9] = [(Point_1,[0]),
                       (Point_1,[1]),
                       (Point_1,[2]),
                       (Point_1,[3])]
faces3[Quadrangle_9] = []

#3D
#linear
Tetrahedron_4 = 'tet4'
geoSupport[Tetrahedron_4] = GeoTet
numberOfNodes[Tetrahedron_4] = 4
mirrorPermutation[Tetrahedron_4] = [0,2,1,3]
dimension[Tetrahedron_4] = 3
linear[Tetrahedron_4] = True
degree[Tetrahedron_4] = 1
faces[Tetrahedron_4] = [(Triangle_3,[0, 2, 1]),
                        (Triangle_3,[0, 1, 3]),
                        (Triangle_3,[1, 2, 3]),
                        (Triangle_3,[2, 0, 3])]
faces2[Tetrahedron_4] = [(Bar_2,[0, 1]),
                         (Bar_2,[1, 2]),
                         (Bar_2,[2, 0]),
                         (Bar_2,[0, 3]),
                         (Bar_2,[1, 3]),
                         (Bar_2,[2, 3])]
faces3[Tetrahedron_4] = [(Point_1,[0]),
                         (Point_1,[1]),
                         (Point_1,[2]),
                         (Point_1,[3])]

#non linear
Pyramid_5  = 'pyr5'
geoSupport[Pyramid_5] = GeoPyr
numberOfNodes[Pyramid_5] = 5
mirrorPermutation[Pyramid_5] = [0,3,2,1,4]
dimension[Pyramid_5] = 3
linear[Pyramid_5] = False
degree[Pyramid_5] = 1
faces[Pyramid_5] = [(Quadrangle_4,[0, 3, 2,1]),
                        (Triangle_3,[0, 1, 4]),
                        (Triangle_3,[1, 2, 4]),
                        (Triangle_3,[2, 3, 4]),
                        (Triangle_3,[3, 0, 4])]
faces2[Pyramid_5] = [(Bar_2, [0, 1]),
                     (Bar_2, [1, 2]),
                     (Bar_2, [2, 3]),
                     (Bar_2, [3, 0]),
                     (Bar_2, [0, 4]),
                     (Bar_2, [1, 4]),
                     (Bar_2, [2, 4]),
                     (Bar_2, [3, 4])]
faces3[Pyramid_5] = [(Point_1,[0]),
                     (Point_1,[1]),
                     (Point_1,[2]),
                     (Point_1,[3]),
                     (Point_1,[4])]

Wedge_6 = 'wed6'
geoSupport[Wedge_6] = GeoWed
numberOfNodes[Wedge_6] = 6
mirrorPermutation[Wedge_6] = [0,2,1,3,5,4]
dimension[Wedge_6] = 3
linear[Wedge_6] = False
degree[Wedge_6] = 1
faces[Wedge_6] = [(Triangle_3,[0, 2, 1]),
                  (Triangle_3,[3, 4, 5]),
                  (Quadrangle_4,[0, 1, 4,3]),
                  (Quadrangle_4,[1, 2, 5,4]),
                  (Quadrangle_4,[0, 3, 5,2])]
faces2[Wedge_6] = [(Bar_2,[0, 1]),
                   (Bar_2,[1, 2]),
                   (Bar_2,[2, 0]),
                   (Bar_2,[0, 3]),
                   (Bar_2,[1, 4]),
                   (Bar_2,[2, 5]),
                   (Bar_2,[3, 4]),
                   (Bar_2,[4, 5]),
                   (Bar_2,[5, 3])]
faces3[Wedge_6] = [(Point_1,[0]),
                   (Point_1,[1]),
                   (Point_1,[2]),
                   (Point_1,[3]),
                   (Point_1,[4]),
                   (Point_1,[5])]

Hexaedron_8 = 'hex8'
geoSupport[Hexaedron_8] = GeoHex
numberOfNodes[Hexaedron_8] = 8
mirrorPermutation[Hexaedron_8] = [0,3,2,1,4,7,6,5]
dimension[Hexaedron_8] = 3
linear[Hexaedron_8] = False
degree[Hexaedron_8] = 1
faces[Hexaedron_8] = [(Quadrangle_4,[3, 0, 4, 7]),
                      (Quadrangle_4,[1, 2, 6, 5]),
                      (Quadrangle_4,[0, 1, 5, 4]),
                      (Quadrangle_4,[2, 3, 7, 6]),
                      (Quadrangle_4,[0, 3, 2, 1]),
                      (Quadrangle_4,[4, 5, 6, 7])]
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
                       (Bar_2,[3,7])]

faces3[Hexaedron_8] = [(Point_1,[0]),
                       (Point_1,[1]),
                       (Point_1,[2]),
                       (Point_1,[3]),
                       (Point_1,[4]),
                       (Point_1,[5]),
                       (Point_1,[6]),
                       (Point_1,[7])]
#quadratic
Tetrahedron_10 = 'tet10'
geoSupport[Tetrahedron_10] = GeoTet
numberOfNodes[Tetrahedron_10] = 10
mirrorPermutation[Tetrahedron_10] = [0,2,1,3,6,5,4,7,9,8]
dimension[Tetrahedron_10] = 3
linear[Tetrahedron_10] = False
degree[Tetrahedron_10] = 2
faces[Tetrahedron_10] = [(Triangle_6,[0, 2, 1, 6, 5, 4]),
                         (Triangle_6,[0, 1, 3, 4, 8, 7]),
                         (Triangle_6,[1, 2, 3, 5, 9, 8]),
                         (Triangle_6,[2, 0, 3, 6, 7, 9])]

faces2[Tetrahedron_10] = [(Bar_3,[0, 1, 4]),
                          (Bar_3,[1, 2, 5]),
                          (Bar_3,[2, 0, 6]),
                          (Bar_3,[0, 3, 7]),
                          (Bar_3,[1, 3, 8]),
                          (Bar_3,[2, 3, 9])]
faces3[Tetrahedron_10] = [(Point_1,[0]),
                          (Point_1,[1]),
                          (Point_1,[2]),
                          (Point_1,[3]),
                          (Point_1,[4])]

Pyramid_13  = 'pyr13'
geoSupport[Pyramid_13] = GeoPyr
numberOfNodes[Pyramid_13] = 13
mirrorPermutation[Pyramid_13] = [0,3,2,1,4,8,7,6,5,9,12,11,10]
dimension[Pyramid_13] = 3
linear[Pyramid_13] = False
degree[Pyramid_13] = 2
faces[Pyramid_13] = [(Quadrangle_8,[0, 3, 2, 1, 8, 7, 6, 5]),
                     (Triangle_6,  [0, 1, 4, 5, 10, 9]),
                     (Triangle_6,  [1, 2, 4, 6, 11, 10]),
                     (Triangle_6,  [2, 3, 4, 7, 12, 11]),
                     (Triangle_6,  [3, 0, 4, 8,  9, 12])]
faces2[Pyramid_13] = [(Bar_2,[0, 1, 5]),
                      (Bar_2,[1, 2, 6]),
                      (Bar_2,[2, 3, 7]),
                      (Bar_2,[3, 0, 8]),
                      (Bar_2,[0, 4, 9]),
                      (Bar_2,[1, 4, 10]),
                      (Bar_2,[2, 4, 11]),
                      (Bar_2,[3, 4, 12])]
faces3[Pyramid_13] = [(Point_1,[0]),
                      (Point_1,[1]),
                      (Point_1,[2]),
                      (Point_1,[3]),
                      (Point_1,[4])]

Wedge_15 = 'wed15'
geoSupport[Wedge_15] = GeoWed
numberOfNodes[Wedge_15] = 15
dimension[Wedge_15] = 3
linear[Wedge_15] = False
degree[Wedge_15] = 2
mirrorPermutation[Wedge_15] = [0,2,1,3,5,4,8,7,6,11,10,9,12,14,13]
faces[Wedge_15] =[(Triangle_6,[0, 2, 1, 8, 7, 6]),
                  (Triangle_6,[3, 4, 5, 9,10,11]),
                  (Quadrangle_8,[0,1,4,3,6,13,9,12]),
                  (Quadrangle_8,[0,3,5,2,12,11,14,8]),
                  (Quadrangle_8,[2,5,4,1,14,10,13,7])]
faces2[Wedge_15] =[(Bar_3,[0,1,6]),
                   (Bar_3,[1,2,7]),
                   (Bar_3,[2,0,8]),
                   (Bar_3,[3,4,9]),
                   (Bar_3,[4,5,10]),
                   (Bar_3,[5,3,11]),
                   (Bar_3,[0,3,12]),
                   (Bar_3,[1,4,13]),
                   (Bar_3,[2,5,14])]
faces3[Wedge_15] = [(Point_1,[0]),
                    (Point_1,[1]),
                    (Point_1,[2]),
                    (Point_1,[3]),
                    (Point_1,[4]),
                    (Point_1,[5])]

Wedge_18 = 'wed18'
geoSupport[Wedge_18] = GeoWed
numberOfNodes[Wedge_18] = 18
dimension[Wedge_18] = 3
linear[Wedge_18] = False
degree[Wedge_18] = 2
mirrorPermutation[Wedge_18] = [0,2,1,3,5,4,8,7,6,11,10,9,12,14,13,17,16,15]
faces[Wedge_18] =[(Triangle_6,[0, 2, 1, 8, 7, 6]),
                  (Triangle_6,[3, 4, 5, 9,10,11]),
                  (Quadrangle_9,[0,1,4,3,6,13,9,12,15]),
                  (Quadrangle_9,[0,3,5,2,12,11,14,8,17]),
                  (Quadrangle_9,[2,5,4,1,14,10,13,7,16])]
faces2[Wedge_18] =[(Bar_3,[0,1,6]),
                   (Bar_3,[1,2,7]),
                   (Bar_3,[2,0,8]),
                   (Bar_3,[3,4,9]),
                   (Bar_3,[4,5,10]),
                   (Bar_3,[5,3,11]),
                   (Bar_3,[0,3,12]),
                   (Bar_3,[1,4,13]),
                   (Bar_3,[2,5,14])]
faces3[Wedge_18] = [(Point_1,[0]),
                    (Point_1,[1]),
                    (Point_1,[2]),
                    (Point_1,[3]),
                    (Point_1,[4]),
                    (Point_1,[5])]

Hexaedron_20 = 'hex20'
geoSupport[Hexaedron_20] = GeoHex
numberOfNodes[Hexaedron_20] = 20
dimension[Hexaedron_20] = 3
linear[Hexaedron_20] = False
degree[Hexaedron_20] = 2
mirrorPermutation[Hexaedron_20] = [0,3,2,1,4,7,6,5,11,10,9,8,15,14,13,12,16,19,18,17]
faces[Hexaedron_20] = [(Quadrangle_8,[3, 0, 4, 7,11,16,15,19]),
                       (Quadrangle_8,[1, 2, 6, 5, 9,18,13,17]),
                       (Quadrangle_8,[0, 1, 5, 4, 8,17,12,16]),
                       (Quadrangle_8,[2, 3, 7, 6,10,19,14,18]),
                       (Quadrangle_8,[0, 3, 2, 1,11,10, 9, 8]),
                       (Quadrangle_8,[4, 5, 6, 7,12,13,14,15])]
faces2[Hexaedron_20] = [(Bar_3,[0,1,8]),
                        (Bar_3,[1,2,9]),
                        (Bar_3,[2,3,10]),
                        (Bar_3,[3,0,11]),
                        (Bar_3,[4,5,12]),
                        (Bar_3,[5,6,13]),
                        (Bar_3,[6,7,14]),
                        (Bar_3,[7,4,15]),
                        (Bar_3,[0,4,16]),
                        (Bar_3,[1,5,17]),
                        (Bar_3,[2,6,18]),
                        (Bar_3,[3,7,19])]
faces3[Hexaedron_20] = [(Point_1,[0]),
                        (Point_1,[1]),
                        (Point_1,[2]),
                        (Point_1,[3]),
                        (Point_1,[4]),
                        (Point_1,[5]),
                        (Point_1,[6]),
                        (Point_1,[7])]

Hexaedron_27 = 'hex27'
geoSupport[Hexaedron_27] = GeoHex
numberOfNodes[Hexaedron_27] = 27
dimension[Hexaedron_27] = 3
linear[Hexaedron_27] = False
degree[Hexaedron_27] = 2
mirrorPermutation[Hexaedron_27] = [0,3,2,1,4,7,6,5,11,10,9,8,15,14,13,12,16,19,18,17,20,23,21,22,24,25,26]
faces[Hexaedron_27] = [(Quadrangle_9,[3, 0, 4, 7,11,16,15,19,20]),
                       (Quadrangle_9,[1, 2, 6, 5, 9,18,13,17,21]),
                       (Quadrangle_9,[0, 1, 5, 4, 8,17,12,16,22]),
                       (Quadrangle_9,[2, 3, 7, 6,10,19,14,18,23]),
                       (Quadrangle_9,[0, 3, 2, 1,11,10, 9, 8,24]),
                       (Quadrangle_9,[4, 5, 6, 7,12,13,14,15,25])]
faces2[Hexaedron_27] = [(Bar_3,[0,1,8]),
                        (Bar_3,[1,2,9]),
                        (Bar_3,[2,3,10]),
                        (Bar_3,[3,0,11]),
                        (Bar_3,[4,5,12]),
                        (Bar_3,[5,6,13]),
                        (Bar_3,[6,7,14]),
                        (Bar_3,[7,4,15]),
                        (Bar_3,[0,4,16]),
                        (Bar_3,[1,5,17]),
                        (Bar_3,[2,6,18]),
                        (Bar_3,[3,7,19])]
faces3[Hexaedron_27] = [(Point_1,[0]),
                        (Point_1,[1]),
                        (Point_1,[2]),
                        (Point_1,[3]),
                        (Point_1,[4]),
                        (Point_1,[5]),
                        (Point_1,[6]),
                        (Point_1,[7])]

ElementsInfo:Dict[str, ElementInformation] = {}

"""Module variable to store ElementInformation for every type of element

The key is the element name,
The value is an instance of ElementInformation
"""

for name,geo in geoSupport.items():
    ei = ElementInformation(name,geo)
    ei.numberOfNodes = numberOfNodes[name]
    ei.linear = linear[name]
    ei.degree =  degree[name]
    ei.faces = faces[name]
    ei.faces2 = faces2[name]
    ei.faces3 = faces3[name]
    ei.mirrorPermutation = mirrorPermutation[name]
    ElementsInfo[name] = ei

def CheckIntegrity(GUI=False):
    """CheckIntegrity function. Tests

    Parameters
    ----------
    GUI : bool, optional
        if True, generate (in some case) an output on a new window, by default False

    Returns
    -------
    str
        ok if all ok
    """
    print(GeoPoint)
    print(GeoPoint==GeoBar)
    print(GeoPoint==1)
    print(GeoPoint!=GeoBar)
    print(hash(GeoPoint))
    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=True))
