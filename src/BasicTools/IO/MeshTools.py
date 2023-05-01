# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Mesh tools
"""

import BasicTools.Containers.ElementNames as EN

# ASCII elements names and tags names

ASCIIName = {}
#ASCIIName[EN.Point_1] = 'Vertices'
ASCIIName[EN.Bar_2] = 'Edges'
ASCIIName[EN.Triangle_3] = 'Triangles'
ASCIIName[EN.Quadrangle_4] = 'Quadrilaterals'
ASCIIName[EN.Tetrahedron_4] = 'Tetrahedra'
ASCIIName[EN.Hexaedron_8] = 'Hexahedra'
ASCIIName[EN.Pyramid_5] = "Pyramids"
ASCIIName[EN.Wedge_6] = "Prisms"
ASCIIName[EN.Hexaedron_20] = 'HexahedraP2'
ASCIIName[EN.Quadrangle_8] = 'QuadrilateralsP2'


ASCIITypes = {}
for types, name in ASCIIName.items():
    ASCIITypes[name] = types


Corners = "Corners"
Ridges = "Ridges"
RequiredVertices = "RequiredVertices"
RequiredEdges = "RequiredEdges"
RequiredTriangles = "RequiredTriangles"

# flag in the file -> (element tag and tagname)
ASCIITags = {}
ASCIITags["Ridges"] = (EN.Bar_2, Ridges)
ASCIITags['RequiredTriangles'] = (EN.Triangle_3, RequiredTriangles)
ASCIITags['RequiredEdges'] = (EN.Bar_2, RequiredEdges)


# binary elemnts number and tags numbers from file "libmesh5.h"
BinaryKeywords = {
    "GmfReserved1": 0,
    "GmfVersionFormatted": 1,
    "GmfReserved2": 2,
    "GmfDimension": 3,
    "GmfVertices": 4,
    "GmfEdges": 5,
    "GmfTriangles": 6,
    "GmfQuadrilaterals": 7,
    "GmfTetrahedra": 8,
    "GmfPentahedra": 9,
    "GmfHexahedra": 10,
    "GmfReserved3": 11,
    "GmfReserved4": 12,
    "GmfCorners": 13,
    "GmfRidges": 14,
    "GmfRequiredVertices": 15,
    "GmfRequiredEdges": 16,
    "GmfRequiredTriangles": 17,
    "GmfRequiredQuadrilaterals": 18,
    "GmfTangentAtEdgeVertices": 19,
    "GmfNormalAtVertices": 20,
    "GmfNormalAtTriangleVertices": 21,
    "GmfNormalAtQuadrilateralVertices": 22,
    "GmfAngleOfCornerBound": 23,
    "GmfTrianglesP2": 24,
    "GmfTrianglesP3": 25,
    "GmfTrianglesP4": 26,
    "GmfQuadrilateralsP2": 27,
    "GmfQuadrilateralsP3": 28,
    "GmfQuadrilateralsP4": 29,
    "GmfTetrahedraP2": 30,
    "GmfTetrahedraP3": 31,
    "GmfTetrahedraP4": 32,
    "GmfHexahedraP2": 33,
    "GmfHexahedraP3": 34,
    "GmfHexahedraP4": 35,
    "GmfReserved17": 36,
    "GmfReserved18": 37,
    "GmfReserved19": 38,
    "GmfReserved20": 39,
    "GmfReserved21": 40,
    "GmfReserved22": 41,
    "GmfReserved23": 42,
    "GmfReserved24": 43,
    "GmfReserved25": 44,
    "GmfReserved26": 45,
    "GmfReserved27": 46,
    "GmfReserved28": 47,
    "GmfReserved29": 48,
    "GmfReserved30": 49,
    "GmfBoundingBox": 50,
    "GmfReserved31": 51,
    "GmfReserved32": 52,
    "GmfReserved33": 53,
    "GmfEnd": 54,
    "GmfReserved34": 55,
    "GmfReserved35": 56,
    "GmfReserved36": 57,
    "GmfReserved37": 58,
    "GmfTangents": 59,
    "GmfNormals": 60,
    "GmfTangentAtVertices": 61,
    "GmfSolAtVertices": 62,
    "GmfSolAtEdges": 63,
    "GmfSolAtTriangles": 64,
    "GmfSolAtQuadrilaterals": 65,
    "GmfSolAtTetrahedra": 66,
    "GmfSolAtPentahedra": 67,
    "GmfSolAtHexahedra": 68,
    "GmfDSolAtVertices": 69,
    "GmfISolAtVertices": 70,
    "GmfISolAtEdges": 71,
    "GmfISolAtTriangles": 72,
    "GmfISolAtQuadrilaterals": 73,
    "GmfISolAtTetrahedra": 74,
    "GmfISolAtPentahedra": 75,
    "GmfISolAtHexahedra": 76,
    "GmfIterations": 77,
    "GmfTime": 78,
    "GmfReserved38": 79}

FieldTypes = {
    "GmfSca": 1,  # scalar size = 1
    "GmfVec": 2,  # vector size = 2
    "GmfSymMat": 3,  # symmetric mat = dim*(dim +1)/2
    "GmfMat": 4  # full matrix = dim*dim
}

BinaryNumber = {}
BinaryNumber[EN.Point_1] = BinaryKeywords["GmfVertices"]
BinaryNumber[EN.Bar_2] = BinaryKeywords["GmfEdges"]
BinaryNumber[EN.Triangle_3] = BinaryKeywords["GmfTriangles"]
BinaryNumber[EN.Quadrangle_4] = BinaryKeywords["GmfQuadrilaterals"]
BinaryNumber[EN.Tetrahedron_4] = BinaryKeywords["GmfTetrahedra"]
BinaryNumber[EN.Hexaedron_8] = BinaryKeywords["GmfHexahedra"]
BinaryNumber[EN.Pyramid_5] = BinaryKeywords["GmfReserved30"]
BinaryNumber[EN.Wedge_6] = BinaryKeywords["GmfPentahedra"]
BinaryNumber[EN.Hexaedron_20] = BinaryKeywords["GmfHexahedraP2"]
BinaryNumber[EN.Quadrangle_8] = BinaryKeywords["GmfQuadrilateralsP2"]

BinaryName = {}
BinaryName[EN.Point_1] = "GmfVertices"
BinaryName[EN.Bar_2] = "GmfEdges"
BinaryName[EN.Triangle_3] = "GmfTriangles"
BinaryName[EN.Quadrangle_4] = "GmfQuadrilaterals"
BinaryName[EN.Tetrahedron_4] = "GmfTetrahedra"
BinaryName[EN.Hexaedron_8] = "GmfHexahedra"
BinaryName[EN.Pyramid_5] = "GmfReserved30"
BinaryName[EN.Wedge_6] = "GmfPentahedra"
BinaryName[EN.Hexaedron_20] = "GmfHexahedraP2"
BinaryName[EN.Quadrangle_8] = "GmfQuadrilateralsP2"


BinaryTypes = {}
for types, number in BinaryNumber.items():
    BinaryTypes[number] = types

BinaryTags = {}
BinaryTags[BinaryKeywords["GmfRidges"]] = (EN.Bar_2, Ridges)
BinaryTags[BinaryKeywords["GmfRequiredEdges"]] = (EN.Bar_2, RequiredEdges)
BinaryTags[BinaryKeywords["GmfRequiredTriangles"]] = (EN.Triangle_3, RequiredTriangles)

BinaryFields = {}
BinaryFields[BinaryKeywords["GmfSolAtVertices"]] = ("SolAtVertices")

def GetTypesForVersion(version:int):
    """get data types from mesh format version

    Parameters
    ----------
    version : int
        format version

    Returns
    -------
    tuple, tuple, tuple
        size, format and type for pos, int and float
    """
    import numpy as np

    posSize = 4
    posFormat = 'i'
    posType = np.int32

    intSize = 4
    intFormat = 'i'
    intType = np.int32

    floatSize = 4
    floatFormat = 'f'
    floatType = np.float32


    if version > 1:
        floatSize = 8
        floatFormat = 'd'
        floatType = np.float64

    if version > 2:
        posSize = 8
        posFormat = 'q'
        posType = np.int64

    if version > 3:
        intSize = 8
        intFormat = 'q'
        intType = np.int64

    return (posSize,posFormat,posType), (intSize,intFormat,intType), (floatSize, floatFormat, floatType)

def CheckIntegrity():
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())  # pragma: no cover
