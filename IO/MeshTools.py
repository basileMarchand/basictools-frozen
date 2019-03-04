# -*- coding: utf-8 -*-



import BasicTools.Containers.ElementNames as EN

## ASCII elements names and tags names

ASCIIName = {}
#ASCIIName[EN.Point_1] = 'Vertices'
ASCIIName[EN.Bar_2]         = 'Edges'
ASCIIName[EN.Triangle_3]    = 'Triangles'
ASCIIName[EN.Quadrangle_4]  = 'Quadrilaterals'
ASCIIName[EN.Tetrahedron_4] = 'Tetrahedra'
ASCIIName[EN.Hexaedron_8]   = 'Hexahedra'
ASCIIName[EN.Hexaedron_20]   = 'HexahedraP2'
ASCIIName[EN.Quadrangle_8]  = 'QuadrilateralsP2'

ASCIITypes = {}
for types,name in ASCIIName.items():
    ASCIITypes[name] = types


Corners = "Corners"
Ridges = "Ridges"
RequiredVertices = "RequiredVertices"
RequiredEdges ="RequiredEdges"
RequiredTriangles ="RequiredTriangles"

# flag in the file -> (element tag and tagname)
ASCIITags = {}
ASCIITags["Ridges"] =  (EN.Bar_2,Ridges)
ASCIITags['RequiredTriangles'] = (EN.Triangle_3,RequiredTriangles)
ASCIITags['RequiredEdges'] = (EN.Bar_2,RequiredEdges)


## binary elemnts number and tags numbers


BinaryKeywords ={
"GmfReserved1": 0 ,
"GmfVersionFormatted": 1 ,
"GmfReserved2": 2 ,
"GmfDimension": 3 ,
"GmfVertices": 4 ,
"GmfEdges": 5 ,
"GmfTriangles": 6 ,
"GmfQuadrilaterals": 7 ,
"GmfTetrahedra": 8 ,
"GmfPentahedra": 9 ,
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
"GmfQuadrilateralsP3": 27,
"GmfQuadrilateralsP4": 29,
"GmfTetrahedraP2": 30,
"GmfTetrahedraP3": 31,
"GmfTetrahedraP4": 32,
"GmfHexahedraP2": 33,
"GmfHexahedraP3": 34,
"GmfHexahedraP4": 35,
"GmfReserved17": 36,
"GmfReserved18": 37,
"GmfReserved19": 37,
"GmfReserved20": 39,
"GmfReserved21": 40,
"GmfReserved22": 41,
"GmfReserved23": 42,
"GmfReserved24": 43,
"GmfReserved25": 44,
"GmfReserved26": 45,
"GmfReserved27": 46,
"GmfReserved28": 47,
"GmfReserved29": 47,
"GmfReserved30": 49,
"GmfBoundingBox": 50,
"GmfReserved31": 51,
"GmfReserved32": 52,
"GmfReserved33": 53,
"GmfEnd": 54,
"GmfReserved34": 55,
"GmfReserved35": 56,
"GmfReserved36": 57,
"GmfReserved37": 57,
"GmfTangents": 59,
"GmfNormals": 60,
"GmfTangentAtVertices": 61,
"GmfSolAtVertices": 62,
"GmfSolAtEdges": 63,
"GmfSolAtTriangles": 64,
"GmfSolAtQuadrilaterals": 65,
"GmfSolAtTetrahedra": 66,
"GmfSolAtPentahedra": 67,
"GmfSolAtHexahedra": 67,
"GmfDSolAtVertices": 69,
"GmfISolAtVertices": 70,
"GmfISolAtEdges": 71,
"GmfISolAtTriangles": 72,
"GmfISolAtQuadrilaterals": 73,
"GmfISolAtTetrahedra": 74,
"GmfISolAtPentahedra": 75,
"GmfISolAtHexahedra": 76,
"GmfIterations": 77,
"GmfTime": 77,
"GmfReserved38": 79}



BinaryNumber = {}
BinaryNumber[EN.Point_1] = BinaryKeywords["GmfVertices"]
BinaryNumber[EN.Bar_2] = BinaryKeywords["GmfEdges"]

BinaryNumber[EN.Triangle_3] = BinaryKeywords["GmfTriangles"]
BinaryNumber[EN.Quadrangle_4] = BinaryKeywords["GmfQuadrilaterals"]
BinaryNumber[EN.Tetrahedron_4] = BinaryKeywords["GmfTetrahedra"]
BinaryNumber[EN.Hexaedron_8] = BinaryKeywords["GmfHexahedra"]


BinaryTypes = {}
for types,number in BinaryNumber.items():
    BinaryTypes[number] = types

#"GmfRequiredVertices": 15,
#"GmfRequiredEdges": 16,
#"GmfRequiredTriangles": 17,

BinaryTags = {}
BinaryTags[BinaryKeywords["GmfRidges"]] =(EN.Bar_2,Ridges)
BinaryTags[BinaryKeywords["GmfRequiredEdges"]] =(EN.Bar_2,RequiredEdges)
BinaryTags[BinaryKeywords["GmfRequiredTriangles"]] =(EN.Triangle_3,RequiredTriangles)

BinaryFields ={}
BinaryFields[BinaryKeywords["GmfSolAtVertices"]] =("SolAtVertices")
#BinaryTags[BinaryKeywords["GmfRequiredEdges"]] =(EN.Bar_2,RequiredEdges)
#BinaryTags[BinaryKeywords["GmfRequiredTriangles"]] =(EN.Bar_2,RequiredTriangles)

def CheckIntegrity():
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover