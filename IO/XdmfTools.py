# -*- coding: utf-8 -*-
import BasicTools.Containers.ElementNames as EN

__author__ = "Felipe Bordeu"


XdmfName = {}
XdmfName[EN.Point_1] = 'Polyvertex'
XdmfName[EN.Bar_2] = 'Polyline'
XdmfName[EN.Triangle_3] = 'Triangle'
XdmfName[EN.Quadrangle_4] = 'Quadrilateral'
XdmfName[EN.Tetrahedron_4] ="Tetrahedron"
XdmfName[EN.Pyramid_5] = 'Pyramid'
XdmfName[EN.Wedge_6] = 'Wedge'
XdmfName[EN.Hexaedron_8] = 'Hexahedron'


XdmfName[EN.Bar_3] = "Edge_3"
XdmfName[EN.Triangle_6] = 'Triangle_6'
XdmfName[EN.Quadrangle_9] = 'Quadrilateral_9'
XdmfName[EN.Quadrangle_8] = 'Quadrilateral_8'
XdmfName[EN.Tetrahedron_10] = 'Tetrahedron_10'
XdmfName[EN.Pyramid_13] = 'Pyramid_13'
XdmfName[EN.Wedge_15] = 'Wedge_13'
XdmfName[EN.Wedge_18] = 'Wedge_18'
XdmfName[EN.Hexaedron_20] = 'Hexahedron_20'

XdmfNumber = {}
XdmfNumber[EN.Point_1] = 0x1
XdmfNumber[EN.Bar_2] = 0x2
XdmfNumber[EN.Triangle_3] = 0x4
XdmfNumber[EN.Quadrangle_4] = 0x5
XdmfNumber[EN.Tetrahedron_4] = 0x6
XdmfNumber[EN.Pyramid_5] = 0x7
XdmfNumber[EN.Wedge_6] = 0x8
XdmfNumber[EN.Hexaedron_8] = 0x9

XdmfNumber[EN.Bar_3] = 0x22
XdmfNumber[EN.Triangle_6] = 0x24
XdmfNumber[EN.Quadrangle_9] = 0x23
XdmfNumber[EN.Quadrangle_8] = 0x25
XdmfNumber[EN.Tetrahedron_10] = 0x26
XdmfNumber[EN.Pyramid_13] = 0x27
XdmfNumber[EN.Wedge_15] = 0x28
XdmfNumber[EN.Wedge_18] = 0x29
XdmfNumber[EN.Hexaedron_20] = 0x30

XdmfNumberToEN = {v:k for k,v in XdmfNumber.items()}
XdmfNameToEN = {v:k for k,v in XdmfName.items()}

class FieldNotFound(ValueError):
     """Exception to treat Field Not found """

     def __init__(self, value):
         self.value = 'Field "' + value + '" not found, Sorry!!'

     def __str__(self):
         return repr(self.value) # pragma: no cover


def CheckIntegrity():

    FieldNotFound('toto');
    return 'OK'

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
