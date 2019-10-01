# -*- coding: utf-8 -*-
""" Ansys Workbench Mechanical input file reader

"""
import numpy as np

import BasicTools.Containers.ElementNames as EN
import BasicTools.Containers.UnstructuredMesh as UM
from BasicTools.IO.ReaderBase import ReaderBase


def ReadAnsys(fileName=None, string=None, out=None, **kwargs):
    reader = AnsysReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string, out=out, **kwargs)
    return reader.output

class AnsysReader(ReaderBase):
    def __init__(self):
        super(AnsysReader, self).__init__()
        self.commentChar = '!'
        self.readFormat = 'r'

    def Read(self, fileName=None, string=None, out=None):
        if fileName is not None:
            self.SetFileName(fileName)
        if string is not None:
            self.SetStringToRead(string)
        self.StartReading()
        res = UM.UnstructuredMesh() if out is None else out

        node_rank_from_id = {}
        element_type_and_rank_from_id = {}
        tagsNames = []
        element_type_ids = dict() # Local to global element type numbering

        while(True):
            line = self.ReadCleanLine()
            if not line:
                break

            if line.startswith('nblock'):
                # nblock, NUMFIELD, Solkey, NDMAX[, NDSEL]
                tokens = line.split(',')
                field_count = int(tokens[1])
                solid_key = tokens[2]
                max_node_count = int(tokens[3])

                res.nodes = np.empty((max_node_count, field_count))
                res.originalIDNodes = np.empty((max_node_count,), dtype=np.int)

                assert(solid_key == '')
                # Skip format line
                line = self.ReadCleanLine()

                node_rank = 0
                while True:
                    line = self.ReadCleanLine()
                    if line.startswith('-1'):
                        break
                    assert(node_rank < max_node_count)
                    tokens = line.split()
                    node_id = int(tokens[0])
                    node_rank_from_id[node_id] = node_rank
                    res.originalIDNodes[node_rank] = node_id
                    res.nodes[node_rank, :] = [float(t) for t in tokens[1:]]
                    node_rank += 1
                continue

            if line.startswith('et,'):
                tokens = line.split(',')
                element_type_ids[int(tokens[1])] = tokens[2]
                continue

            if line.startswith('eblock'):
                # eblock, NUM_NODES, Solkey[,,count]
                tokens = line.split(',')
                # int(tokens[1]) -> num_nodes: Cannot be trusted
                solid_key = tokens[2]
                max_element_count = int(tokens[4])

                # Skip format line
                line = self.ReadCleanLine()

                if solid_key == 'solid':
                    self.ReadSolidEblock(res, element_type_ids, node_rank_from_id, element_type_and_rank_from_id, max_element_count)
                else:
                    self.ReadNonSolidEblock(res, element_type_ids, node_rank_from_id, element_type_and_rank_from_id, max_element_count)
                continue

            if line.startswith('CMBLOCK'):
                tokens = line.split(',')
                tag_name = tokens[1]
                kind = tokens[2]
                item_count = int(tokens[3])

                items_per_line = 8
                full_line_count = item_count // items_per_line
                remainder = item_count % items_per_line
                line_count = full_line_count if remainder == 0 \
                        else full_line_count + 1

                # Skip format line
                line = self.ReadCleanLine()

                items = list()
                for i in range(line_count):
                    line = self.ReadCleanLine()
                    items.extend((int(t) for t in line.split()))

                if kind == 'NODE':
                    tag = res.nodesTags.CreateTag(tag_name)
                    tag.SetIds([node_rank_from_id[n] for n in items])
                else:
                    assert(kind == 'ELEMENT')
                    for e in items:
                        t, r = element_type_and_rank_from_id[e]
                        res.AddElementToTag(r, tag_name)
                continue

            # Ugly hack to handle nodal forces
            if line.startswith('type,'):
                tokens = line.split(',')
                element_type_id = element_type_ids[int(tokens[1])]
                if element_type_id == '201':
                    # FOLLW201 is a one-node 3d element used to apply nodal forces
                    line = self.ReadCleanLine()
                    tokens = line.split(',')
                    assert(len(tokens) == 3 and tokens[0] == 'en')
                    node_id = int(tokens[2])
                    # Figure out a nodal tag name from the element id
                    tag_name = 'AtElem_' + tokens[1]
                    tag = res.nodesTags.CreateTag(tag_name)
                    tag.SetIds([node_rank_from_id[node_id]])
                continue

        self.EndReading()
        res.PrepareForOutput()
        self.output = res
        return res

    def ReadSolidEblock(self, res, element_type_ids, node_rank_from_id, element_type_and_rank_from_id, max_element_count):
        element_rank = 0
        while True:
            line = self.ReadCleanLine()
            if line.startswith('-1'):
                break
            assert(element_rank < max_element_count)
            tokens = line.split()
            values = [int(t) for t in tokens]
            material_id = values[0] # unused
            element_id = values[10]
            element_node_count = values[8]
            if element_node_count > 8:
                overflow = self.ReadCleanLine()
                values.extend((int (t) for t in overflow.split()))
            nodes = values[11:11+element_node_count]
            element_type_id = element_type_ids[values[1]]
            internal_element_type, unique_nodes = \
                    internal_element_type_from_ansys[element_type_id](nodes)
            connectivity = [node_rank_from_id[n] for n in unique_nodes]
            elements = res.GetElementsOfType(internal_element_type)
            internal_count = elements.AddNewElement(connectivity, element_id)
            internal_rank = internal_count - 1
            element_type_and_rank_from_id[element_id] = \
                    (internal_element_type, internal_rank)
            auto_tag = 'eblock_{}'.format(values[1])
            res.AddElementToTagUsingOriginalId(element_id, auto_tag)
            element_rank += 1

    def ReadNonSolidEblock(self, res, element_type_ids, node_rank_from_id, element_type_and_rank_from_id, max_element_count):
        element_rank = 0
        while True:
            line = self.ReadCleanLine()
            if line.startswith('-1'):
                break
            assert(element_rank < max_element_count)
            tokens = line.split()
            values = [int(t) for t in tokens]
            element_id = values[0]
            element_properties = values[1:4]
            # Ensure that all 3 properties are identical
            assert(element_properties[:-1] == element_properties[1:])
            element_type_id = element_type_ids[element_properties[0]]
            nodes = values[5:]
            internal_element_type, unique_nodes = \
                    internal_element_type_from_ansys[element_type_id](nodes)
            connectivity = [node_rank_from_id[n] for n in unique_nodes]
            elements = res.GetElementsOfType(internal_element_type)
            internal_count = elements.AddNewElement(connectivity, element_id)
            internal_rank = internal_count - 1
            element_type_and_rank_from_id[element_id] = \
                    (internal_element_type, internal_rank)
            auto_tag = 'eblock_{}'.format(element_properties[0])
            res.AddElementToTagUsingOriginalId(element_id, auto_tag)
            element_rank += 1


def discriminate_surf154(nodes):
    # SURF154: EN.Quadrangle_4 or EN.Quadrangle_9
    # May degenerate to EN.Triangle_3 or EN.Triangle_6
    # Node numbering: ijklmnop
    repeated_kl = nodes[2] == nodes[3]
    from itertools import compress
    if len(nodes) == 4:
        if repeated_kl:
            internal_element_type = EN.Triangle_3
            unique_nodes = nodes[:-1]
        else:
            internal_element_type = EN.Quadrangle_4
            unique_nodes = nodes
    else:
        assert(len(nodes) == 8)
        if repeated_kl:
            assert(nodes[2] == nodes[6])
            internal_element_type = EN.Triangle_6
            unique_nodes = list(compress(nodes, (1, 1, 1, 0, 1, 1, 0, 1)))
        else:
            internal_element_type = EN.Quadrangle_8
            unique_nodes = nodes
    return internal_element_type, unique_nodes

def discriminate_solid185(nodes):
    # SOLID185: EN.Hexaedron_8
    # May degenerate to EN.Wedge_6, EN.Pyramid_5 or EN_Tetrahedron_4
    # Node numbering: ijklmnop
    repeated_kl = nodes[2] == nodes[3]
    repeated_mn = nodes[4] == nodes[5]
    repeated_op = nodes[6] == nodes[7]
    from itertools import compress
    if repeated_op:
        if repeated_mn:
            if repeated_kl:
                internal_element_type = EN.Tetrahedron_4
                unique_nodes = compress(nodes, (1, 1, 1, 0, 1))
            else:
                internal_element_type = EN.Pyramid_5
                unique_nodes = compress(nodes, (1, 1, 1, 1, 1))
        else:
            internal_element_type = EN.Wedge_6
            unique_nodes = compress(nodes, (1, 1, 1, 0, 1, 1, 1))
    else:
        internal_element_type = EN.Hexaedron_8
        unique_nodes = nodes
    return internal_element_type, list(unique_nodes)

def discriminate_solid186(nodes):
    # SOLID186: EN.Hexaedron_20
    # May degenerate to wedge, pyramid or tetrahedron
    implemented = False
    assert(implemented)

def discriminate_solid187(nodes):
    # SOLID187: EN.Tetrahedron_10
    return EN.Tetrahedron_10, nodes


internal_element_type_from_ansys = {
        '154': discriminate_surf154,
        '185': discriminate_solid185,
        '187': discriminate_solid187
        }


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".ansys", AnsysReader)


def CheckIntegrity():
    __teststring = u"""
nblock,3,,24
(1i9,3e20.9e3)
      551     7.784032421E-02     6.661491953E-02     2.000000000E-01
     1691     8.991484887E-02     6.820190681E-02     1.903279204E-01
     1944     7.455965603E-02     6.107661302E-02     1.906569403E-01
     2111     9.395317480E-02     5.598013490E-02     1.886333684E-01
     2218     8.604121057E-02     6.526570146E-02     1.798559428E-01
     2233     8.587154519E-02     5.451195787E-02     1.957137802E-01
     4975     8.387758654E-02     6.740841317E-02     1.951639602E-01
     4976     7.619999012E-02     6.384576627E-02     1.953284701E-01
     4978     8.185593470E-02     6.056343870E-02     1.978568901E-01
    11353     8.223725245E-02     6.463925992E-02     1.904924303E-01
    11355     9.193401184E-02     6.209102085E-02     1.894806444E-01
    11357     8.797802972E-02     6.673380414E-02     1.850919316E-01
    11358     8.789319703E-02     6.135693234E-02     1.930208503E-01
    12938     8.030043330E-02     6.317115724E-02     1.852564415E-01
    12939     8.021560061E-02     5.779428544E-02     1.931853602E-01
    13639     8.999719269E-02     6.062291818E-02     1.842446556E-01
    13640     8.991236000E-02     5.524604638E-02     1.921735743E-01
    13703     8.595637788E-02     5.988882967E-02     1.877848615E-01
    32654    -1.274217310E+01     3.840702614E+01    -1.612452772E+01
    37816    -1.334739993E+01     3.905933630E+01    -1.603912279E+01
    37856    -1.364793556E+01     3.797073871E+01    -1.607483826E+01
    39901    -1.425619495E+01     3.874218666E+01    -1.561731625E+01
    42378    -1.371046977E+01     3.892752029E+01    -1.493156312E+01
    62174    -1.331961824E+01     3.837984154E+01    -1.550789885E+01
-1
et,1,185
eblock,19,solid,,340744
(19i9)
        1        1        1        1        0        0        0        0        8        0        1    37816    39901    62174    62174    42378    42378    42378    42378
        1        1        1        1        0        0        0        0        8        0        2    37816    37856    62174    62174    39901    39901    39901    39901
        1        1        1        1        0        0        0        0        8        0        3    32654    62174    37816    37816    37856    37856    37856    37856
-1
CMBLOCK,FewNodes,NODE,      2
(8i10)
     37856     37816
et,2,187
eblock,19,solid,,3
(19i9)
        1        2        1        1        0        0        0        0       10        0        1     1691     2233     1944     2218    11358    12939    11353    11357
    13703    12938
        1        2        1        1        0        0        0        0       10        0        2     1691     2111     2233     2218    11355    13640    11358    11357
    13639    13703
        1        2        1        1        0        0        0        0       10        0        3      551     1691     2233     1944     4975    11358     4978     4976
    11353    12939
-1
et,3,154
eblock,10,,,2
(15i9)
    10915        3        3        3       12      551     1691     2233     2233     4975    11358     2233     4976
    10916        3        3        3       12     1691     2233     1944     1944    11358    12939     1944    11357
-1
"""

    res = ReadAnsys(string=__teststring)

    print("----")
    print('coords: {}'.format(res.nodes))
    print('node ids: {}'.format(res.originalIDNodes))
    print('tet4: {}'.format((res.GetElementsOfType('tet4').connectivity)))
    print('tet10: {}'.format((res.GetElementsOfType('tet10').connectivity)))
    print('tri6: {}'.format((res.GetElementsOfType('tri6').connectivity)))
    print('node tags: {}'.format(res.GetNodalTag('FewNodes')))
    print('node sets: {}'.format(res.GetNodalTag('FewNodes').GetIds()))

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
