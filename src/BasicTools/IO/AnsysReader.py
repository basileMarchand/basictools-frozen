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
                # Skip Format line
                line = self.ReadCleanLine()

                node_rank = 0
                while True:
                    line = self.ReadCleanLine()
                    if line.startswith('-1'):
                        break
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
                # entry_count = int(tokens[1])
                solid_key = tokens[2]
                max_element_count = int(tokens[4])

                assert(solid_key == 'solid')
                # Skip Format line
                line = self.ReadCleanLine()

                element_rank = 0
                while True:
                    line = self.ReadCleanLine()
                    if line.startswith('-1'):
                        break
                    assert(element_rank < max_element_count)
                    tokens = line.split()
                    values = [int(t) for t in tokens]
                    material_id = values[0]
                    element_id = values[10]
                    element_node_count = values[8]
                    nodes = values[11:11+element_node_count]
                    element_type_id = element_type_ids[values[1]]
                    internal_element_type, unique_nodes = \
                            ansys_element_types[element_type_id](nodes)
                    connectivity = [node_rank_from_id[n] for n in unique_nodes]
                    elements = res.GetElementsOfType(internal_element_type)
                    internal_count = elements.AddNewElement(connectivity, element_id)
                    internal_rank = internal_count - 1
                    element_type_and_rank_from_id[element_id] = \
                            (internal_element_type, internal_rank)
                    element_rank += 1
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

                # Skip Format line
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
                        t.AddElementToTag(r, tag_name)
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

def discriminate_solid185(nodes):
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

#Ansys Element types:
# SOLID185: EN.Hexaedron_8, may degenerate to EN.Wedge_6, EN.Pyramid_5 or EN_Tetrahedron_4
# SOLID186: EN.Hexaedron_20, may degenerate to wedge, pyramid or tetrahedron
# SOLID187: EN.Tetrahedron_10
ansys_element_types = {'185': discriminate_solid185}

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".ansys", AnsysReader)

def CheckIntegrity():

    __teststring = u"""
nblock,3,,6
(1i9,3e20.9e3)
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
"""
    res = ReadAnsys(string=__teststring)

    print("----")
    print(res.nodes)
    print(res.originalIDNodes)
    print(res.GetElementsOfType('tet4').connectivity)
    print(res.GetNodalTag('FewNodes'))
    print(res.GetNodalTag('FewNodes').GetIds())

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
