# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" Ansys Workbench Mechanical input file reader

"""
import numpy as np

import BasicTools.Containers.ElementNames as EN
import BasicTools.Containers.UnstructuredMesh as UM
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Helpers.ParserHelper import LocalVariables
from BasicTools.NumpyDefs import PBasicIndexType


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
        self.encoding = "latin-1"

    def Read(self, fileName=None, string=None, out=None):
        if fileName is not None:
            self.SetFileName(fileName)
        if string is not None:
            self.SetStringToRead(string)

        with self.GetIterator() as iterator:
            session = Session(iterator, out)

        result = session.GetMesh()
        result.PrepareForOutput()
        self.output = result
        return result

    def GetIterator(self):
        return InputContextManager(self)


class InputContextManager:
    def __init__(self, reader):
        self.reader = reader

    def __enter__(self):
        self.reader.StartReading()
        return InputIterator(self.reader)

    def __exit__(self, exc_type, exc_value, traceback):
        self.reader.EndReading()
        return False


class InputIterator:
    def __init__(self, reader):
        self.reader = reader

    def __iter__(self):
        return self

    def __next__(self):
        line = self.reader.ReadCleanLine()
        if not line:
            raise StopIteration()
        return line

class Session:
    def __init__(self, iterator, out=None):
        self.iterator = iterator
        self.result = UM.UnstructuredMesh() if out is None else out
        self.block_count = 0

        # Nodes
        self.result.nodes = np.empty((0, 3), dtype=np.double)
        self.result.originalIDNodes = np.empty((0,), dtype=PBasicIndexType)
        self.node_count = 0
        self.node_rank_from_id = dict()

        # Elements
        self.element_type_and_rank_from_id = dict()
        self.element_type_ids = dict() # Local to global element type numbering
        self.current_element_type = None
        self.elementSelectoin = None

        # Variables
        self.substitutions = LocalVariables(prePostChars=('',''))

        commands = {
                '*set': self.ParseAssignment,
                'nblock': self.ParseNodeBlock,
                'eblock': self.ParseElementBlock,
                'en': self.ParseUnblockedElement,
                'et': self.ParseElementTypeDefinition,
                'type': self.ParseElementTypeSelection,
                'CMBLOCK': self.ParseTagDefinition,
                'esel': self.ParseEselDefinition,
                'cm': self.ParseCMDefinition
                }

        def Pass(args): pass

        for line in self.iterator:
            tokens = line.split(',')
            keyword = tokens[0]
            arguments = tokens[1:]
            commands.get(keyword, Pass)(arguments)

    def ParseAssignment(self, args):
        self.substitutions.SetVariable(args[0], args[1])

    def ParseNodeBlock(self, args):
        # Arguments: NUMFIELD,[Solkey,NDMAX[,NDSEL]]
        field_count = int(args[0])
        solid_key = args[1] if len(args) > 1 else ''
        max_new_node_count = int(args[2]) if len(args) > 2 else 1

        assert(field_count == 3)
        assert(solid_key == '')

        expected_node_count = self.node_count + max_new_node_count
        oldNodes = self.result.nodes
        self.result.nodes = np.zeros((expected_node_count,3))
        self.result.nodes[0:oldNodes.shape[0],:] = oldNodes

        oldOriginalIDNodes = self.result.originalIDNodes
        self.result.originalIDNodes = np.zeros(expected_node_count,dtype=int)
        self.result.originalIDNodes[0:oldOriginalIDNodes.shape[0]] = oldOriginalIDNodes

        # Skip format line
        line = next(self.iterator)

        while True:
            line = next(self.iterator)
            if line.startswith('-1'):
                break
            tokens = line.split()
            node_id = int(tokens[0])
            self.node_rank_from_id[node_id] = self.node_count
            self.result.originalIDNodes[self.node_count] = node_id
            self.result.nodes[self.node_count, :] = [float(t) for t in tokens[1:]]
            self.node_count += 1

    def ParseElementBlock(self, args):
        # Arguments: NUM_NODES,Solkey[,,count]
        # int(args[0]) -> num_nodes: Cannot be trusted
        solid_key = args[1]
        max_element_count = int(args[3])

        # Skip format line
        line = next(self.iterator)

        if solid_key == 'solid':
            self.ReadSolidEblock(max_element_count)
        else:
            self.ReadNonSolidEblock(max_element_count)

        self.block_count += 1

    def ParseUnblockedElement(self, args):
        et = self.current_element_type
        assert(self.element_type_ids[et] in ('170', '201'))
        element_id = int(args[0])
        node_id = int(self.substitutions.Apply(args[1]))
        node_rank = self.node_rank_from_id[node_id]
        elements = self.result.GetElementsOfType(EN.Point_1)
        internal_count = elements.AddNewElement([node_rank], element_id)
        internal_rank = internal_count - 1
        auto_etag = 'et_{}'.format(et)
        elements.AddElementToTag(internal_rank, auto_etag)
        # Figure out a nodal tag name from the element id
        auto_ntag = 'elem_{}'.format(element_id)
        tag = self.result.nodesTags.CreateTag(auto_ntag)
        tag.SetIds([self.node_rank_from_id[node_id]])

    def ParseElementTypeDefinition(self, args):
        et = int(self.substitutions.Apply(args[0]))
        self.element_type_ids[et] = args[1]

    def ParseElementTypeSelection(self, args):
        et = int(self.substitutions.Apply(args[0]))
        self.current_element_type = et

    def ParseCMDefinition(self,args):
        #cm,PEAU_MAT_COMPOSITE_SKIN_SPAR_B_3,elem
        if args[1].lower() != "elem":
            return
        tagname = args[0]
        for selection in self.elementSelection:
            for names, data in self.result.elements.items():
                if selection in data.tags:
                    data.tags.CreateTag(tagname,False).AddToTag(data.tags[selection].GetIds())

    def ParseEselDefinition(self,args):
        #esel,s,type,,2
        #esel,a,type,,3,18
        #esel,a,type,,22,29
        if args[0].lower() == "none":
            self.elementSelection = []
            return

        if args[0].lower() == "all":
            self.elementSelection = None
            return

        if args[0].lower() == "s":
            self.elementSelection = []
        elif args[0].lower() == "a":
            pass

        if args[1].lower() =="type":
            VMIN = int(args[3])
            if len(args) > 4:
                VMAX = int(args[4])
            else:
                VMAX = VMIN

            if len(args) > 5:
                VINC = int(args[5])
            else:
                VINC = 1

            for et in range(VMIN,VMAX+1,VINC):
                self.elementSelection.append(f"et_{et}")


    def ParseTagDefinition(self, args):
        tag_name = args[0].strip()
        kind = args[1]
        item_count = int(args[2])

        items_per_line = 8
        full_line_count = item_count // items_per_line
        remainder = item_count % items_per_line
        line_count = full_line_count if remainder == 0 \
                else full_line_count + 1

        # Skip format line
        line = next(self.iterator)

        items = list()
        for i in range(line_count):
            line = next(self.iterator)
            items.extend((int(t) for t in line.split()))

        if kind == 'NODE':
            tag = self.result.nodesTags.CreateTag(tag_name)
            tag.SetIds([self.node_rank_from_id[n] for n in items])
        else:
            assert(kind in ('ELEMENT','ELEM') )
            for element_id in items:
                element_type, rank = self.element_type_and_rank_from_id[element_id]
                container = self.result.GetElementsOfType(element_type)
                container.AddElementToTag(rank, tag_name)

    def GetMesh(self):
        return self.result

    def ReadEblock(self, data_parser, max_element_count):
        element_rank = 0
        while True:
            line = next(self.iterator)
            if line.startswith('-1'):
                break
            assert(element_rank < max_element_count)
            tokens = line.split()
            values = [int(t) for t in tokens]

            element_id, et, real_constant, _, nodes = \
                    data_parser(values, self.iterator)

            element_type_id = self.element_type_ids[et]
            internal_element_type, unique_nodes = \
                    internal_element_type_from_ansys[element_type_id](nodes)
            connectivity = [self.node_rank_from_id[n] for n in unique_nodes]
            elements = self.result.GetElementsOfType(internal_element_type)
            internal_count = elements.AddNewElement(connectivity, element_id)
            internal_rank = internal_count - 1
            self.element_type_and_rank_from_id[element_id] = \
                    (internal_element_type, internal_rank)
            auto_etags = (
                    'et_{}'.format(et),
                    'rc_{}'.format(real_constant),
                    'EB_{}'.format(self.block_count))
            for t in auto_etags:
                elements.AddElementToTag(internal_rank, t)
            element_rank += 1

    def ReadSolidEblock(self, max_element_count):
        def SolidDataParser(values, it):
            material_id = values[0]
            et = values[1]
            real_constant = values[2]
            element_id = values[10]
            element_node_count = values[8]
            if element_node_count > 8:
                overflow = next(it)
                values.extend((int (t) for t in overflow.split()))
            nodes = values[11:11+element_node_count]
            return (element_id, et, real_constant, material_id, nodes)

        self.ReadEblock(SolidDataParser, max_element_count)

    def ReadNonSolidEblock(self, max_element_count):
        def NonSolidDataParser(values, _):
            element_id = values[0]
            et = values[1]
            real_constant = values[2]
            material_id = values[3]
            nodes = values[5:]
            return (element_id, et, real_constant, material_id, nodes)

        self.ReadEblock(NonSolidDataParser, max_element_count)

def discriminate_tri_or_quad(nodes):
    # SURF154/TARGE170/CONTA174: EN.Quadrangle_4 or EN.Quadrangle_9
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

def discriminate_shell181(nodes):
    return discriminate_tri_or_quad(nodes)


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

# FOLLW201 is a one-node 3d element used to apply nodal forces

internal_element_type_from_ansys = {
        '170': discriminate_tri_or_quad,
        '174': discriminate_tri_or_quad,
        '154': discriminate_tri_or_quad,
        '181': discriminate_shell181,
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
*set,tid,3
et,tid,154
eblock,10,,,2
(15i9)
    10915        3        4        3       12      551     1691     2233     2233     4975    11358     2233     4976
    10916        3        5        3       12     1691     2233     1944     1944    11358    12939     1944    11357
-1
"""

    res = ReadAnsys(string=__teststring)

    print("----")
    print('coords: {}'.format(res.nodes))
    print('node ids: {}'.format(res.originalIDNodes))
    print('tet4: {}'.format((res.GetElementsOfType('tet4').connectivity)))
    print('tet10: {}'.format((res.GetElementsOfType('tet10').connectivity)))
    print('tri6: {}'.format((res.GetElementsOfType('tri6').connectivity)))
    node_tag = res.GetNodalTag('FewNodes')
    print('node set {}: {}'.format(node_tag, node_tag.GetIds()))
    for t in ('et_1', 'et_2', 'et_3', 'rc_4', 'rc_5', 'EB_0', 'EB_1', 'EB_2'):
        print('element set {}: {}'.format(t, res.GetElementsInTag(t)))

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
