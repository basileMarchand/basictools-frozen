# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Ansys Workbench Mechanical input file reader
"""

import numpy as np

import BasicTools.Containers.ElementNames as EN
import BasicTools.Containers.UnstructuredMesh as UM
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Helpers.ParserHelper import LocalVariables
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType
from BasicTools.IO.AnsysTools import PermutationAnsysToBasicTools

def ReadAnsys(fileName=None, string=None, out=None, **kwargs):
    """Function API for reading an Ansys result file

    Parameters
    ----------
    fileName : str, optional
        name of the file to be read, by default None
    string : str, optional
        data to be read as a string instead of a file, by default None
    out : UnstructuredMesh, optional
        output unstructured mesh object containing reading result, by default None

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    reader = AnsysReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string, out=out, **kwargs)
    return reader.output


class AnsysReader(ReaderBase):
    """Ansys Reader class
    """

    def __init__(self):
        super(AnsysReader, self).__init__()
        self.commentChar = '!'
        self.readFormat = 'r'
        self.encoding = "latin-1"

    def Read(self, fileName=None, string=None, out=None):
        """Function that performs the reading of an Ansys result file

        Parameters
        ----------
        fileName : str, optional
            name of the file to be read, by default None
        string : str, optional
            data to be read as a string instead of a file, by default None
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
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
        self.CSDict = {}    #  to store co

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
                'cm': self.ParseCMDefinition,
                '*Get':self.ParseGet,
                'N': self.ParseN,
                'CS': self.ParseCS,
                'EMODIF': self.PaserEMODIF,

                }

        def Pass(args): pass

        for multiline in self.iterator:
            lines = multiline.split('$')
            for line in lines:
                line = line.strip()
                tokens = line.split(',')
                keyword = tokens[0]
                arguments = tokens[1:]
                commands.get(keyword, Pass)(arguments)


    def ParseAssignment(self, args):
        self.substitutions.SetVariable(args[0], args[1])

    def ParseGet(self, args):
        #https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_cmd/Hlp_C_GET.html#get.prep.node
        par = str(args[0])        # The name of the resulting parameter. See *SET for name restrictions.
        entity = str(args[1])     # Entity keyword. Valid keywords are NODE, ELEM, KP, LINE, AREA, VOLU, etc., as shown for Entity = in the tables below.
        ENTNUM = int(args[2])     # The number or label for the entity (as shown for ENTNUM = in the tables below). In some cases, a zero (or blank) ENTNUM represents all entities of the set.
        Item1 =  str(args[3])     # The name of a particular item for the given entity. Valid items are as shown in the Item1 columns of the tables below.
        IT1NUM =  str([4])     # The number (or label) for the specified Item1 (if any). Valid IT1NUM values are as shown in the IT1NUM columns of the tables below. Some Item1 labels do not require an IT1NUM value.
        if len(args) >5 :
            Item2 =  str(args[5])     # A second set of item labels and numbers to further qualify the item for which data are to be retrieved. Most items do not require this level of information.
            IT2NUM =  str([6])

        if entity == "NODE" and ENTNUM == 0 and Item1 == "NUM" and IT1NUM == "MAXD":
            #normaly I must use the max id node
            self.substitutions.SetVariable(par,self.result.nodes.shape[0])
        else:
            return

    def ParseN(self,args):
        NID= (self.substitutions.Apply(args[0]))
        NX = float(self.substitutions.Apply(args[1]))
        NY = float(self.substitutions.Apply(args[2]))
        NZ = float(self.substitutions.Apply(args[3]))
        self.CSDict["_Internal_N_"+NID] = [NX,NY,NZ]

    def ParseCS(self,args):
        KCN = int(args[0]) # Arbitrary reference number assigned to this coordinate system. Must be greater than 10. A coordinate system previously defined with this number will be redefined.
        KCS = args[1] # Coordinate system type: 0 or CART - Cartesian
                      #                         1 or CYLIN - Cylindrical (circular or elliptical)
                      #                         2 or SPHE   — Spherical (or spheroidal)
                      #                         3 or TORO   — Toroidal
        NORIG = self.CSDict["_Internal_N_"+self.substitutions.Apply(args[2])] # Node defining the origin of this coordinate system. If NORIG = P, graphical picking is enabled and all remaining command fields are ignored (valid only in the GUI).
        NXAX = self.CSDict["_Internal_N_"+self.substitutions.Apply(args[3])]  # Node defining the positive x-axis orientation of this coordinate system.
        NXYPL = self.CSDict["_Internal_N_"+self.substitutions.Apply(args[4])] # Node defining the x-y plane (with NORIG and NXAX) in the first or second quadrant of this coordinate system.
        if len(args) > 5:
            PAR1 = args[5]  # Used for elliptical, spheroidal, or toroidal systems. If KCS = 1 or 2, PAR1 is the ratio of the ellipse Y-axis radius to X-axis radius (defaults to 1.0 (circle)). If KCS = 3, PAR1 is the major radius of the torus.
            PAR2 = args[6]  # Used for spheroidal systems. If KCS = 2, PAR2 = ratio of ellipse Z-axis radius to X-axis radius (defaults to 1.0 (circle)).
        if KCS == "0" or KCS == "CART":
            origin =np.array(NORIG)
            V0 = np.array(NXAX)
            V0 /= np.linalg.norm(V0)
            V1 = np.array(NXYPL)
            V1 -= np.dot(V0,V1)*V0
            V1 /= np.linalg.norm(V1)
            self.CSDict[KCN] = [origin, V0,V1]
        else:
            raise

    def PaserEMODIF(self, args):
        #https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_cmd/Hlp_C_EMODIF.html
        IEL = args[0] #Modify nodes and/or attributes for element number IEL. If ALL, modify all selected elements [ESEL]. If IEL = P, graphical picking is enabled and all remaining command fields are ignored (valid only in the GUI). A component name may also be substituted for IEL.
        STLOC = args[1]
        if "V0" not in self.result.elemFields:
            self.result.elemFields["origin"] = np.zeros((self.result.GetNumberOfElements(),3), dtype=PBasicFloatType )
            self.result.elemFields["V0"] = np.zeros((self.result.GetNumberOfElements(),3), dtype=PBasicFloatType )
            self.result.elemFields["V0"][:,0:3] = [1,0,0]
            self.result.elemFields["V1"] = np.zeros((self.result.GetNumberOfElements(),3), dtype=PBasicFloatType )
            self.result.elemFields["V1"][:,0:3] = [0,1,0]

        if STLOC == "ESYS":
            I1 = int(args[2])
            element_type, rank = self.element_type_and_rank_from_id[int(IEL)]
            container = self.result.GetElementsOfType(element_type)
            self.result.ComputeGlobalOffset()

            origin, V0, V1 = self.CSDict[I1]
            self.result.elemFields["origin"][container.globaloffset + rank,:] =  origin
            self.result.elemFields["V0"][container.globaloffset + rank,:] =  V0
            self.result.elemFields["V1"][container.globaloffset + rank,:] =  V1

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
        self.element_type_and_rank_from_id[element_id] = (EN.Point_1, internal_rank)
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
    # May degenerate to EN.Wedge_15, EN.Pyramid_13 or EN.Tetrahedron_10
    # Node numbering: ijklmnopqrstuvwxyzab
    repeated_kl = nodes[2] == nodes[3]
    repeated_mn = nodes[4] == nodes[5]
    repeated_op = nodes[6] == nodes[7]
    from itertools import compress
    if repeated_op:
        if repeated_mn:
            if repeated_kl:
                assert(nodes[4] == nodes[6])
                assert(nodes[6] == nodes[12])
                assert(nodes[6] == nodes[13])
                assert(nodes[6] == nodes[14])
                assert(nodes[6] == nodes[15])
                assert(nodes[18] == nodes[19])
                assert(nodes[2] == nodes[10])
                internal_element_type = EN.Tetrahedron_10
                unique_nodes = list(compress(nodes, (1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0)))
                nodes = [ nodes[per] for per in PermutationAnsysToBasicTools[EN.Tetrahedron_10] ]
            else:
                assert(nodes[4] == nodes[6])
                assert(nodes[6] == nodes[12])
                assert(nodes[6] == nodes[13])
                assert(nodes[6] == nodes[14])
                assert(nodes[6] == nodes[15])
                internal_element_type = EN.Pyramid_13
                unique_nodes = list(compress(nodes, (1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1)))
        else:
            assert(nodes[6] == nodes[14])
            assert(nodes[18] == nodes[19])
            assert(nodes[2] == nodes[10])
            internal_element_type = EN.Wedge_15
            unique_nodes = list(compress(nodes, (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0)))
    else:
        internal_element_type = EN.Hexaedron_20
        unique_nodes = nodes
    return internal_element_type, unique_nodes

def discriminate_solid187(nodes):
    internal_element_type = EN.Tetrahedron_10
    nodes = [ nodes[per] for per in PermutationAnsysToBasicTools[EN.Tetrahedron_10] ]
    return internal_element_type, nodes

# FOLLW201 is a one-node 3d element used to apply nodal forces

internal_element_type_from_ansys = {
        '170': discriminate_tri_or_quad,
        '174': discriminate_tri_or_quad,
        '154': discriminate_tri_or_quad,
        '181': discriminate_shell181,
        '185': discriminate_solid185,
        '186': discriminate_solid186,
        '187': discriminate_solid187
        }


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".ansys", AnsysReader)


def CheckIntegrity():
    __teststring = u"""
nblock,3,,44
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
     6782     7.261344291E-04    -6.938080072E-04     7.640761699E-02
     6805     2.186338673E-04    -8.078558619E-04     7.646145635E-02
     7308     6.305392899E-04    -2.000856285E-04     7.594041518E-02
     7309     1.108636472E-04    -3.288991121E-04     7.599344534E-02
     7859     1.725123756E-04    -7.844719866E-04     7.599047369E-02
     7860     7.591255417E-04    -7.164507578E-04     7.593581544E-02
     8108     1.476547648E-04    -3.576023233E-04     7.647352256E-02
     8111     6.989548354E-04    -2.490909499E-04     7.643098938E-02
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
    45628     4.723841482E-04    -7.508319345E-04     7.643453667E-02
    45632     7.426299854E-04    -7.051293825E-04     7.617171621E-02
    45633     7.125446322E-04    -4.714494786E-04     7.641930319E-02
    45769     1.955731215E-04    -7.961639242E-04     7.622596502E-02
    45770     1.831443160E-04    -5.827290926E-04     7.646748946E-02
    47969     3.707014685E-04    -2.644923703E-04     7.596693026E-02
    47970     6.948324158E-04    -4.582681931E-04     7.593811531E-02
    47971     6.647470626E-04    -2.245882892E-04     7.618570228E-02
    47973     1.416880114E-04    -5.566855494E-04     7.599195951E-02
    47974     1.292592060E-04    -3.432507177E-04     7.623348395E-02
    50074     4.658189586E-04    -7.504613722E-04     7.596314456E-02
    50908     4.233048001E-04    -3.033466366E-04     7.645225597E-02
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
et,4,186
eblock,19,solid,,4
(19i9)
        3        4        1        3        0        0        0        0       20        0        1     6782     7860     7859     6805     8111     7308     7309     8108
    45632    50074    45769    45628    47971    47969    47974    50908    45633    47970    47973    45770
        3        4        1        3        0        0        0        0       20        0        2     6782     7860     7859     7859     8111     8111     8111     8111
    45632    50074     7859    45628     8111     8111     8111     8111    45633    47970    47973    47973
        3        4        1        3        0        0        0        0       20        0        3     6782     7860     7859     6805     8111     8111     8111     8111
    45632    50074    45769    45628     8111     8111     8111     8111    45633    47970    47973    45770
        3        4        1        3        0        0        0        0       20        0        4     6782     7860     7859     7859     8111     7308     7309     7309
    45632    50074     7859    45628    47971    47969     7309    50908    45633    47970    47973    47973
-1
"""

    res = ReadAnsys(string=__teststring)

    print("----")
    print('coords: {}'.format(res.nodes))
    print('node ids: {}'.format(res.originalIDNodes))
    print('hex20: {}'.format((res.GetElementsOfType('hex20').connectivity)))
    print('pyr13: {}'.format((res.GetElementsOfType('pyr13').connectivity)))
    print('tet4: {}'.format((res.GetElementsOfType('tet4').connectivity)))
    print('tet10: {}'.format((res.GetElementsOfType('tet10').connectivity)))
    print('tri6: {}'.format((res.GetElementsOfType('tri6').connectivity)))
    print('wed15: {}'.format((res.GetElementsOfType('wed15').connectivity)))
    node_tag = res.GetNodalTag('FewNodes')
    print('node set {}: {}'.format(node_tag, node_tag.GetIds()))
    for t in ('et_1', 'et_2', 'et_3', 'et_4', 'rc_4', 'rc_5', 'EB_0', 'EB_1', 'EB_2', 'EB_3'):
        print('element set {}: {}'.format(t, res.GetElementsInTag(t)))

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
