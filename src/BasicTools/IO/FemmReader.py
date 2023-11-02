# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Fem file reader
"""

import numpy as np
import re


from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType

import BasicTools.Containers.UnstructuredMesh as UM
import BasicTools.Containers.ElementNames as ElementNames

from BasicTools.IO.ReaderBase import ReaderBase


def ReadFemm(fileName=None,string=None):
    """Function for reading a Femm (ans) result file

    Parameters
    ----------
    fileName : str, optional
        name of the file to be read, by default None
    string : str, optional
        data to be read as a string instead of a file, by default None

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    obj = FemmReader()
    if fileName is not None:
        obj.SetFileName(fileName)

    if string is not None:
        obj.SetStringToRead(string)

    return obj.Read()

keysToIgnore = [ "[Format]",
                 "[Comment]",
                 "[Frequency]",
                 "[Precision]",
                 "[MinAngle]",
                 "[DoSmartMesh]",
                 "[Depth]",
                 "[ProblemType]",
                 "[Coordinates]",
                 "[ACSolver]",
                 "[PrevType]",
                 "[PrevSoln]",
                 "[Comment]",
                 "[PointProps]",
                 "[NumPoints]",
                 "[NumSegments]",
                 "[NumArcSegments]",
                 "[NumHoles]",
                 "[NumBlockLabels",    ]

class FemmReader(ReaderBase):
    """Geof Reader class
    """
    def __init__(self):
        super().__init__()
        self.readFormat = 'r'
        self.rescaleToMeters = True

    def Read(self, fileName=None,string=None, blockAsFields = False):

        import BasicTools.Containers.UnstructuredMesh as UM
        if fileName is not None:
            self.SetFileName(fileName)

        if string is not None:
            self.SetStringToRead(string)

        self.StartReading()

        nbIP  = 0

        res = UM.UnstructuredMesh()

        def ReadBlock():
            data = {}
            l = self.ReadCleanLine().strip()
            if not l.startswith("<Begin"):
                raise Exception(f"Line ({l}) is not the start of a block")
            while(True):
                l = self.ReadCleanLine().strip()
                if l.startswith("<End"):
                    break
                key, value = re.findall(r'(?:[^\s="]|"(?:\\.|[^"])*")+', l)
                #key, value = l.split("=")
                key = key.strip(" <>")
                value = value.strip(' "')
                data[key] = value
            return data

        factorsFromUnit = {"millimeters":0.001, "centimeters":0.01, "meters":1, "inches": 2.54*0.01}

        factor = 1.
        while(True):
            l = self.ReadCleanLine()
            if l == None :
                break
            if l[0] != "[":
                continue
            if any([ l.startswith(x) for x in  keysToIgnore ]):
                continue

            if l.startswith("[LengthUnits]"):
                factor = factorsFromUnit[(l.split("=")[1].strip())]
                continue

            if l.startswith("[BlockProps]"):
                self.BlockProps = []
                nbBlocs = int(l.split("=")[1].strip())
                for bn in range(nbBlocs):
                    self.BlockProps.append( ReadBlock() )
                continue

            if l.startswith("[BdryProps]"):
                self.BdryProps = []
                nbBlocs = int(l.split("=")[1].strip())
                for bn in range(nbBlocs):
                    self.BdryProps.append( ReadBlock() )
                continue

            if l.startswith("[CircuitProps]"):
                self.CircuitProps = []
                nbBlocs = int(l.split("=")[1].strip())
                for bn in range(nbBlocs):
                    self.CircuitProps.append( ReadBlock() )
                continue

            if l.startswith("[Solution]"):
                s = self.ReadCleanLine().split()
                nbNodes = int(s[0])
                self.PrintDebug("Reading "+str(nbNodes)+ " Nodes ")

                res.nodes = np.empty((nbNodes,3), dtype=PBasicFloatType)
                res.originalIDNodes= np.empty((nbNodes,),dtype=PBasicIndexType)
                for x in range(nbNodes):
                    s  = self.ReadCleanLine().split()
                    res.originalIDNodes[x] = x
                    res.nodes[x,:] = list(map(float,s[0:3]))

                res.nodeFields["A"] = np.copy(res.nodes[:,2])
                res.nodes[:,2] = 0

                s  = self.ReadCleanLine().split()
                nbElements = int(s[0])
                tris = res.GetElementsOfType(ElementNames.Triangle_3)
                bars = res.GetElementsOfType(ElementNames.Bar_2)
                tris.Allocate(nbElements)
                conn = tris.connectivity
                tris.originalIds = np.arange(nbElements)
                tri_faces = ElementNames.faces[ElementNames.Triangle_3]

                tags = np.empty(nbElements,dtype=PBasicIndexType)


                for x in range(nbElements):
                    l  = self.ReadCleanLine()
                    s = l.split()
                    conn[x,:]  = list(map(int,s[0:3]))
                    tags[x] = int(s[3])
                    for b, val in enumerate(map(int,s[4:7])):
                        if val == -1:
                            continue
                        n = bars.AddNewElement(conn[x,:][tri_faces[b][1]],0) -1
                        bars.tags.CreateTag("Bdry"+str(val),False).AddToTag(n)



                uniqueTags = np.unique(tags)
                if blockAsFields :
                    block = np.zeros(tris.cpt+bars.cpt, dtype=PBasicIndexType)
                    triOffset = res.ComputeGlobalOffset()[ElementNames.Triangle_3]
                    block[triOffset:triOffset+tris.cpt] = tags
                    res.elemFields["Blocks"] = block
                else:
                    for tn in uniqueTags:
                        tris.tags.CreateTag("Block"+str(tn)).SetIds(np.where(tags==tn)[0])

                continue

            raise Exception(f" don't know how to treat this line: {l}" )

        self.EndReading()
        if factor != 1. and self.rescaleToMeters:
            print("Changing units of the file to meters")
            res.nodes *= factor
        res.PrepareForOutput()
        return res


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".ans",FemmReader,)


def CheckIntegrity():

    data = """[Format]      =  4.0
[Frequency]   =  0
[Precision]   =  1e-008
[MinAngle]    =  20
[DoSmartMesh] =  0
[Depth]       =  100
[LengthUnits] =  millimeters
[ProblemType] =  planar
[Coordinates] =  cartesian
[ACSolver]    =  0
[PrevType]    =  0
[PrevSoln]    =  ""
[Comment]     =  ""
[PointProps]  = 0
[BdryProps]   = 1
  <BeginBdry>
    <BdryName> = "A=0"
    <BdryType> = 0
    <A_0> = 0
    <A_1> = 0
    <A_2> = 0
    <Phi> = 0
    <c0> = 0
    <c0i> = 0
    <c1> = 0
    <c1i> = 0
    <Mu_ssd> = 0
    <Sigma_ssd> = 0
    <innerangle> = 0
    <outerangle> = 0
  <EndBdry>
[BlockProps]  = 1
  <BeginBlock>
    <BlockName> = "Silicon Core Iron"
    <Mu_x> = 7000
    <Mu_y> = 7000
    <H_c> = 0
    <H_cAngle> = 0
    <J_re> = 0
    <J_im> = 0
    <Sigma> = 0
    <d_lam> = 0
    <Phi_h> = 0
    <Phi_hx> = 0
    <Phi_hy> = 0
    <LamType> = 0
    <LamFill> = 1
    <NStrands> = 0
    <WireD> = 0
    <BHPoints> = 0
  <EndBlock>
  [CircuitProps]  = 1
  <BeginCircuit>
    <CircuitName> = "A"
    <TotalAmps_re> = 1
    <TotalAmps_im> = 0
    <CircuitType> = 1
  <EndCircuit>
  [NumPoints] = 3
64.999999999999986	0	0	2
12.5	0	0	1
42.250000000000007	1.9000000000000004	0	0
[Solution]
4
0	151.6746	0	-1
2.6470867654124621	151.65149922376801	0	-1
-2.6470867654124017	151.65149922376759	0	-1
1.2931770755502037	148.18341546539696	0	-1
2
3	1	0	0	-1	5	-1	0
0	2	3	0	5	-1	-1	0

    """

    res = ReadFemm(string=data)
    print(res)

    from BasicTools.Helpers.Tests import WriteTempFile
    newFileName = WriteTempFile(filename="AnsFileTest.ans",content=data)
    import BasicTools.Containers.UnstructuredMesh as UM
    mesh = ReadFemm(fileName=newFileName)
    assert( mesh.GetNumberOfElements() == 4)
    assert( mesh.GetNumberOfNodes() == 4)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover

