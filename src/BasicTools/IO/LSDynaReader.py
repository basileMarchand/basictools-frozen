# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" LSDyna file reader.
    Documentation of the format:
    http://ftp.lstc.com/anonymous/outgoing/jday/manuals/DRAFT_Vol_I.pdf
"""
import numpy as np

from BasicTools.IO.ReaderBase import ReaderBase

from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType


LSDynaNumber = {
    4:'tet4' # 4-point tetraedron
}


def ReadLSDyna(fileName=None,string=None,out=None,printNotRead=True):
    """Function API for reading a LSDyna mesh file

    Parameters
    ----------
    fileName : str, optional
        name of the file to be read, by default None
    string : str, optional
        data to be read as a string instead of a file, by default None
    out : UnstructuredMesh, optional
        output unstructured mesh object containing reading result, by default None
    printNotRead : bool, optional
        if True, prints in console the lines dot understood by the reader, by default True

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    reader = LSDynaReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out,printNotRead=printNotRead)
    return reader.output


def LineToListNoQuote(text):
    return [s.strip() for s in text.split()]


def ListToNumber(list):
    """hack where the element type is defined by inspecting the number of identical
    columns at the end of the connectivity list.
    """
    l = len(list)
    val = list[-1]
    for i in range(l-2, -1, -1):
        if list[i] != val:
            break
    return i+2

class LSDynaReader(ReaderBase):
    """LSDyna Reader class
    """
    def __init__(self):
        super(LSDynaReader,self).__init__()
        self.commentChar= "$"
        self.readFormat = 'r'

    def Read(self, fileName=None,string=None,out=None,printNotRead=True):
        """Function that performs the reading of a LSDyna mesh file

        Parameters
        ----------
        fileName : str, optional
            name of the file to be read, by default None
        string : str, optional
            data to be read as a string instead of a file, by default None
        out : UnstructuredMesh, optional
            output unstructured mesh object containing reading result, by default None
        printNotRead : bool, optional
            if True, prints in console the lines dot understood by the reader, by default True

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        import BasicTools.Containers.UnstructuredMesh as UM
        if fileName is not None:
            self.SetFileName(fileName)

        if string is not None:
            self.SetStringToRead(string)

        self.StartReading()

        if out is None:
            res = UM.UnstructuredMesh()
        else:
            res = out

        filetointernalid = {}

        oidToElementContainer = {}
        oidToLocalElementNumber = {}
        l = self.ReadCleanLine()

        originalIds = []
        nodes= []

        nodeSetCounter = 0

        while(True):

            #premature EOF
            if l is None:
                print("ERROR premature EOF: please check the integrity of your .k file") # pragma: no cover
                break # pragma: no cover
            #if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue

            if l.find("*ELEMENT")>-1:
                while(True):
                    l  = self.ReadCleanLine()
                    if l.find("*") > -1:
                        break
                    s = LineToListNoQuote(l)
                    n = ListToNumber(s[2:])
                    try:
                        nametype = LSDynaNumber[n]
                    except KeyError:
                        raise("Elements with "+str(n)+"vertices not compatible with reader")
                    conn = [x for x in  map(int,s[2:6])]
                    elements = res.GetElementsOfType(nametype)
                    oid = int(s[0])
                    cpt = elements.AddNewElement(conn,oid)
                    oidToElementContainer[oid] = elements
                    oidToLocalElementNumber[oid] = cpt

                    elTag = "canonical:"+s[1]
                    elements.tags.CreateTag(elTag,False).AddToTag(cpt-1)
                continue

            if l.find("*NODE")>-1:
                dim = 3
                s = None
                cpt = 0
                while(True):
                    l  = self.ReadCleanLine()
                    if l.find("*") > -1:
                        break
                    s = LineToListNoQuote(l)
                    oid = int(s[0])
                    filetointernalid[oid] = cpt
                    originalIds.append(oid)
                    nodes.append(list(map(float,s[1:dim+1])))
                    cpt += 1
                continue

            if l.find("*SET_NODE_LIST")>-1:
                tag = res.GetNodalTag(str(nodeSetCounter))
                nodeSetCounter += 1
                self.ReadCleanLine()
                while(True):
                    l  = self.ReadCleanLine()
                    if l.find("*") > -1:
                        break
                    s = np.array(l.split(), dtype=int)
                    for oid in s[s>0]:
                        tag.AddToTag(int(oid))
                continue

            if l.find("*END")>-1:
                self.PrintDebug("End file")
                break

            # case not treated
            if printNotRead == True:
                self.PrintDebug("line starting with <<"+l[:20]+">> not considered in the reader")
            l = self.ReadCleanLine()
            continue

        self.EndReading()

        res.nodes = np.array(nodes,dtype=PBasicFloatType)
        res.nodes.shape = (cpt,dim)
        res.originalIDNodes = np.array(originalIds,dtype=PBasicIndexType)

        updateIDsFunc = lambda x: filetointernalid[x]

        for tag in res.nodesTags:
            tag.SetIds(np.vectorize(updateIDsFunc)(tag.GetIds()))


        for _, data in res.elements.items():
            data.tighten()
            data.connectivity = np.vectorize(updateIDsFunc)(data.connectivity)

        res.PrepareForOutput()

        self.output = res
        return res


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".k",LSDynaReader)

def CheckIntegrity():

    data = u"""
$# LS-DYNA Keyword file created by LS-PrePost(R) V4.5.0
*KEYWORD
*TITLE
*ELEMENT_SOLID
$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
       1 3000001       1       2       3       4       4       4       4       4
*SET_NODE_LIST
$#     sid       da1       da2       da3       da4    solver
         1       0.0       0.0       0.0       0.0MECH
$#    nid1      nid2      nid3      nid4      nid5      nid6      nid7      nid8
         1         2         3         4         0         0         0         0
*SET_NODE_LIST
$#     sid       da1       da2       da3       da4    solver
         1       0.0       0.0       0.0       0.0MECH
$#    nid1      nid2      nid3      nid4      nid5      nid6      nid7      nid8
         1         3         0         0         0         0         0         0
*NODE
$#   nid               x               y               z      tc      rc
       1             0.0             0.0             0.0       0       0
       2             0.6             0.0             0.0       0       0
       3             0.6             0.7             0.0       0       0
       4             0.2             0.2             0.3       0       0
       5             0.3             0.0             0.1       0       0
*END
"""

    res = ReadLSDyna(string=data)
    print(res)

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
