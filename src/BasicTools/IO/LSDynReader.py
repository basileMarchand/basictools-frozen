# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" LSDyna file reader

"""
import numpy as np

from BasicTools.IO.ReaderBase import ReaderBase

from BasicTools.NumpyDefs import PBasicIndexType



LSDynNumber = {
  4:'tet4' # 4-point tetraedron
}



def ReadLSDyn(fileName=None,string=None,out=None,printNotRead=True):
    reader = LSDynReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out,printNotRead=printNotRead)
    return reader.output


def LineToListNoQuote(text):
    return [s.strip() for s in text.split()]


def ListToNumber(list):
    l = len(list)
    val = list[-1]
    for i in range(l-2, -1, -1):
        if list[i] != val:
            break
    return i+2

class LSDynReader(ReaderBase):
  def __init__(self):
        super(LSDynReader,self).__init__()
        self.commentChar= "$"
        self.readFormat = 'r'

  def Read(self, fileName=None,string=None,out=None,printNotRead=True):
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
    FENames = {}

    oidToElementContainer = {}
    oidToLocalElementNumber = {}
    l = self.ReadCleanLine()
    while(True):

      #premature EOF
      if l is None:
          print("ERROR premature EOF: please check the integrity of your geof file") # pragma: no cover
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
              nametype = LSDynNumber[n]
          except KeyError:
              raise("Elements with "+str(n)+"vertices not compatible with reader")
          conn = [x for x in  map(int,s[2:6])]
          elements = res.GetElementsOfType(nametype)
          oid = int(s[0])
          cpt = elements.AddNewElement(conn,oid)
          oidToElementContainer[oid] = elements
          oidToLocalElementNumber[oid] = cpt
        continue

      if l.find("*NODE")>-1:
        dim     = 3
        res.originalIDNodes= np.empty((0,1),dtype=PBasicIndexType)
        originalIds = []
        nodes= []
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
        res.nodes = np.array(nodes,dtype=float)
        res.nodes.shape = (cpt,dim)
        res.originalIDNodes = np.array(originalIds,dtype=int)

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

    for k in range(len(elements.connectivity)):
      for j in range(len(elements.connectivity[k])):
        elements.connectivity[k][j] = filetointernalid[elements.connectivity[k][j]]

    res.PrepareForOutput()

    self.output = res
    return res


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".k",LSDynReader)


def CheckIntegrity():

    data = u"""
$# LS-DYNA Keyword file created by LS-PrePost(R) V4.5.0
*KEYWORD
*TITLE
*ELEMENT_SOLID
$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8
       1 3000001       1       2       3       4       4       4       4       4
*NODE
$#   nid               x               y               z      tc      rc
       1             0.0             0.0             0.0       0       0
       2             0.6             0.0             0.0       0       0
       3             0.6             0.7             0.0       0       0
       4             0.2             0.2             0.3       0       0
       5             0.3             0.0             0.1       0       0
*END
    """

    res = ReadLSDyn(string=data)
    print(res)

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
