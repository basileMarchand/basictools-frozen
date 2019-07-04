# -*- coding: utf-8 -*-

""" Geo fiel reader (Zset mesh file)

"""

import numpy as np
import struct

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
import BasicTools.Containers.UnstructuredMesh as UM
from BasicTools.IO.ZsetTools import GeofNumber,PermutationZSetToBasicTools, nbIntegrationsPoints

def ReadGeo(fileName=None,out=None,readElset=True,readFaset=True):
    reader = GeoReader()
    reader.SetFileName(fileName)
    reader.Read(fileName=fileName, out=out,readElset=readElset,readFaset=readFaset)
    return reader.output

class GeoReader(ReaderBase):
   def __init__(self):
        super(GeoReader,self).__init__()
        self.readFormat = 'rb'

   def _getTag(self):
       nb = struct.unpack(">i", self.rawread(4))[0]
       rawtag = self.rawread(nb)
       return rawtag[0:-1].decode("utf-8")

   def readMagic(self):
           magic = self.rawread(13)
           if not (b'Z7BINARYGEOF\x00' == magic ) :
               raise(Exception("Bad file"))

   def readInt32(self):
       return struct.unpack(">i", self.rawread(4))[0]

   def readInts32(self,cpt):
       return np.array(struct.unpack(">"+"i"*cpt, self.rawread(4*cpt)))

   def ReadMetaData(self):
       return self.Read(onlyMeta = True)

   def Read(self, fileName=None,out=None,readElset=True,readFaset=True, onlyMeta = False):

       if fileName is not None:
          self.SetFileName(fileName)

       self.StartReading()

       res = UM.UnstructuredMesh()
       metadata = {}

       oidToElementContainer = {}
       oidToLocalElementNumber = {}

       cpt = 0
       while True:
           self.readMagic()
           tag = self._getTag()
           if tag == u"node":
               nbNodes = self.readInt32()
               metadata["nbNodes"] = int(nbNodes)
               dims = self.readInt32()

               if onlyMeta :
                   self.filePointer.seek((dims*8+4)*nbNodes,1 )
               else:
                   res.nodes = np.zeros((nbNodes,3))
                   res.originalIDNodes = np.empty((nbNodes,),dtype=np.int)

                   for i in range(nbNodes):
                       data = self.rawread((dims*8+4))
                       data = struct.unpack((">i"+"d"*dims), data)
                       res.originalIDNodes[i] = data[0]
                       res.nodes[i,0:dims] = data[1:]

           elif tag == u"element":
               n_elem = self.readInt32()
               metadata['nbElements'] = int(n_elem)
               IPPerElement = np.empty(metadata['nbElements'],dtype= np.int)

               n_grp = self.readInt32()
               for i in range(n_grp):
                   n_in_grp = self.readInt32()
                   ltype =  self._getTag()
                   nametype = GeofNumber[ltype]
                   elements = res.GetElementsOfType(nametype)
                   elements.Allocate(n_in_grp)
                   n_node = self.readInt32()

                   nintegpoints =  nbIntegrationsPoints[ltype]
                   if onlyMeta :
                       self.filePointer.seek((4)*(n_node+2)*n_in_grp,1 )
                       for j in range(n_in_grp):
                           IPPerElement[cpt] = nintegpoints
                           cpt += 1
                   else:
                       perm = None
                       if ltype in PermutationZSetToBasicTools:
                           perm = PermutationZSetToBasicTools[ltype]


                       for j in range(n_in_grp):
                           idd = self.readInt32()
                           rank = self.readInt32()
                           conn = self.readInts32(n_node)
                           if perm:
                               conn =  [conn[x] for x in perm ]
                           elements.connectivity[j,:] = conn
                           elements.originalIds[j] = rank

                           oidToElementContainer[rank] = elements
                           oidToLocalElementNumber[rank] = j
                           IPPerElement[cpt] = nintegpoints
                           cpt += 1

           elif tag == u"nset":
               n_nset = self.readInt32()
               for i in range(n_nset):
                   name = self._getTag()
                   size = self.readInt32()

                   if  onlyMeta :
                       self.filePointer.seek(size*4,1)
                   else:
                       ids = np.fromfile(self.filePointer,count=size, dtype=np.int32).byteswap()
                       res.nodesTags.CreateTag(name).SetIds(ids)

           elif tag == u"elset":
               n_nset = self.readInt32()
               for i in range(n_nset):
                   name = self._getTag()

                   size = self.readInt32()

                   if onlyMeta :
                       self.filePointer.seek(size*4,1)
                   else:
                       ids = np.fromfile(self.filePointer,count=size, dtype=np.int32).byteswap()
                       for oid in ids:
                           oidToElementContainer[oid].tags.CreateTag(name,False).AddToTag(oidToLocalElementNumber[oid])

           elif tag == u"bset":
               n_bset = self.readInt32()
               for i in range(n_bset):
                   fasetName = self._getTag()
                   eltype = self._getTag()
                   n_el = self.readInt32()

                   if eltype == "faset":
                       if onlyMeta :
                           for j in range(n_el):
                               how = self.readInt32()
                               self.filePointer.seek(how*4+5,1)
                       else:
                           for j in range(n_el):
                               how = self.readInt32()
                               conn = np.fromfile(self.filePointer,count=how, dtype=np.int32).byteswap()
                               rawname = self.rawread(5)
                               bsetelemtype = rawname[0:2].decode("utf-8")

                               nametype = GeofNumber[bsetelemtype]

                               perm = None
                               if bsetelemtype in PermutationZSetToBasicTools:
                                   conn =  [conn[x] for x in PermutationZSetToBasicTools[bsetelemtype] ]

                               elements = res.GetElementsOfType(nametype)
                               localId = elements.AddNewElement(conn,-1)
                               elements.GetTag(fasetName).AddToTag(localId-1)

                   else:
                       raise(Exception("Error I dont know how to treat bset of type "   + str(eltype) ))

           elif tag == "ipset":
               d = self.readInt32()
               for i in range(d):
                   name = self.readTag()
                   l = self.readInt32()
                   for j in range(l):
                       nb = self.readInt32()
                       nb = self.readInt32()
                       self.filePointer.seek(4*nb)
           elif tag == "return":
               if onlyMeta:
                   metadata['nbIntegrationPoints'] = np.sum(IPPerElement)
                   metadata['IPPerElement'] = IPPerElement
                   return metadata
               else:
                   self.output  = res
                   return res
           else:
               print(res)
               raise(Exception("Tag '"+ str(tag )+ "' not treated"))

       self.EndReading()







from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".geo",GeoReader)


def CheckIntegrity():

    from BasicTools.TestData import GetTestDataPath
    fileName  = GetTestDataPath()+"cube.geo"

    reader = GeoReader()
    reader.SetFileName(fileName)
    print(reader.ReadMetaData())
    res = ReadGeo(fileName=fileName)
    print(res)
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
