# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import Tag as Tag
from OTTools.IO.WriterBase import WriterBase as WriterBase
import numpy as np
import OTTools.FE.ElementNames as EN

gmshName = {}
gmshName[EN.Bar_2]         = '1'
gmshName[EN.Triangle_3]    = '2'
gmshName[EN.Quadrangle_4]  = '3'
gmshName[EN.Tetrahedron_4] = '4'
gmshName[EN.Hexaedron_8]   = '5'
gmshName[EN.Wedge_6]       = '6'
gmshName[EN.Pyramid_5]     = '7'
gmshName[EN.Point_1]   = '15'



def WriteMeshToGmsh(filename,mesh, useOriginalId=False):    
    OW = GmshWriter()
    OW.Open(filename)
    OW.Write(mesh,useOriginalId = useOriginalId)
    OW.Close()
    
class GmshWriter(WriterBase):
    def __init__(self):
        super(GmshWriter,self).__init__()
    
    def __str__(self):
        res  = 'GmshWriter : \n'
        res += '   FileName : '+self.fileName+'\n'
        
    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Write(self,meshObject,useOriginalId=False):
        
        self.filePointer.write("$MeshFormat\n");
        self.filePointer.write("2.2 0 8\n");
        self.filePointer.write("$EndMeshFormat\n");
        self.filePointer.write("$Nodes\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} \n".format(numberofpoints) )
        #
        posn = meshObject.GetPosOfNode()
        if useOriginalId:
           for n in xrange(numberofpoints):
               self.filePointer.write("{} ".format(int(meshObject.originalIDNodes[n])))
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] ) 
        else:
           for n in xrange(numberofpoints):
               self.filePointer.write("{} ".format(n+1) )
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] ) 
        self.filePointer.write("$EndNodes\n");
        self.filePointer.write("$Elements\n");
        self.filePointer.write("{}\n".format(meshObject.GetNumberOfElements()))

        cpt = 1
        for ntype, data in meshObject.elements.iteritems(): 
	  if ntype == 'tet4':
            elemtype = gmshName[ntype]
            for i in xrange(data.GetNumberOfElements() ):
                if useOriginalId:
                    self.filePointer.write("{} {} {} {} {} ".format(data.originalIds[i],elemtype,2,1,1) )
                    self.filePointer.write(" ".join([str(meshObject.originalIDNodes[x]) for x in data.connectivity[i,:].flatten()]))
                else:
                    self.filePointer.write("{} {} {} {} {} ".format(data.originalIds[i],elemtype,2,1,1) )
                    self.filePointer.write(" ".join([str(x+1) for x in data.connectivity[i,:].flatten()]))
                cpt += 1
                self.filePointer.write("\n")
        
        #meshObject.PrepareForOutput();

        celtags = meshObject.GetNamesOfCellTags()
        elements = meshObject.elements.keys()
        tagcounter = 2
        for tagname in celtags:
          for elementContainer in elements:
	    elemtype = gmshName[elementContainer]
            try:
              for connectivity in meshObject.elements[elementContainer].connectivity[meshObject.elements[elementContainer].tags[tagname].id-1]:
                if useOriginalId:
                  print "yet to implement"
                else:
                  self.filePointer.write("{} {} {} {} {} ".format(cpt,elemtype,2,tagcounter,tagcounter) )
                  self.filePointer.write(" ".join([str(x+1) for x in connectivity]))
                cpt += 1
                self.filePointer.write("\n")
            except KeyError:
	      continue
	  tagcounter += 1
	    
        self.filePointer.write("$EndElements\n")
         
def CheckIntegrity():
    import OTTools.FE.UnstructuredMesh as UM
 
    from OTTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=np.int)
    
    tag = Tag("coucou")
    tag.AddToTag(0)
    mymesh.nodesTags["coucou"] = tag

    tris = mymesh.GetElementsOfType('tri3')
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)
    
    mymesh.AddElementToTagUsingOriginalId(3,"Tag1") 
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3") 
    
    OW = GmshWriter()
    OW.Open(tempdir+"Test_GmshWriter.geof")
    OW.Write(mymesh, useOriginalId=False)
    OW.Close()
    return "ok"
    
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover  
