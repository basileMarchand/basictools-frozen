# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import Tag as Tag
from OTTools.IO.WriterBase import WriterBase as WriterBase
import numpy as np
import OTTools.FE.ElementNames as EN

GeofName = {}
#0d
GeofName[EN.Point_1] = "l2d1"

#1d
GeofName[EN.Bar_2] = "l2d2"
#2d
GeofName[EN.Triangle_3] = "c2d3"
GeofName[EN.Triangle_6] = "c2d6"
#3d
GeofName[EN.Tetrahedron_10] = "c3d10"
GeofName[EN.Tetrahedron_4] = "c3d4"
GeofName[EN.Quadrangle_4] = "c2d4"
GeofName[EN.Hexaedron_8] = "c3d8"

GeofName[EN.Wedge_6] = "c3d6"

def WriteMeshToGeof(filename,mesh, useOriginalId=False):    
    OW = GeofWriter()
    OW.Open(filename)
    OW.Write(mesh,useOriginalId = useOriginalId)
    OW.Close()
    
class GeofWriter(WriterBase):
    def __init__(self):
        super(GeofWriter,self).__init__()
    
    def __str__(self):
        res  = 'GeofWriter : \n'
        res += '   FileName : '+self.fileName+'\n'
        
    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Write(self,meshObject,useOriginalId=False):
        
        self.filePointer.write("% This file has been writen by the python routine GeofWriter of the OTTools package\n") 
        self.filePointer.write("% For any question about this routine, please contact SAFRAN TECH Pole M&S Team OT\n") 

        self.filePointer.write("***geometry\n");
        self.filePointer.write("  **node\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} {} \n".format(numberofpoints,3) )
        #
        posn = meshObject.GetPosOfNode()
        if useOriginalId:
           for n in xrange(numberofpoints):
               self.filePointer.write("{} ".format(int(meshObject.originalIDNodes[n])))
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] ) 
        else:
           for n in xrange(numberofpoints):
               self.filePointer.write("{} ".format(n) )
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] ) 
        #
        self.filePointer.write("  **element\n")
        self.filePointer.write("{}\n".format(meshObject.GetNumberOfElements()))

        cpt =0;
        for ntype, data in meshObject.elements.iteritems(): 
            elemtype = GeofName[ntype]
            #npe = data.GetNumberOfNodesPerElement()
            for i in xrange(data.GetNumberOfElements() ):
                if useOriginalId:
                    self.filePointer.write("{} {} ".format(data.originalIds[i],elemtype) )
                    self.filePointer.write(" ".join([str(meshObject.originalIDNodes[x]) for x in data.connectivity[i,:].flatten()]))
                else:
                    self.filePointer.write("{} {} ".format(cpt,elemtype) )
                    self.filePointer.write(" ".join([str(x) for x in data.connectivity[i,:].flatten()]))
                cpt += 1;
                self.filePointer.write("\n")

        self.filePointer.write(" ***group \n")
                
        for tagname in meshObject.nodesTags:
            self.filePointer.write("  **nset {} \n".format(tagname))
            data = np.zeros((meshObject.GetNumberOfNodes(),1),dtype=np.int)
            if useOriginalId:
                self.filePointer.write(" ".join([str(meshObject.originalIDNodes[x]) for x in meshObject.nodesTags[tagname].id]))
            else:
                self.filePointer.write(" ".join([str(x) for x in meshObject.nodesTags[tagname].id]))
            self.filePointer.write("\n")
        
        meshObject.PrepareForOutput();

        celtags = meshObject.GetNamesOfCellTags()
        for tagname in celtags:
            self.filePointer.write("  **elset {} \n".format(tagname))
            data = meshObject.GetElementsInTag(tagname,useOriginalId=useOriginalId)
            self.filePointer.write(" ".join([str(x) for x in data]))
            self.filePointer.write("\n")
         


        self.filePointer.write("****return \n")

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
    
    OW = GeofWriter()
    OW.Open(tempdir+"Test_GeoWriter.geof")
    OW.Write(mymesh, useOriginalId=True)
    OW.Close()
    return "ok"
    
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   
