# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import Tag as Tag
from OTTools.IO.WriterBase import WriterBase as WB
import numpy as np


GeofName = {}
#0d
GeofName["bar2"] = "line"
#1d
#2d
GeofName["tri"] = "t3"
GeofName["tri6"] = "t6"
#3d
GeofName["tet10"] = "c3d10"

def WriteMeshToGeof(filename,mesh, useOriginalId=False):    
    OW = GeofWriter()
    OW.Open(filename)
    OW.Write(mesh,useOriginalId = useOriginalId)
    OW.Close()
    
class GeofWriter(WB):
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
               self.filePointer.write("{} ".format(meshObject.originalIDNodes[n]  ) )
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
        self.filePointer.write("  **elset name \n")
        
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

if __name__ == '__main__':
    import OTTools.FE.UnstructuredMesh as UM
    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=np.int)
    
    tag = Tag("coucou")
    tag.AddToTag(0)
    mymesh.nodesTags["coucou"] = tag

    tris = mymesh.GetElementsOfType('tri')
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)
    
    mymesh.AddElementToTagUsingOriginalId(3,"Tag1") 
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3") 
    
    OW = GeofWriter()
    OW.Open("TestOutputGeo.geof")
    OW.Write(mymesh, useOriginalId=True)
    OW.Close()
    
    #import OTTools.IO.AscReader as AR
    #m =  AR.ReadAsc(fileName='C:\\Users\\D584808\\Documents\\Projects\\Python\\Topotools\\SUPPORT_VERIN_DATA1.ASC')
    #WriteMeshToGeof('FromASCReader.geof',m, useOriginalId=True)