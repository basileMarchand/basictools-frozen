# -*- coding: utf-8 -*-
import OTTools.FE.UnstructuredMesh  as UM
import numpy as np



gmshNumber = {}
#0d
gmshNumber['1'] = "bar2"
gmshNumber['2'] = "tri3"
gmshNumber['3'] = "quad4"
gmshNumber['4'] = "tet4"
gmshNumber['5'] = "hex8"
gmshNumber['6'] = "wed6"
gmshNumber['7'] = "pyr5"
gmshNumber['15'] = "point1"

#1d
#2d
#gmshNumber[] = "tri3"
#gmshNumber[] = "tri6"
#3d
#gmshNumber[] = "tet4"


def ReadGmsh(fileName=None,string=None):
    import shlex
    from cStringIO import StringIO

    if fileName is not None:
        f = open(fileName, 'r')
        string = f.read()
        f.close()
    
    string = StringIO(string)
    
    res = UM.UnstructuredMesh()
    
    filetointernalid = {}
    #filetointernalidElem =  {}
    for line in string:
        l = line.strip('\n').lstrip().rstrip()
        if len(l) == 0: continue

        if l.find("$MeshFormat")>-1 :     
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue     
                if l.find("$EndMeshFormat") > -1:
                    break
            continue
        
        if l.find("$Nodes")>-1 : 
            line = string.readline()
            l = line.strip('\n').lstrip().rstrip()
            
            nbNodes = int(l.split()[0])
            print("Reading "+str(nbNodes)+ " Nodes")
            res.nodes = np.empty((nbNodes,3))
            res.originalIDNodes= np.empty((nbNodes,),dtype=np.int)
            cpt =0;
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("$EndNodes") > -1:
                    break
                s = l.split()
                #print(s)
                #print(res.originalIDNodes)
                oid = int(s[0])
                filetointernalid[oid] = int(cpt)
                res.originalIDNodes[int(cpt)] = int(s[0])
                res.nodes[cpt,:] = map(float,s[1:])
                cpt +=1
            continue

        if l.find("$Elements")>-1 : 
            line = string.readline()
            l = line.strip('\n').lstrip().rstrip()
            
            nbElements = int(l.split()[0])
            print("Reading "+str(nbElements)+ " Elements")
            #res.nodes = np.empty((nbNodes,dim))
            #res.originalIDNodes= np.empty((nbNodes,))
            cpt =0;
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("$EndElements") > -1:
                    if nbElements != cpt:
                        print("File problem!! number of elements read not equal to the total number of elemetns")
                        print(nbElements)
                        print(cpt)
                    break
                s = l.split()
                
                oid = int(s[0])
                gmshElemType = s[1]
                nametype = gmshNumber[gmshElemType]
                
                ntags = int(s[2])
                
                conn = [filetointernalid[x] for x in  map(int,s[(ntags+3):]) ]
                elements = res.GetElementsOfType(nametype)    
                elements.AddNewElement(conn,oid)
                
                if ntags >=0 :
                    res.AddElementToTagUsingOriginalId(oid,"PhyTag"+str(s[3]))    
                if ntags >=1 :
                    res.AddElementToTagUsingOriginalId(oid,"GeoTag"+str(s[4]))    
                
                for n in range(2, ntags):
                    res.AddElementToTagUsingOriginalId(oid,"ExtraTag"+str(n-2)+"_"+str(s[3+n]))    
                    
                cpt +=1
            continue
        print("ignoring line : " + l )
    return res

data = """ 
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 30 0 0
2 30 0 75
3 30 -2.5 0
$EndNodes
$Elements
2
1 15 2 223 1 2
2 15 2 227 2 3
$EndElements
"""

 
if __name__ == '__main__':

    #res = ReadGmsh(string=data)
    res = ReadGmsh(fileName="Soudage.msh")
    print("----")
    print(res.nodes)
    print(res.originalIDNodes)
    print(res.GetElementsOfType('bar2').connectivity)
    
    from OTTools.IO.GeofWriter import GeofWriter
    OW = GeofWriter()
    OW.Open("Soudage.geof")
    OW.Write(res, useOriginalId=True)
    OW.Close()

    from OTTools.IO.XdmfWriter import WriteMeshToXdmf
    
    WriteMeshToXdmf("Soudage.xdmf",res)