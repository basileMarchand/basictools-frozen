# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import numpy as np
import OTTools.FE.UnstructuredMesh as UM
import OTTools.FE.ElementNames as EN

AscNumber = {}

AscNumber['2006'] = EN.Triangle_6
AscNumber['3010'] = EN.Tetrahedron_10
AscNumber['1002'] = EN.Bar_2

    
def ReadAsc(fileName=None,string=None):
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

        if l.find("BEGIN_NODES")>-1 : 
            
            nbNodes = int(l.split()[1])
            print("Reading "+str(nbNodes)+ " Nodes")
            dim = int(l.split()[2])
            res.nodes = np.empty((nbNodes,dim))
            res.originalIDNodes= np.empty((nbNodes,))
            cpt =0;
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("END_NODES") > -1:
                    break
                s = l.split()
                #print(s)
                #print(res.originalIDNodes)
                oid = int(s[0])
                filetointernalid[oid] = cpt
                res.originalIDNodes[cpt] = int(s[0])
                res.nodes[cpt,:] = map(float,s[6:])
                cpt +=1
            continue

        if l.find("BEGIN_ELEMENTS")>-1 : 
            
            nbElements = int(l.split()[1])
            print("Reading "+str(nbElements)+ " Elements")
            #res.nodes = np.empty((nbNodes,dim))
            #res.originalIDNodes= np.empty((nbNodes,))
            cpt =0;
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("END_ELEMENTS") > -1:
                    if nbElements != cpt:
                        print("File problem!! number of elements read not equal to the total number of elemetns")
                        print(nbElements)
                        print(cpt)
                    break
                s = l.split()
                
                nametype = AscNumber[s[1]];
                
                conn = [filetointernalid[x] for x in  map(int,s[5:]) ]
                
                # for some types we need permutation
                if nametype == EN.Triangle_6:
                    conn = [ conn[per] for per in [0, 2, 4, 1, 3, 5] ]
                elif nametype == EN.Tetrahedron_10:
                    conn = [ conn[per] for per in [0, 2, 4, 9, 1,3 ,5,6,7,8] ]
                    
                elements = res.GetElementsOfType(nametype)    
                oid = int(s[0])

                elements.AddNewElement(conn,oid)
                cpt +=1
            continue
        if l.find("BEGIN_GROUPS")>-1 : 
            
            nbgroups = int(l.split()[1])
            print("Reading "+str(nbgroups)+ " Groups")
            #res.nodes = np.empty((nbNodes,dim))
            #res.originalIDNodes= np.empty((nbNodes,))
            cpt =0;
            while(True):
                
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("END_GROUPS") > -1:
                    if(nbgroups != (cpt)):
                        print("File problem!! number of groups read not equal to the total number of groups")
                        print(nbgroups)
                        print(cpt)
                    break
                s = shlex.split(l)
                print("Reading Group " + s[1]) 
                if s[2] == '1' :
                    #node group
                    tag = res.GetNodalTag(s[1])
                    tag.id = np.array( [filetointernalid[x] for x in  map(int,s[7:]) ] ,dtype=np.int)
                else:
                    #element group
                    for x in range(7,len(s)) :
                        Oid = int(s[x])
                        #print(Oid)
                        res.AddElementToTagUsingOriginalId(Oid,s[1])
                cpt +=1
            continue
        print("Ignoring line : " + l )
    return res
    


def CheckIntegrity():
        
    __checkintegritydata = """   
BEGIN_NODES 6 3
1 0 0 0 0 0 295.673175860532 28.0704731415872 346.109138100075
2 0 0 0 0 0 295.105225 28.260575 345.628395
3 0 0 0 0 0 295.180501114015 25.8084581250318 344.876373186428
4 0 0 0 0 0 295.3886425 28.1693925 345.8617875
5 0 0 0 0 0 295.231426751792 27.0671220355891 345.153077365196
6 0 0 0 0 0 295.629817604236 26.9160040797623 345.45435766729
END_NODES
BEGIN_ELEMENTS 1
21 2006 0 0 0  1 2 3 4 5 6
END_ELEMENTS
BEGIN_GROUPS 1
2 M2D 2 0 "PART_ID 2"  ""  "PART built in Visual-Environment" 21
END_GROUPS
"""
    res = ReadAsc(string=__checkintegritydata)
    print(res)
    if res.GetElementsOfType(EN.Triangle_6).originalIds[0] != 21:
        raise
    return 'ok'
 
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   
