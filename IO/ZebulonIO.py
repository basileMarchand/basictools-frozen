# -*- coding: utf-8 -*-

# Read Zebulon sparse matrices and vectors
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import time

def ReadMat(fileName, symetry=True , returnReorderOnly=False):

    if returnReorderOnly==True:

        f = open(fileName)
        with open(fileName, "r") as f:
            for line in f:
                try:
                    splittedLine = line.split()
                    if splittedLine[0]=='order':
                        numReorder = int(splittedLine[1])
                        numReorderLines = int(numReorder/10.)+4
                    if splittedLine[0]=='column_pointer':
                        break
                except IndexError:# pragma: no cover
                    continue
        f.close()
        reorder = []
        f = open(fileName)
        for i, line in enumerate(f):
            if i==1:
                sizeMat = int(line.split()[0])
            elif i > 2 and i < numReorderLines:
                splitted = line.split()
                for sp in splitted:
                    reorder.append(int(sp))
            elif i > numReorderLines:
                break
        f.close()

        return reorder

    else:

        start = time.time()
        skippedLines = [0,0,0,0]
        index = 0
        with open(fileName, "r") as f:
            for line in f:
                try:
                    splittedLine = line.split()
                    if splittedLine[0]=='order':
                        numReorder = int(splittedLine[1])
                        numReorderLines = int(numReorder/10.)+4+skippedLines[index]
                    if splittedLine[0]=='column_pointer':
                        index = 1
                        numColPtr = int(splittedLine[1])
                        numColPtrLines = int(numColPtr/10.)+numReorderLines+2+skippedLines[index]
                    if splittedLine[0]=='not_null':
                        index = 2
                        numNotNull = int(splittedLine[1])
                        numNotNullLines = int(numNotNull/10.)+numColPtrLines+2+skippedLines[index]
                    if splittedLine[0]=='upper_part':
                        index = 3
                        numUpperPart = int(splittedLine[1])
                        numUpperPartLines = int(numUpperPart/10.)+numNotNullLines+2+skippedLines[index]
                except IndexError:
                    skippedLines[index] += 1
                    continue

        f.close()

        reorder = []
        col = []
        row = []
        data = []

        f = open(fileName)

        for i, line in enumerate(f):
            if i==1:
                sizeMat = int(line.split()[0])
            elif i > 2 and i < numReorderLines:
                splitted = line.split()
                for sp in splitted:
                    reorder.append(int(sp))
            elif i > numReorderLines and i < numColPtrLines:
                splitted = line.split()
                for sp in splitted:
                    col.append(int(sp))
            elif i > numColPtrLines and i < numNotNullLines:
                splitted = line.split()
                for sp in splitted:
                    row.append(int(sp))
            elif i > numNotNullLines  and i < numUpperPartLines:
                splitted = line.split()
                for sp in splitted:
                    data.append(float(sp))
            elif i > numUpperPartLines:
                break
        f.close()

        A = csr_matrix((data, row, col), shape=(sizeMat, sizeMat), dtype = float)

        if symetry == True:
          A = A + A.T
          A.setdiag(A.diagonal()/2)

        return A


def ReadVec(fileName, dtype = float):

    V = []
    with open(fileName, "r") as f:
        for line in f:
            V.append(dtype(float(line)))

    return np.array(V)


def WriteVec(data, fileName):
    with open(fileName, 'w') as myFile:
        for dat in data:
            myFile.write(str(dat)+'\n')
    myFile.close()


def ReadInp(fileName=None,string=None):
    from cStringIO import StringIO
    from collections import OrderedDict as OD


    if fileName is not None:
        f = open(fileName, 'r')
        string = f.read()
        f.close()

    string = StringIO(string)


    res = OD();
    pobj= [res]
    cobj = pobj[-1];

    plevel = [5]
    clevel = plevel[-1]

    for line in string:
        l = line.strip('\n').lstrip().rstrip()
        if len(l) == 0 : continue
        if l[0] == "%" : continue
        l = l.split("%")[0]
        #print(str(clevel)+"---------------------------" + l)
        if l.find("****return")>-1 :
            pobj = []
            pobj.append(res)
            cobj = pobj[-1]
            continue

        while l.find("*"*(clevel))>-1 :
            pobj.pop()
            pobj.pop()
            cobj = pobj[-1]
            plevel.pop()
            clevel = plevel[-1]


        #print(cobj)

        for sublevel in xrange(4,0,-1):
          if l.find("*"*(sublevel))>-1 :
            #print(sublevel),
            #print(cobj)
            l = l[sublevel:]
            data = l.split()
            #print(data)
            if not cobj.has_key("*"*(sublevel)):
                cobj["*"*(sublevel)] = OD()
            pobj.append(cobj["*"*(sublevel)])
            cobj = pobj[-1]
            cobj[data[0]] = OD()
            pobj.append(cobj[data[0]])
            cobj = pobj[-1]
            plevel.append(sublevel)
            clevel = plevel[-1]
            if len(data)> 1:
                cobj["header"]  = data[1:]
            else:
                cobj["header"]  = ['']
            clevel = sublevel
            break
        else :
             #print(l)
            # Treat the cdata (no stating start)
            if not cobj.has_key("CDATA"):
                cobj["CDATA"] = []
            cobj["CDATA"].append(l)


    return res['****']

def WriteInp(data,output= None):
    import sys
    if output is None:
        output = sys.stdout

    for key4, value4 in data.iteritems():
        output.write("****" + key4),
        output.write(" ".join(value4["header"])+"\n")
        for i in xrange(0,4):
          #print(i)
          #print(value4)
          if value4.has_key('*'*i):
            for key3, value3 in value4['*'*i].iteritems():
              output.write("   "+ "*"*i + key3+" ")
              output.write(" ".join(value3["header"])+"\n")

              if value3.has_key("CDATA"):
                 #print(" ".join(value3["CDATA"]))
                 output.write(" \n".join(map(str,value3["CDATA"]))+"\n")

              for j in xrange(2,0,-1):
                 if value3.has_key('*'*j):
                   for key2, value2 in value3['*'*j].iteritems():
                      output.write("      "+"*"*j + key2+" ")
                      output.write(" ".join(map(str,value2["header"]))+"\n")
                      if value2.has_key("CDATA"):
                          #print(" ".join(value3["CDATA"]))
                          output.write(" \n".join(map(str,value2["CDATA"]))+"\n")
                      for k in xrange(1,0,-1):
                          if value2.has_key('*'*k):
                             for key1, value1 in value2['*'*k].iteritems():
                                 output.write("         "+"*"*k + key1+" ")
                                 output.write(" ".join(map(str,value1["header"]))+"\n")
                                 if value1.has_key("CDATA"):
                                     output.write(" \n".join(map(str,value1["CDATA"]))+"\n")
        output.write('****return\n')

def CheckIntegrity():
    import BasicTools.IO.ZebulonIO as ZIO
    import BasicTools.Helpers.Tests as T
    import BasicTools.TestData as T2
    dataPath = T2.GetTestDataPath()
    ZIO.ReadMat(dataPath+'Zmatrix')
    ZIO.ReadMat(dataPath+'Zmatrix', True)
    vec = ZIO.ReadVec(dataPath+'Zvector')
    tempPath = T.TestTempDir.GetTempPath()
    ZIO.WriteVec(vec, tempPath+'Zvector')

    res = ReadInp(T2.GetTestDataPath() + 'calcul1.inp')
    res['simulate']['***']['solver']['*']['ratio']['header'][1] = 1e-4
    res['simulate']['***']['test']['**']['load']['header'][1] = 6
    WriteInp(res)

    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
