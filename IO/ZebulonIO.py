# -*- coding: utf-8 -*-
""" Read Zebulon sparse matrices and vectors
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import time
from collections import OrderedDict as OD

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
    from collections import OrderedDict as OD


    if fileName is not None:
        string = open(fileName, 'r')
    elif string is not None:
        from io import StringIO
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
            try:
              pobj.pop()
              pobj.pop()
              cobj = pobj[-1]
              plevel.pop()
              clevel = plevel[-1]
            except IndexError:
              break

        for sublevel in range(4,0,-1):
          if len(l)>=sublevel and l[:sublevel]=="*"*(sublevel):
          #if l.find("*"*(sublevel))>-1 :
            l = l[sublevel:]
            data = l.split()
            #print(data)
            if "*"*(sublevel) not in cobj:
                cobj["*"*(sublevel)] = OD()

            pobj.append(cobj["*"*(sublevel)])
            cobj = pobj[-1]
            #if data[0] not in cobj:
            cobj[data[0]] = OD()
            
            pobj.append(cobj[data[0]])
            cobj = pobj[-1]
            plevel.append(sublevel)
            clevel = plevel[-1]
            if len(data)> 1:
              cobj["header"]  = data[1:]
              #cobj[data[1]]  = data[1:]
            else:
              cobj["header"]  = ['']
            clevel = sublevel
            break
        else:
            #print(l)
            # Treat the cdata (no stating start)
            if "CDATA" not in cobj:
                cobj["CDATA"] = []
            cobj["CDATA"].append(l)

    return res['****']



def WriteInp(data,output= None):
    import sys
    if output is None:
        output = sys.stdout

    for key4, value4 in data.items():
        output.write("****" + key4),
        output.write(" ".join(value4["header"])+"\n")
        for i in range(0,4):
          #print(i)
          #print(value4)
          if '*'*i in value4:
            for key3, value3 in value4['*'*i].items():
              output.write("   "+ "*"*i + key3+" ")
              output.write(" ".join(value3["header"])+"\n")

              if "CDATA" in value3:
                 #print(" ".join(value3["CDATA"]))
                 output.write(" \n".join(map(str,value3["CDATA"]))+"\n")

              for j in range(2,0,-1):
                 if '*'*j in value3:
                   for key2, value2 in value3['*'*j].items():
                      output.write("      "+"*"*j + key2+" ")
                      output.write(" ".join(map(str,value2["header"]))+"\n")
                      if "CDATA" in value2:
                          #print(" ".join(value3["CDATA"]))
                          output.write(" \n".join(map(str,value2["CDATA"]))+"\n")
                      for k in range(1,0,-1):
                          if '*'*k in value2:
                             for key1, value1 in value2['*'*k].items():
                                 output.write("         "+"*"*k + key1+" ")
                                 output.write(" ".join(map(str,value1["header"]))+"\n")
                                 if "CDATA" in value1:
                                     output.write(" \n".join(map(str,value1["CDATA"]))+"\n")
        output.write('****return\n')



def ReadInp2(fileName=None,string=None):
    """
    Reads an Zebulon inp files and return the results in nested lists
    """

    if fileName is not None:
      string = open(fileName, 'r')
    elif string is not None:
      from io import StringIO
      string = StringIO(string)


    string0 = ''
    for line in string:
      l = line.strip('\n').lstrip().rstrip()
      if (len(l) == 0) or (l[0] == "%"): continue
      l = l.split("%")[0]
      for ll in l:
        string0 += ll
      string0 += '\n'

    string = string0.strip('\n').lstrip().rstrip().split('****return')[:-1]
    for i in range(len(string)):
      string[i] = string[i].strip('\n').lstrip().rstrip()[4:].strip('\n').lstrip().rstrip().split('\n***')
      for j in range(len(string[i])):
        string[i][j] = string[i][j].strip('\n').lstrip().rstrip().split('\n**')
        for k in range(len(string[i][j])):
          string[i][j][k] = string[i][j][k].strip('\n').lstrip().rstrip().split('\n*')
          for l in range(len(string[i][j][k])):
            string[i][j][k][l] = string[i][j][k][l].strip('\n').lstrip().rstrip().split('\n')
            for m in range(len(string[i][j][k][l])):
              string[i][j][k][l][m] = string[i][j][k][l][m].lstrip().rstrip().split()
    return string


def WriteInp2(data,output= None):
    """
    Writes an Zebulon inp files constructed by the function ReadInp2
    """
    import sys
    if output is None:
        output = sys.stdout

    current4star = data[0][0][0][0]
    for i in range(len(data)):
      for j in range(len(data[i])):
        for k in range(len(data[i][j])):
          for l in range(len(data[i][j][k])):
            if data[i][0][0][0] != current4star:
              output.write('****return\n\n')
              current4star = data[i][0][0][0]
            output.write(Nstar(i, j, k, l)+ValueToString(data[i][j][k][l])+'\n')

    output.write('****return\n')
            



def GetFromInp(data,dic):
    """
    returns a list of lists containing the elements of each line of the inp file "data" read by ReadInp2()
    respecting the conditions in dic: for example, for dic = {'4':'simulate', '2':'model'}, the function returns
    all the lines being in a **** section starting by "simulate" and a ** section starting by "model" 
    """
    res = []

    for i in range(len(data)):
        for j in range(len(data[i])):
            for k in range(len(data[i][j])):
                for l in range(len(data[i][j][k])):
                  count = 0
                  if '4' in dic and data[i][0][0][0][0][0] == dic['4']:
                    count += 1
                  if '3' in dic and data[i][j][0][0][0][0] == dic['3']:
                    count += 1
                  if '2' in dic and data[i][j][k][0][0][0] == dic['2']:
                    count += 1
                  if '1' in dic and data[i][j][k][l][0][0] == dic['1']:
                    count += 1
                  if count == len(dic):
                    res.append(data[i][j][k][l])
    return res


def Nstar(i, j, k, l):
    if l!=0:
      return "      *"
    if k!=0:
      return "    **"
    if j!=0:
      return "  ***"
    return "****"


def NstarDigit(i, j, k, l):
    if l!=0:
      return 1
    if k!=0:
      return 2
    if j!=0:
      return 3
    return 4


def ValueToString(value):
    string = ''
    for i in range(len(value)):
      for j in range(len(value[i])):
        string += value[i][j]+' '
      string = string[:-1]+'\n'
    return string[:-1]


def GetCycleTables(data):
    """
    return a dictionary containing the infos of an inp file "data" read by ReadInp2()
    concerning all the infos respecting ****calcul, ***table and **cycle
    """
    tableData = GetFromInp(data,{'4':'calcul', '3':'table', '2':'cycle'})
    
    cycleTables = {}
    for t in range(int(len(tableData)/3)):
      cycleTables[tableData[3*t][0][1]] = {}
      cycleTables[tableData[3*t][0][1]]['tInit'] = tableData[3*t][0][2]
      cycleTables[tableData[3*t][0][1]]['tEnd'] = tableData[3*t][0][3]
      tempTime = []
      for i in range(1, len(tableData[3*t+1])):
        tempTime += tableData[3*t+1][i]
      tempValue = []
      for i in range(1, len(tableData[3*t+2])):
        tempValue += tableData[3*t+2][i]
      cycleTables[tableData[3*t][0][1]]['time'] = np.array(tempTime, dtype = float)
      cycleTables[tableData[3*t][0][1]]['value'] = np.array(tempValue, dtype = float)

    return cycleTables


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

    # check ReadInp and WriteInp
    res = ReadInp(T2.GetTestDataPath() + 'calcul1.inp')
    res['simulate']['***']['solver']['*']['ratio']['header'][1] = 1e-4
    res['simulate']['***']['test']['**']['load']['header'][1] = 6
    WriteInp(res)

    # check ReadInp2 and WriteInp2
    res = ReadInp2(T2.GetTestDataPath() + 'calcul2.inp')
    GetFromInp(res,{'4':'calcul', '2':'impose_nodal_dof'})
    GetFromInp(res,{'4':'calcul', '3':'linear_solver'})
    GetCycleTables(res)
    WriteInp2(res)

    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
