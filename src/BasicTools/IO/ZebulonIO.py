# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" Read Zebulon sparse matrices and vectors
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import time
from collections import OrderedDict as OD
import os

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
    """
    dictionary-based reader: WARNING: multiple entry keys may be overwritten !!
    """
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
    """
    dictionary-based reader: WARNING: multiple entry keys may be missing !!
    """
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


def StringReader(inputString, folder):
    """
    convert string from Zebulon input file into cleaned string, robust with respect to (potentially nested) '@include'
    """

    appendedString = ''

    for line in inputString:
      l = line.strip('\n').lstrip().rstrip()
      if (len(l) == 0) or (l[0] == "%"): continue
      l = l.split("%")[0]
      if l.split()[0]=='@include':
        sep = os.sep
        lsep = len(sep)
        if folder[-lsep:] == sep or len(folder) == 0:
          sep = ""
        appendedString += StringReader(open(folder + sep + l.split()[1], 'r'), folder)
        continue
      for ll in l:
        appendedString += ll
      appendedString += '\n'

    return appendedString


def ReadInp2(fileName=None,string=None,startingNstar=4):
    """
    Reads an Zebulon inp files and return the results in nested lists
    """

    if fileName is not None:
      string = open(fileName, 'r')
    elif string is not None:
      from io import StringIO
      string = StringIO(string)

    string0 = StringReader(inputString = string, folder = os.path.dirname(fileName))

    string = string0.strip('\n').lstrip().rstrip().split(startingNstar*"*"+'return')[:-1]
    for i in range(len(string)):
      string[i] = string[i].strip('\n').lstrip().rstrip()[startingNstar:].strip('\n').lstrip().rstrip().split('\n***')
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
    Returns a list of lists containing the elements of each line of the inp
    file "data" read by ReadInp2() respecting the conditions in dic: for
    example, for dic = {'4':['calcul'], '2':['file', 'temperature']}, the
    function returns all the lines being in a ``****`` section starting by
    "simulate" and a ``**`` section starting by "file temperature" (robust with
    respect to the number of spaces between 'file' and 'temperature')
    """

    res = []
    for i in range(len(data)):
        for j in range(len(data[i])):
            for k in range(len(data[i][j])):
                for l in range(len(data[i][j][k])):
                  count = 0
                  if '4' in dic and all(data[i][0][0][0][0][p] == val for p, val in enumerate(dic['4'])):
                    count += 1
                  if '3' in dic and all(data[i][j][0][0][0][p] == val for p, val in enumerate(dic['3'])):
                    count += 1
                  if '2' in dic and all(data[i][j][k][0][0][p] == val for p, val in enumerate(dic['2'])):
                    count += 1
                  if '1' in dic and all(data[i][j][k][l][0][p] == val for p, val in enumerate(dic['1'])):
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



def GetInputTimeSequence(data):
    """
    returns the list of the time steps to be computed (from Zebulon input file)
    """

    incrementSequence = []
    timeSequence      = []
    dTimeSequence     = []

    dataIncrement = GetFromInp(data,{'3':['resolution'], '1':['increment']})
    dataTime      = GetFromInp(data,{'3':['resolution'], '1':['time']})
    dataDtime     = GetFromInp(data,{'3':['resolution'], '1':['dtime']})

    if len(dataTime) > 0:
      count = 0
      for line in dataTime[0]:
        for t in line:
          if count > 0:
            timeSequence.append(float(t))
          count += 1

    if len(dataDtime) > 0:
      if len(dataTime) > 0:
        print("both time and dtime defined !!!")
        return
      count = 0
      for line in dataDtime[0]:
        for t in line:
          if count > 0:
            dTimeSequence.append(float(t))
          count += 1
      timeSequence.append(0)
      for i, dt in enumerate(dTimeSequence):
        timeSequence.append(timeSequence[i]+dt)

    if len(dataIncrement) > 0:
      count = 0
      for line in dataIncrement[0]:
        for t in line:
          if count > 0:
            incrementSequence.append(int(t))
          count += 1

    reconstructedTimeSequence = []

    lt = len(timeSequence)
    if timeSequence[0] != 0:
      reconstructedTimeSequence.append(0.)

    if lt == 1:
      for j in range(incrementSequence[0]):
        reconstructedTimeSequence.append((j+1)*(timeSequence[0])/incrementSequence[0])
    else:
      reconstructedTimeSequence.append(timeSequence[0])
      for i in range(lt-1):
        for j in range(incrementSequence[i]):
          #print(i, j, len(timeSequence), len(incrementSequence))
          reconstructedTimeSequence.append(timeSequence[i]+(j+1)*(timeSequence[i+1]-timeSequence[i])/incrementSequence[i])


    return reconstructedTimeSequence



def GetTables(data):
    """
    Returns a dictionary containing the infos of an inp file "data" read by
    ReadInp2() concerning all the infos respecting ``****calcul``, ``***table``
    and ``**name`` or ``**cycle``.
    """

    cyclicTableData = GetFromInp(data,{'4':['calcul'], '3':['table'], '2':['cycle']})
    noncyclicTableData = GetFromInp(data,{'4':['calcul'], '3':['table'], '2':['name']})

    tables = {}
    tableData = cyclicTableData
    for t in range(int(len(tableData)/3)):
      tables[tableData[3*t][0][1]] = {}
      tables[tableData[3*t][0][1]]['tInit'] = tableData[3*t][0][2]
      tables[tableData[3*t][0][1]]['tEnd']  = tableData[3*t][0][3]
      tempTime = []
      for i in range(1, len(tableData[3*t+1])):
        tempTime += tableData[3*t+1][i]
      tempValue = []
      for i in range(1, len(tableData[3*t+2])):
        tempValue += tableData[3*t+2][i]
      tables[tableData[3*t][0][1]]['time']  = np.array(tempTime, dtype = float)
      tables[tableData[3*t][0][1]]['value'] = np.array(tempValue, dtype = float)

      if len(tempTime) != len(tempValue):
        print("WARNING ! length of time table and value table for "+ tableData[3*t][0][1] + "are different")

    tableData = noncyclicTableData
    for t in range(int(len(tableData)/3)):
      tables[tableData[3*t][0][1]] = {}
      tempTime = []
      for i in range(1, len(tableData[3*t+1])):
        tempTime += tableData[3*t+1][i]
      tempValue = []
      for i in range(1, len(tableData[3*t+2])):
        tempValue += tableData[3*t+2][i]
      tables[tableData[3*t][0][1]]['time']  = np.array(tempTime, dtype = float)
      tables[tableData[3*t][0][1]]['value'] = np.array(tempValue, dtype = float)

      if len(tempTime) != len(tempValue):
        print("WARNING ! length of time table and value table for "+ tableData[3*t][0][1] + " are different")

    return tables




def GetBoundaryConditions(data):
    """
    Returns a dictionary containing the infos of an inp file "data" read by ReadInp2()
    concerning all the infos respecting ``****calcul``, ``***bc``
    """
    bcData = GetFromInp(data,{'4':['calcul'], '3':['bc']})
    bcs = {}
    for i in range(1, len(bcData)):
      if bcData[i][0][0] not in bcs:
        bcs[bcData[i][0][0]] = []
      bc = [b for b in bcData[i][1:]]
      bcs[bcData[i][0][0]].append(bc)

    return bcs



def GetLoadings(data):
    """
    returns a dictionary containing the infos of an inp file "data" read by ReadInp2()
    concerning all the infos respecting bc and temperature
    """
    loadings = GetBoundaryConditions(data)

    temperatureLoading = GetParameterFiles(data, parameterName = 'temperature')
    if temperatureLoading:
        loadings['temperature'] = [[["ALLNODE", temperatureLoading]]]


    return loadings




def GetInitDofValues(data):
    """
    returns a dictionary containing the infos of an inp file "data" read by ReadInp2()
    concerning all the infos respecting ``****calcul``, ``***init_dof_value``
    """
    initData = GetFromInp(data,{'4':['calcul'], '3':['init_dof_value']})
    if initData == []:
        return ('uniform', 0.)
    else:
        return (initData[0][0][2], initData[0][0][3])


def GetParameterFiles(data, parameterName = None):
    """
    returns a dictionary containing the infos of an inp file "data" read by ReadInp2()
    concerning all the infos respecting ``****calcul``, ``***parameter``, ``**file``
    !! only with 'cycle_conversion' time table

    extract floats from string: solution from https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string
    """
    if parameterName == None:
      print("Warning, not supported for the moment when parameterName is not specified (dictionary keys may be overwritten)")
      paraFilesData = GetFromInp(data,{'4':['calcul'], '3':['parameter'], '2':['file']})
    else:
      paraFilesData = GetFromInp(data,{'4':['calcul'], '3':['parameter'], '2':['file', parameterName]})

    import re
    numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
    rx = re.compile(numeric_const_pattern, re.VERBOSE)

    parameterFiles = {}
    for i in range(1, len(paraFilesData)):
      if paraFilesData[i][0][0] != 'cycle_conversion':
        parameterFiles[paraFilesData[i][0][0]] = [p for p in paraFilesData[i][0][1:]]
      else:
        parameterFiles['cycle_conversion'] = {}
        parameterFiles['cycle_conversion']['tInit'] = paraFilesData[i][0][1]
        parameterFiles['cycle_conversion']['tEnd']  = paraFilesData[i][0][2]
        parameterFiles['cycle_conversion']['tStep'] = paraFilesData[i][0][3]
        parameterFiles['cycle_conversion']['timeTable'] = []
        parameterFiles['cycle_conversion']['fileTable'] = []
        for j in range(1, len(paraFilesData[i])):
          line = " ".join(paraFilesData[i][j])
          timeStamp = float(rx.findall(line)[0])
          fileName  = line.split("file ")[-1].split(" ")[0]
          parameterFiles['cycle_conversion']['timeTable'].append(timeStamp)
          parameterFiles['cycle_conversion']['fileTable'].append(fileName)
    return parameterFiles


def GetMaterialFiles(data):

    matList = {}

    tempMatList = GetFromInp(data,{'3':['material'], '2':['elset']})

    if tempMatList == []:
        matList['ALLELEMENT'] = GetFromInp(data,{'3':['material'], '1':['file']})[0][0][1]
    else:
        tempData = GetFromInp(data,{'3':['material'], '2':['elset']})
        for data in tempData:
            if data[0][0] == 'elset':
               curElSet = data[0][1]
            elif data[0][0] == 'file':
                matList[curElSet] = data[0][1]

    return matList



def GetProblemType(data):

    calcul = GetFromInp(data,{'4':['calcul']})[0][0]
    if len(calcul) > 1:
        return calcul[1]
    else:
        return "mechanical"


def ReadBinaryFile(fileName):
    return np.fromfile(fileName, dtype=np.float32).byteswap()



def GetDensity(materialFileName):
    res = ReadInp2(materialFileName, startingNstar=3)
    return float(GetFromInp(res,{'3':['behavior'], '2':['coefficient']})[0][1][1])


def GetBehavior(materialFileName):
    res = ReadInp2(materialFileName, startingNstar=3)
    return GetFromInp(res,{'3':['behavior']})[0][0][1]


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
    WriteInp2(res)
    GetFromInp(res,{'4':['calcul'], '3':['parameter'], '2':['file', 'temperature']})
    GetFromInp(res,{'4':['calcul'], '2':['impose_nodal_dof']})
    GetFromInp(res,{'4':['calcul'], '3':['linear_solver']})
    GetProblemType(res)

    GetTables(res)
    GetBoundaryConditions(res)
    GetParameterFiles(res, parameterName = 'temperature')
    GetInputTimeSequence(res)
    print(GetMaterialFiles(res))

    res = ReadInp2(T2.GetTestDataPath() + 'calcul3.inp')
    GetTables(res)
    GetBoundaryConditions(res)
    GetLoadings(res)
    GetInputTimeSequence(res)
    print(GetInitDofValues(res))
    GetMaterialFiles(res)
    print(GetDensity(T2.GetTestDataPath() + 'elas'))
    print(GetBehavior(T2.GetTestDataPath() + 'elas'))
    GetProblemType(res)

    res = ReadInp2(T2.GetTestDataPath() + 'mat', startingNstar=3)
    GetFromInp(res,{'3':['behavior', 'thermal'], '2':['conductivity', 'isotropic']})
    GetFromInp(res,{'3':['behavior', 'thermal'], '2':['coefficient']})



    ReadBinaryFile(T2.GetTestDataPath() + 'UtExample/cube.node')
    return 'ok'


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
