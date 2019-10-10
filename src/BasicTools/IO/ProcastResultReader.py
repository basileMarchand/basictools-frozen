# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

import numpy as np

import BasicTools.Containers.UnstructuredMesh  as UM

def ReadResult(nbNodes, fileName=None, string=None):


    if fileName is not None:
       string = open(fileName, 'r')
    elif string is not None:
       from io import StringIO
       string = StringIO(string)

    timeSteps   = []
    temperature = np.empty((nbNodes,0),dtype=np.double)


    while (True) :
       line = string.readline()
       if line == "":
          break
    #for line in string:
       l = line.strip('\n').lstrip().rstrip()
       if len(l) == 0: continue

       if l.find("2         4         1         5         2         1")>-1:
           line = string.readline()
           l = line.strip('\n').lstrip().rstrip()
           s = l.split()
           timeSteps.append(int(s[3]))
           newTemp = np.empty((nbNodes,1),dtype=np.double)

           line = string.readline(); line = string.readline(); l = line.strip('\n').lstrip().rstrip()
           while(True):
               s = l.split()
               index = int(s[0])-1
               line = string.readline()
               l = line.strip('\n').lstrip().rstrip()
               s = l.split()
               newTemp[index] = float(s[0])
               line = string.readline()
               l = line.strip('\n').lstrip().rstrip()
               if len(l) == 0: continue
               if l.find("-1") > -1:
                   break
           temperature = np.hstack((temperature, newTemp))
           continue

    return timeSteps, temperature



def CheckIntegrity():

    __teststring = u"""
    -1
    55
TEMPERATURE
ProCAST RESULTS, ESI



         2         4         1         5         2         1
         2         1         0         0
  0.00000E+00
         1
  1.53000E+03
         2
  1.53000E+03
         3
  1.53000E+03
         4
  1.53000E+03
    -1
    -1
    55
TEMPERATURE
ProCAST RESULTS, ESI



         2         4         1         5         2         1
         2         1         0      2140
  2.12800E+03
         1
  5.41073E+02
         2
  5.41557E+02
         3
  5.42009E+02
         4
  5.42137E+02
    -1

"""
    nbNodes = 4
    timeSteps, temperature = ReadResult(nbNodes, string=__teststring)
    print(timeSteps)
    print(temperature)

    from BasicTools.Helpers.Tests import WriteTempFile
    filename = WriteTempFile(filename="procastResultTest",content=__teststring)
    timeSteps, temperature = ReadResult(nbNodes, fileName=filename)



    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
