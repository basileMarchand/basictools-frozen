# Read Zebulon sparse matrices and vectors
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import time

def ReadMat(fileName, returnReorderOnly=False):

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
        A = A + A.T
        A.setdiag(A.diagonal()/2)

        return A


def ReadVec(fileName):

    V = []
    with open(fileName, "r") as f:
        for line in f:
            V.append(float(line))

    return np.array(V)


def WriteVec(data, fileName):
    with open(fileName, 'w') as myFile:
        for dat in data:
            myFile.write(str(dat)+'\n')
    myFile.close()
    

def CheckIntegrity():
    import OTTools.IO.ZebulonIO as ZIO
    import OTTools.Helpers.Tests as T
    dataPath = T.GetTestDataPath()
    ZIO.ReadMat(dataPath+'Zmatrix')
    ZIO.ReadMat(dataPath+'Zmatrix', True)
    vec = ZIO.ReadVec(dataPath+'Zvector')
    tempPath = T.TestTempDir.GetTempPath()
    ZIO.WriteVec(vec, tempPath+'Zvector')
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover   
    