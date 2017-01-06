# -*- coding: utf-8 -*-

import numpy as np
from OTTools.IO.ReaderBase import ReaderBase
import OTTools.FE.ElementNames as EN
import OTTools.FE.UnstructuredMesh  as UM


def ReadGCode(fileName=None,string=None ):
    reader = GReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read()
    return reader.output

class GReader(ReaderBase):
    def __init__(self):
        super(GReader,self).__init__()


    def Read(self):
        extrudeur = 0

        self.StartReading()
        res = UM.UnstructuredMesh()
        res.nodes = np.empty((0,3), float)

        currentposx = 0.
        currentposy = 0.
        currentposz = 0.
        nodes = []
        thicknes = []
        firsttime = True
        for line in self.filePointer:
            l = line.strip('\n').lstrip().rstrip()
            #empty line
            if len(l) == 0: continue
            #comment
            if l[0] == ";": continue
            l = l.split(';')[0]
            st = l.split()
            if st[0] == "G0" or st[0] == "G1":
                thi = 0
                for s in st:
                    if s[0] == "X":
                        currentposx = float(s[1:])
                    if s[0] == "Y":
                        currentposy = float(s[1:])
                    if s[0] == "Z":
                        currentposz = float(s[1:])
                    if s[0] == "E":
                        val = float(s[1:])
                        if val > extrudeur :
                            thi = 1
                        else:
                            thi = 0
                        extrudeur = val

                if firsttime:
                    firsttime = False
                else :
                    thicknes.append(thi)
                nodes.append(currentposx)
                nodes.append(currentposy)
                nodes.append(currentposz)

                continue
            if st[0] == "G92":
                for s in st:
                    if s[0] == "E":
                        extrudeur = float(s[1:])
                continue

            if st[0] == "M106":
                continue

            print("ignoring line " + str(l) )
        self.EndReading()
        res.nodes = np.reshape(nodes,newshape=(len(nodes)/3, 3))
        res.originalIDNodes = np.arange(res.GetNumberOfNodes())
        elems = res.GetElementsOfType(EN.Bar_2)
        elems.Allocate(res.GetNumberOfNodes()-1)
        elems.connectivity[:,0] = xrange(res.GetNumberOfNodes()-1)
        elems.connectivity[:,1] = xrange(1,res.GetNumberOfNodes())
        res.elemFields = {}
        G = np.array(thicknes)

        G.shape =  (res.GetNumberOfNodes()-1,1)
        res.elemFields['OnOff'] = G
        self.output = res

def CheckIntegrity():
    from OTTools.TestData import GetTestDataPath

    res = ReadGCode(fileName=GetTestDataPath()+ "GCodeTest.gcode" )

    print res

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
