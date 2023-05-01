# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Post file reader
"""
import re
import numpy as np

from BasicTools.IO.ReaderBase import ReaderBase

class PostReader(ReaderBase):
    def __init__(self, fileName=None):
        super(PostReader,self).__init__(fileName=fileName)
        self.refsAsAField = True

    def Read(self):
        """Function that performs the reading of a post file

        Returns
        -------
        dict
            read data
        """
        self.StartReading()
        res = {}

        while(True):
            line = self.filePointer.readline()
            if line == "" :
                break

            l = line.strip('\n')
            if len(l) == 0:
                continue

            if "#volume integrate" in line :
                continue

            if l[0] == "#":
                if "<" in l:
                    s = l.split("<")
                else:
                    s = re.sub('#',' ',re.sub(' +',' ',l)).split()

                if len(s) == 1:
                    l = "time"+ l
                    s = l.split("#")
                names = [re.sub("[>#]" ,"" ,x).lstrip().rstrip() for x in s ]
                self.PrintVerbose("Reading : "+ " ".join(names))
                data = np.empty((0,len(s)),dtype=np.float32)

                while(True):
                    line = self.filePointer.readline()
                    if line == "" :
                        break
                    l = line.strip('\n')
                    if len(l) == 0:
                        break

                    l= l.replace("VI","")
                    data = np.vstack((data,np.fromstring(l,sep=" ") ))

                #print(data)
                for i in range(len(names)):
                    res[names[i]] = data[:,i]
            else:
                self.Print("dont know how t treat " + str(l))# pragma: no cover
        return res

def CheckIntegrity():
    string = u"""
# Volume
1.000000000000000e+00 2.536237549886229e+05

# time <ener>
1.000000000000000e+00 -nan

# time <enerintegII>
1.000000000000000e+00 1.054915810911460e-01

#volume integrate
# time enerinteg
1.000000000000000e+00 VI 4.644104676325955e-02

"""
    theReader = PostReader()
    theReader.SetStringToRead(string)
    res = theReader.Read()
    print(res)
    return 'OK'


if __name__ == '__main__':
    from BasicTools.Helpers.BaseOutputObject import BaseOutputObject# pragma: no cover
    BaseOutputObject.SetGlobalDebugMode(True)# pragma: no cover
    print(CheckIntegrity())# pragma: no cover
