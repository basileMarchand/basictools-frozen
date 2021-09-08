# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import time

class Timer():
    almanac = {}

    def __init__(self, name=None):
        self.name = name
        self.starttime = 0
        self.stoptime = 0
        if not name in Timer.almanac:
            # cumulative time and number of time called
            Timer.almanac[name] = [0.,0]

    def __enter__(self):
        self.Start()

    def __exit__(self, type, value, traceback):
        self.Stop()

    def __str__(self):
        if self.starttime  == 0 and self.stoptime==0:
            return self.PrintTimes()
        res = ""
        val = self.almanac[self.name]
        res += "\n" + str(self.name) + ": " + str(time.time()-self.starttime if self.stoptime == 0 else self.stoptime-self.starttime ) +"  ("+str(val[1])+") : " + '{:6.3e}'.format(val[0]) + " s (mean {:6.3e} s/call)".format(val[0]/val[1] )
        return res

    @classmethod
    def PrintTimes(cls):
        res = ""
        for name, val in cls.almanac.items():
            if name is None: continue
            res += "\n" + str(name) + ":   ("+str(val[1])+") : " + '{:6.3e}'.format(val[0]) + " s (mean {:6.3e} s/call)".format(val[0]/val[1] )
        return res

    def Start(self):
        self.starttime =  time.time()
        return self

    def Stop(self):
        self.stoptime = time.time()
        data = Timer.almanac[self.name]
        data[0] += self.GetDiffTime()
        data[1] += 1

    def Reset(self):
        Timer.almanac = {}

    def GetDiffTime(self):
        return self.stoptime - self.starttime

def CheckIntegrity(GUI=False):

    from BasicTools.Helpers.Timer import Timer
    with Timer("os, sys  import Time"):
        import os
        import sys

    with Timer("Time to Solve"):
        print('toto')

    print(Timer.PrintTimes())
    Timer().Reset()

    with Timer("Time of 1 print"):
        print('toto')

    a = Timer("3 grouped prints").Start()
    print("1 Mississippi")
    print("2 Mississippi")
    print("3 Mississippi")
    a.Stop()

    a = Timer("3 independent prints").Start()
    print("1 Mississippi")
    a.Stop()
    a.Start()
    print("2 Mississippi")
    a.Stop()
    a.Start()
    print("3 Mississippi")
    a.Stop()

    print(Timer.PrintTimes())

    return "ok"

if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))
