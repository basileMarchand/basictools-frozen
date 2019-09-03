# -*- coding: utf-8 -*-

import time

class Timer():
    almanac = {}

    def __init__(self, name=None):
        self.name = name
        self.starttime = 0
        if not name in Timer.almanac:
            # cumulative time and number of time called
            Timer.almanac[name] = [0,0]

    def __enter__(self):
        self.starttime =  time.time()

    def __exit__(self, type, value, traceback):
        stoptime = time.time()
        data = Timer.almanac[self.name]
        data[0] += stoptime-self.starttime
        data[1] += 1

    def __str__(self):
        res = ""
        for name, val in self.almanac.items():
            if name is None: continue
            res += "\n" + str(name) + " ("+str(val[1])+") : " + '{:6.3e}'.format(val[0]) + " s (mean {:6.3e} s/call)".format(val[0]/val[1] )
        return res

    def Start(self):
        self.starttime =  time.time()
        return self

    def Stop(self):
        stoptime = time.time()
        data = Timer.almanac[self.name]
        data[0] += stoptime-self.starttime
        data[1] += 1

    def Reset(self):
        Timer.almanac = {}

def CheckIntegrity(GUI=False):

    from BasicTools.Helpers.Timer import Timer
    with Timer("Numpy import Time"):
        import numpy as np

    with Timer("Time to Solve"):
        print('toto')

    print(Timer())
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

    print(Timer())

    return "ok"

if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))
