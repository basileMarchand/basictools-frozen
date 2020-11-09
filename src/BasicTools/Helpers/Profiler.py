# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import cProfile
import pstats
import io
import numpy as np


class Profiler():

    def __init__(self):
        self.pr = cProfile.Profile()
        self.s = io.StringIO()
        self.ps = None
        self.discardedPortion = 0.2
        self.totTimes = None
        self.cumTimes = None


    def Start(self):
        self.pr.enable()


    def Stop(self):
        self.pr.disable()
        self.ps = pstats.Stats(self.pr, stream=self.s).sort_stats('cumulative')
        self.ps.reverse_order()
        self.ps.print_stats()


    def SetDiscardedPortion(self, discardedPortion):
        self.discardedPortion = discardedPortion


    def SortStats(self):

        lines = self.s.getvalue().split('\n')[5:-3]
        functionNames = [l[46:][-30:] for l in lines]
        parsing = np.array([l.split()[:5] for l in lines])

        tottimes = np.array(parsing[:,1], dtype = float)
        cumtimes = np.array(parsing[:,3], dtype = float)

        tottimes = tottimes / np.sum(tottimes)
        cumtimes = cumtimes / np.sum(cumtimes)

        tottimesArgSort = np.argsort(tottimes)
        cumtimesArgSort = np.argsort(cumtimes)

        cumsumtottimes = np.cumsum(tottimes[tottimesArgSort])
        cumsumcumtimes = np.cumsum(cumtimes[cumtimesArgSort])

        from BasicTools.Helpers import Search as S

        tottimesArgSortInv = tottimesArgSort[S.BinarySearch(cumsumtottimes, self.discardedPortion)+1:][::-1]
        cumtimesArgSortInv = cumtimesArgSort[S.BinarySearch(cumsumcumtimes, self.discardedPortion)+1:][::-1]


        tottimes = tottimes[tottimesArgSortInv]
        tottimes = np.hstack((tottimes, 1 - np.sum(tottimes)))

        cumtimes = cumtimes[cumtimesArgSortInv]
        cumtimes = np.hstack((cumtimes, 1 - np.sum(cumtimes)))

        functionNamestottime = [functionNames[i] for i in tottimesArgSortInv] + ["rest"]
        functionNamescumtime = [functionNames[i] for i in cumtimesArgSortInv] + ["rest"]


        from collections import OrderedDict

        self.totTimes = OrderedDict(zip(tottimes, functionNamestottime))
        self.cumTimes = OrderedDict(zip(cumtimes, functionNamescumtime))



    def __str__(self):

        from BasicTools.Helpers.TextFormatHelper import TFormat
        string = TFormat.InBlue("Profiler")+\
        " discarding "+str(int(100.*self.discardedPortion))+"% of the smallest functions\n"
        string += TFormat.InRed("total times:\n")+str(self.totTimes)+"\n"
        string += TFormat.InRed("cumulated times:\n")+str(self.cumTimes)
        return string


        #return self.s.getvalue()[:-1] + "   ncalls  tottime  percall  cumtime  percall filename:lineno(function)"


def CheckIntegrity(GUI=False):

    import time

    p = Profiler()
    p.Start()
    p.SetDiscardedPortion(0.2)

    time.sleep(0.002)
    print("toto")

    p.Stop()
    p.SortStats()
    print(p)



    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))

