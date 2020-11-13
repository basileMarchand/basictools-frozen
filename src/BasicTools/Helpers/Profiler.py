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

    def __init__(self, discardedPortion = 0.2, timerSignificantDigits = 3):
        self.pr = cProfile.Profile()
        self.s = io.StringIO()
        self.ps = None
        self.discardedPortion = discardedPortion
        self.timerSignificantDigits = timerSignificantDigits
        self.totTimes = None
        self.cumTimes = None
        self.sumTotTime = None
        self.sumCumTime = None


    def Start(self):
        self.pr.enable()


    def Stop(self):

        self.pr.disable()
        self.ps = pstats.Stats(self.pr, stream=self.s).sort_stats('cumulative')
        self.ps.reverse_order()
        self.ps.print_stats()


        lines = self.s.getvalue().split('\n')[5:-3]
        functionNames = [l[46:][-60:] for l in lines]
        parsing = np.array([l.split()[:5] for l in lines])

        tottimes = np.array(parsing[:,1], dtype = float)
        cumtimes = np.array(parsing[:,3], dtype = float)
        self.sumTotTime = np.sum(tottimes)
        self.sumCumTime = np.sum(cumtimes)

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

        functionNamestottime = [functionNames[i] for i in tottimesArgSortInv] + ["REMAINING"]
        functionNamescumtime = [functionNames[i] for i in cumtimesArgSortInv] + ["REMAINING"]


        from collections import OrderedDict

        self.totTimes = OrderedDict(zip(tottimes, functionNamestottime))
        self.cumTimes = OrderedDict(zip(cumtimes, functionNamescumtime))



    def __str__(self):

        from BasicTools.Helpers.TextFormatHelper import TFormat
        string = TFormat.InBlue("Profiler")+" discarding "
        string += str(int(100.*self.discardedPortion))+"% of the smallest functions\n"
        string += TFormat.InRed("Total times:\n")

        for part, name in self.totTimes.items():
            string += TFormat.InGreen(TFormat.Center(name,width=70))+": "+str(round(100.*part,2))
            string += "% ("+str(round(self.sumTotTime*part,self.timerSignificantDigits))+"s)\n"
        string += TFormat.InRed("Cumulated  times:\n")
        for part, name in self.cumTimes.items():
            string += TFormat.InGreen(TFormat.Center(name,width=70))+": "+str(round(100.*part,2))
            string += "% ("+str(round(self.sumCumTime*part,self.timerSignificantDigits))+"s)\n"
        return string





def CheckIntegrity(GUI=False):

    import time

    p = Profiler(0.2)

    p.Start()
    time.sleep(0.002)
    p.Stop()

    print(p)


    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))

