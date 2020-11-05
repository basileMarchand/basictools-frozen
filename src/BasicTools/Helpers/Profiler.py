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


    def Start(self):
        self.pr.enable()


    def Stop(self):
        self.pr.disable()
        self.ps = pstats.Stats(self.pr, stream=self.s).sort_stats('cumulative')
        self.ps.reverse_order()
        self.ps.print_stats()


    def PlotStats(self, fileName, partRest = 0.2):

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

        tottimesArgSortInv = tottimesArgSort[S.BinarySearch(cumsumtottimes, partRest)+1:][::-1]
        cumtimesArgSortInv = cumtimesArgSort[S.BinarySearch(cumsumcumtimes, partRest)+1:][::-1]


        tottimes = tottimes[tottimesArgSortInv]
        tottimes = np.hstack((tottimes, 1 - np.sum(tottimes)))

        cumtimes = cumtimes[cumtimesArgSortInv]
        cumtimes = np.hstack((cumtimes, 1 - np.sum(cumtimes)))

        functionNamestottime = [functionNames[i] for i in tottimesArgSortInv] + ["rest"]
        functionNamescumtime = [functionNames[i] for i in cumtimesArgSortInv] + ["rest"]


        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2,1)

        axs[0].pie(tottimes, labels=functionNamestottime, autopct='%1.1f%%', normalize = True)
        axs[0].axis('equal')
        axs[0].set_title("tottime")

        axs[1].pie(cumtimes, labels=functionNamescumtime, autopct='%1.1f%%', normalize = True)
        axs[1].axis('equal')
        axs[1].set_title("cumtime")

        plt.savefig(fileName)



    def __str__(self):
        return self.s.getvalue()[:-1] + "   ncalls  tottime  percall  cumtime  percall filename:lineno(function)"


def CheckIntegrity(GUI=False):

    import time

    p = Profiler()
    p.Start()

    time.sleep(0.002)
    print("toto")

    p.Stop()
    print(p)

    p.PlotStats("test")


    return "ok"


if __name__ == '__main__':

    print(CheckIntegrity( GUI=True))

