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
        self.ps = pstats.Stats(self.pr, stream=self.s).sort_stats(pstats.SortKey.CUMULATIVE)
        self.ps.reverse_order()
        self.ps.print_stats()
        
        
    def PlotStats(self, fileName, numberOfFunctions = 4):
        lines = self.s.getvalue().split('\n')[5:-3]
        functionNames = [l[46:][-20:] for l in lines]
        parsing = np.array([l.split()[:5] for l in lines])
                
        tottimes = np.array(parsing[:,1], dtype = float)
        cumtimes = np.array(parsing[:,3], dtype = float)

        sumtottimes = np.sum(tottimes)
        sumcumtimes = np.sum(cumtimes)

        tottimesOrder = tottimes.argsort()[::-1][:numberOfFunctions]
        cumtimesOrder = cumtimes.argsort()[::-1][:numberOfFunctions]

        tottimes = tottimes[tottimesOrder]
        sumpartialtottimes = np.sum(tottimes)
        tottimes = np.hstack((tottimes, sumtottimes - sumpartialtottimes))
        functionNamestottime = [functionNames[i] for i in tottimesOrder] + ["rest"]


        cumtimes = cumtimes[cumtimesOrder]
        sumpartialcumtimes = np.sum(cumtimes)
        cumtimes = np.hstack((cumtimes, sumcumtimes - sumpartialcumtimes))
        functionNamescumtime = [functionNames[i] for i in cumtimesOrder] + ["rest"]

        
        tottime = tottimes/np.sum(tottimes)
        cumtime = cumtimes/np.sum(cumtimes)
        
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2,1)
        
        axs[0].pie(tottime, labels=functionNamestottime, autopct='%1.1f%%')
        axs[0].axis('equal')
        axs[0].set_title("tottime")
        
        axs[1].pie(cumtime, labels=functionNamescumtime, autopct='%1.1f%%')
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

