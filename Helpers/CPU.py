# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import os

def GetNumberOfAvailableCpus():
   SLURM_JOB_CPU_PER_NODE = os.environ.get('SLURM_JOB_CPU_PER_NODE')
   if SLURM_JOB_CPU_PER_NODE is None:
       return 1

class CPU():
    """ Class to help doning the multithreading without using to many cpus
    """
    cpudispo = GetNumberOfAvailableCpus()

    def __init__(self, nbCPUNeeded=-1, WithError=True):
        self.nbCPUNeeded = nbCPUNeeded
        self.nbCPUAllocated = 0
        self.withError = WithError

    def __enter__(self):
        if self.nbCPUNeeded == -1:
            if CPU.cpudispo == 0:
                if self.withError:
                    raise Exception("No more cpus avilable")
                else:
                    self.nbCPUAllocated = 0
            else:
                self.nbCPUAllocated = self.cpudispo
                CPU.cpudispo =0
            return self


        if self.nbCPUNeeded <= CPU.cpudispo:
            self.nbCPUAllocated = self.nbCPUNeeded
            CPU.cpudispo -= self.nbCPUNeeded
        else:
            if CPU.cpudispo == 0:
                if self.withError:
                    raise Exception("No more cpus avilable")
                else:
                    self.nbCPUAllocated = 0
            else:
                self.nbCPUAllocated = self.cpudispo
                CPU.cpudispo =0

        return self

    def __exit__(self, type, value, traceback):
        CPU.cpudispo += self.nbCPUAllocated
        self.nbCPUAllocated = 0

def CheckIntegrity():
    try:
       with CPU(8) as cpu:
        print("cpu ask 8 Allocated  " + str( cpu.nbCPUAllocated))
        print("CPU dispo " + str( CPU.cpudispo))
        with CPU() as cpu2:
            print("cpu2 ask max Allocated  " + str( cpu2.nbCPUAllocated))
            print("CPU2 dispo " + str( CPU.cpudispo))

            with CPU(7,False) as cpu3:
                print("cpu3 ask 7 Allocated  " + str( cpu3.nbCPUAllocated))
                print("CPU3 dispo " + str( CPU.cpudispo))
                CPU()
    except Exception as e:
         print(e)

   #     pass
    print("CPU dispo " + str( CPU.cpudispo))

    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity())
