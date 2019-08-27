# -*- coding: utf-8 -*-
import os
import BasicTools.Helpers.ParserHelper as PH

def GetNumberOfAvailableCpus():
   SLURM_JOB_CPU_PER_NODE = os.environ.get('SLURM_JOB_CPU_PER_NODE')
   if SLURM_JOB_CPU_PER_NODE is None:
       import subprocess
       import platform
       try:
           out = subprocess.check_output(["/usr/bin/squeue","-h","-t","r","-u",str(os.environ.get("USER")),"-w",platform.node(),"-o","'%C'"], shell=False ).decode("utf-8","ignore")
           out = out.strip().strip("'")
           return PH.ReadInt(out)
       except:
           return 1

class CPU():
    """ Class to help doing the multithreading without using to many cpus
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
       with CPU(2) as cpu:
        print("CPU ask 2 Allocated :" + str( cpu.nbCPUAllocated))
        print("CPU dispo " + str( CPU.cpudispo))
        with CPU() as cpu2:
            print("cpu2 ask max available, Allocated : " + str( cpu2.nbCPUAllocated))
            print("CPU2 dispo " + str( CPU.cpudispo))

            with CPU(7,False) as cpu3:
                print("cpu3 ask 7 Allocated  " + str( cpu3.nbCPUAllocated))
                print("CPU3 dispo " + str( CPU.cpudispo))
    except Exception as e:
         print(e)

   #     pass
    print("CPU dispo " + str( CPU.cpudispo))

    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity())
