# -*- coding: utf-8 -*-

class MPIInterface():
    """
    Class to test if the mpi inteface is avilable
    This class work even if the module mpi4py is not avilable

    """
    def __init__(self):
        try:
            from  mpi4py import MPI
            self.MPI = MPI
            self.mpiOK = True
            self.parallel = True
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
            if self.size == 1:
                self.parallel = False
        except:
            self.MPI = None
            self.mpiOK = False
            self.parallel = False

            self.comm = None
            self.rank = 0
            self.size = 1

    @classmethod
    def IsParallel(cls):
        """
        return True if the main program was lauched in a mpi enviroment (mpirun)
        """
        return MPIInterface().parallel

    @classmethod
    def Rank(cls):
        """
        return the rank of the current process (0 if not in a mpi enviroment)
        """
        return MPIInterface().rank

    @classmethod
    def Size(cls):
        """
        return the size of the mpi enviroment (0 if not in a mpi enviroment)
        """
        return MPIInterface().size

    def __str__(self):
        """
        return a string representation of the class
        """
        res = "MPIInterface : " + ("ok" if self.mpiOK else "no")+ "\n"
        res += 'Total size is ' + str(self.size)+ "\n"
        res += 'My rank is ' + str(self.rank)
        return res



def CheckIntegrity(GUI=False):
    obj = MPIInterface()
    print(str(obj.__str__() ))
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity(GUI=True)) # pragma: no cover

