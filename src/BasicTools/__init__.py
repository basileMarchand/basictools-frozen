#

__all__ = ["Helpers","IO","FE","TensorTools","Containers","Actions","Linalg"];

__cython_src__ =  ["FE/CythonIntegration.pyx",
                   "FE/EigenSolver.pyx",
                   "FE/NativeIntegration.pyx",
                   "FE/WeakFormNumerical.pyx"]

__cpp_src__ = ["./FE/src_cpp/NativeIntegration.cpp"]

__name__ = "BasicTools"
__author__ = "<copyright holder>"
__copyright__ = "Copyright (c) <year>, <copyright holder>"
__license__ = "New BSD License"
__version__ = "1.0"

if __name__ == '__main__':# pragma: no cover
    print(" Cam Tools Python Modules.")
    print(" Safran All Right Reserved. 2016")
    print("")

    import BasicTools.Helpers.Tests
    BasicTools.Helpers.Tests.TestAll(extraToolsBoxs= ["BasicTools"])
