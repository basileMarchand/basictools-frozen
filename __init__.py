#

__all__ = ["Helpers","IO","FE","T","Containers","Actions","Linalg"];

__cython_src__ =  ["FE/CythonIntegration.pyx",
                   "FE/EigenSolver.pyx",
                   "FE/NativeIntegration.pyx",
                   "FE/WeakFormNumerical.pyx"]

__cpp_src__ = ["./FE/src_cpp/NativeIntegration.cpp"]

if __name__ == '__main__':# pragma: no cover
    print(" Cam Tools Python Modules.")
    print(" Safran All Right Reserved. 2016")
    print("")

    import BasicTools.Helpers.Tests
    BasicTools.Helpers.Tests.TestAll(extraToolsBoxs= ["BasicTools"])
