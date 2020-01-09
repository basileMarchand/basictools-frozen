# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

__all__ = ["Helpers","IO","FE","TensorTools","Containers","Actions","Linalg"];

__cython_src__ =  ["FE/CythonIntegration.pyx",
                   "FE/EigenSolver.pyx",
                   "FE/NativeIntegration.pyx",
                   "FE/WeakFormNumerical.pyx"]

__cpp_src__ = ["./FE/src_cpp/NativeIntegration.cpp"]

__name__ = "BasicTools"
__copyright_holder__ = "Safran"
__copyright_years__ = "2016-2020"
__copyright__ = f"{__copyright_years__}, {__copyright_holder__}"
__license__ = "BSD 3-Clause License"
__version__ = "1.1"


def main():
    print(f" {__name__} version {__version__}")
    print(f" Copyright (c) {__copyright__}")
    print("")

    import BasicTools.Helpers.Tests
    BasicTools.Helpers.Tests.TestAll( \
            extraToolsBoxs=["BasicTools"],
            coverage={"active": False})
