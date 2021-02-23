# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

_test = ["Actions",
           "Containers",
           "FE",
           "Helpers",
           'ImplicitGeometry',
           "IO",
           "Linalg",
           "TensorTools"];


__name__ = "BasicTools"
__copyright_holder__ = "Safran"
__copyright_years__ = "2016-2021"
__copyright__ = "{}, {}".format(__copyright_years__,__copyright_holder__)
__license__ = "BSD 3-Clause License"
__version__ = "1.5"


def main():
    print(" {} version {}".format(__name__,__version__))
    print(" Copyright (c) {}".format(__copyright__))
    print("")

    import BasicTools.Helpers.Tests
    BasicTools.Helpers.Tests.TestAll( \
            extraToolsBoxs=["BasicTools"],
            coverage={"active": False})
