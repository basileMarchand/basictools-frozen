# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

_test = ["Actions",
         "Containers",
         "Bridges",
         "FE",
         "Helpers",
         'ImplicitGeometry',
         "IO",
         "Linalg",
         "TensorTools",
         "Bridges"]

__name__ = "BasicTools"
__copyright_holder__ = "Safran"
__copyright_years__ = "2016-2022"
__copyright__ = "{}, {}".format(__copyright_years__,__copyright_holder__)
__license__ = "BSD 3-Clause License"
__version__ = "1.9.6"


def main():
    print(" {} version {}".format(__name__,__version__))
    print(" Copyright (c) {}".format(__copyright__))
    print("")

    from BasicTools.Helpers.Tests import RunTests
    import sys
    sys.exit(len(RunTests()))
