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
__copyright_years__ = "2016-2023"
__copyright__ = "{}, {}".format(__copyright_years__,__copyright_holder__)
__license__ = "BSD 3-Clause License"
__version__ = "1.9.11"

def Preload(subSubModule= "BasicTools"):
    import time
    st = time.time()
    from BasicTools.Helpers.Tests import __tryImport
    try:
        __tryImport(subSubModule,{},stopAtFirstError=True,modulesToTreat=[],modulesToSkip=[])
        pass
    except Exception as e:
        print(e)
        print('Error Loading module : ' + subSubModule + ' (Current folder'+os.getcwd()+')'  )
        print('-*-*-*-*-*-*> missing CheckIntegrity()??? <*-*-*-*-*-*--'  )
        raise
    #from BasicTools.Helpers.Tests import RunTests
    #RunTests()
    if subSubModule == "BasicTools":
        import BasicTools.Containers.NativeTransfer
        import BasicTools.Linalg.NativeEigenSolver
        import BasicTools.FE.Integrators.NativeIntegration
        import BasicTools.FE.WeakForms.NativeNumericalWeakForm
        import BasicTools.Containers.NativeUnstructuredMesh
        import BasicTools.FE.Numberings.NativeDofNumbering
        import BasicTools.FE.Spaces.NativeSpace
        import BasicTools.Containers.NativeFilters

    return time.time()-st

def main():
    print(" {} version {}".format(__name__,__version__))
    print(" Copyright (c) {}".format(__copyright__))
    print("")


    from BasicTools.Helpers.Tests import RunTests
    import sys
    sys.exit(len(RunTests()))
