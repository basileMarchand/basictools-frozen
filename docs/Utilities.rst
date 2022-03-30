**********
Utilities
**********
BasicTools provides some generic functionalities not related to finite elements problems.
All this functionalities are present in the folder ``BasicTools.Helpers`` of the library.

Configuration
#############

BasicTools has some module level variables accessible to the user.
This variables can be used to alter the behavior of the library.
The following template can be used to configure BasicTools for a specific application.
please refers to the module help for more information.


.. code-block:: python

    # -*- coding: utf-8 -*-

    #import BasicTools.Helpers.BaseOutputObject as BOO
    #print("BOO.useDifferentilatime =", BOO.useDifferentialTime)
    #default True
    # BOO.useDifferentialTime = True
    # BOO.useDifferentialTime = False

    #print("BOO.useFroze_itDecorator =", BOO.useFroze_itDecorator)
    #default False
    #BOO.useFroze_itDecorator = False
    #BOO.useFroze_itDecorator = True

    #import BasicTools.Actions.OpenInParaView as OIP
    #print("OIP.paraviewExec =", OIP.paraviewExec)
    #default "paraview"
    #OIP.paraviewExec = "paraview"

    #import BasicTools.FE.SymWeakForm as SWF
    #print("SWF.testcharacter =", SWF.testcharacter)
    #default "'"
    #SWF.testcharacter = "'"

    #import BasicTools.FE.DofNumbering as DN
    #print("DN.numberingAlgorithm =", DN.numberingAlgorithm)
    #default  "NumpyBase"
    #DN.numberingAlgorithm = "NumpyBase"
    #DN.numberingAlgorithm = "DictBase"
    #DN.numberingAlgorithm = "CppBase"

    #print("DN.__cacheSize__ =", DN.__cacheSize__)
    #default 1
    #DN.__cacheSize__ = 1


    #import BasicTools.FE.Integration as I
    #print("I.UseCpp =", I.UseCpp)
    #default True
    #I.UseCpp = True
    #I.UseCpp = False

    #print("I.UseMultiThread =", I.UseMultiThread)
    #default True
    #I.UseMultiThread =  True
    #I.UseMultiThread =  False

    #print("I.MultiThreadThreshold =", I.MultiThreadThreshold)
    #default 100
    #I.MultiThreadThreshold = 100

    #import BasicTools.Linalg.LinearSolver as IS
    #default "EigenCG"
    #IS.defaultIfError = "EigenCG"
    #IS.defaultIfError = "CG"


