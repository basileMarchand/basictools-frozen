import os
import numpy as np
import collections


def TryCommand(string):
    try:
        os.system(string)
    except:
        True 


def ReadVariableFromPythonFile(fileName):
    dic = {}
    with open (fileName, "r") as myfile:
        data = myfile.readlines()
        for d in data:
            try:
                exec(d, dic, dic)
            except SyntaxError:
                continue
    dic.pop('__builtins__', None)
    return dic


def GetIterable(x):
    if isinstance(x, collections.Iterable):
        return x
    else:
        return (x,)

def AssertFloat(value):
    arguments = GetIterable(value)
    for arg in arguments:
        assert (type(arg) == float or type(arg) == np.float64), "value is not of type float or np.float64"
    return


def CheckIntegrity():

    TryCommand("ls")
    AssertFloat(1.)
    AssertFloat((1., 2.))
    AssertFloat([1., 2.])

    import BasicTools.TestData as BasicToolsTestData
    ReadVariableFromPythonFile(BasicToolsTestData.GetTestDataPath()+"test.py")

    return 'ok'

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover

