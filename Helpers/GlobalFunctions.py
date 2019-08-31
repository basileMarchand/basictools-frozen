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
        if type(arg) == np.ndarray:
            arg = arg.flat[0]
        assert (type(arg) == float or type(arg) == np.float64), "value is not of type float or np.float64"
    return


def AssertString(value):
    arguments = GetIterable(value)
    for arg in arguments:
        if type(arg) == np.ndarray:
            arg = arg.flat[0]
        assert (type(arg) == str), "value is not of type string"
    return


def AssertType(myType, value):
    """
    Check of all the elements of "value" are of any of the types in "myType"
    myType and value can be singletons or iterables
    """
    values = GetIterable(value)
    types = GetIterable(myType)
    for val in values:
        answer = []
        for typ in types:
            answer.append(type(val) == typ)
        assert (True in answer), str(value)+" is not of any of the types "+str(types)
    return


def AssertSameLength(value):
    arguments = GetIterable(value)
    length0 = np.array(arguments[0]).flatten().shape[0]
    for i, arg in enumerate(arguments):
        assert (np.array(arguments[i]).flatten().shape[0] == length0), "values do not have same length"
    return


def CheckIntegrity():

    TryCommand("ls")
    AssertFloat(1.)
    AssertFloat((1., 2.))
    AssertFloat([1., 2.])
    AssertFloat(np.array([1., 2.]))
    AssertString(("A", "toto"))
    AssertType([str, float], ("A", "toto"))

    AssertSameLength((np.array([1., 2.]), [3., 4]))

    import BasicTools.TestData as BasicToolsTestData
    ReadVariableFromPythonFile(BasicToolsTestData.GetTestDataPath()+"test.py")

    return 'ok'

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover

