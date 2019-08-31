import os
import numpy as np


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



def AssertFloat(value):
    assert (type(value) == float or type(value) == np.float64), "value is not of type float or np.float64"


def CheckIntegrity():

  TryCommand("ls")
  AssertFloat(1.)

  import BasicTools.TestData as BasicToolsTestData
  ReadVariableFromPythonFile(BasicToolsTestData.GetTestDataPath()+"test.py")

  return 'ok'

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover

