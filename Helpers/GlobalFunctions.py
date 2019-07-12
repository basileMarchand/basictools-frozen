import os


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
  return dic
