import os


def TryCommand(string):
  try:
    os.system(string)
  except:
    True 
