import json


def ReadJson(fileName=None,text=None):

  if fileName is not None:
    with open(fileName) as data_file:
        return json.load(data_file)

  if text is not None:
        return json.loads(text)

  return None

def CheckIntegrity():

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
