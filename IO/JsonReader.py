import json


def ReadJson(file=None,text=None):

  if file is not None:
    with open(file) as data_file:
        return json.load(data_file)

  if text is not None:
        return json.loads(text)

  return None

def CheckIntegrity():

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
