import json


def ReadJson(file=None,text=None):

  if file is not None:
    with open(file) as data_file:
        return json.load(data_file)

  if text is not None:
        return json.loads(text)

  return None
