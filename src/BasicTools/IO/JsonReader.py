#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import json


def ReadJson(fileName=None,text=None):

  if fileName is not None:
    with open(fileName) as data_file:
        return json.load(data_file)

  if text is not None:
        return json.loads(text)

  return None

def CheckIntegrity():

    text = u"""{"menu": {
  "id": "file",
  "value": "File",
  "popup": {
    "menuitem": [
      {"value": "New", "onclick": "CreateNewDoc()"},
      {"value": "Open", "onclick": "OpenDoc()"},
      {"value": "Close", "onclick": "CloseDoc()"}
    ]
  }
}}
"""
    print(ReadJson(text=text))
    from BasicTools.Helpers.Tests import WriteTempFile
    filename = WriteTempFile(filename="test.jso",content=text)

    if filename is not None:
        print(ReadJson(fileName=filename))
    else:
        return "not OK"# pragma: no cover

    if ReadJson() is not None:
        return "not OK"# pragma: no cover

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover

