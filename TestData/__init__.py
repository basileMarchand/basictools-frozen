""" 
Help function to access the testDataPath of the library
"""
def GetTestDataPath():
    import os
    return os.path.dirname(os.path.abspath(__file__))+os.sep