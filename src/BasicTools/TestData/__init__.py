#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
""" 
Help function to access the testDataPath of the library
"""
def GetTestDataPath():
    import os
    return os.path.dirname(os.path.abspath(__file__))+os.sep