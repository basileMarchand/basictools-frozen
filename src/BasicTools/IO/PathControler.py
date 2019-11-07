# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

from os import path
import os
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Helpers.Tests import TestTempDir

class TemporalChdir(BaseOutputObject):
    def __init__(self,targetPath):
        super(TemporalChdir,self).__init__()
        self.originalPath = os.getcwd()+os.sep
        self.targetPath = targetPath

    def __enter__(self):
        os.chdir(self.targetPath)

    def __exit__(self, type, value, traceback):
        os.chdir(self.originalPath)


class PathControler(BaseOutputObject):
    currentDirectory = os.getcwd()+os.sep
    workingDirectory = os.getcwd()+os.sep

    def __init__(self):
        super(PathControler,self).__init__()

    @staticmethod
    def SetCurrentDirectory(folder):
        PathControler.currentDirectory = path.abspath(path.expanduser(folder))+os.sep

    @staticmethod
    def SetCurrentDirectoryUsingFile(folder):
        PathControler.currentDirectory = path.dirname(path.abspath(path.expanduser(folder)))+os.sep

    @staticmethod
    def SetWorkingDirectory(folder):
        PathControler.workingDirectory = path.abspath(path.expanduser(folder))+os.sep

    @staticmethod
    def SetWorkingDirectoryUsingFile(file):
        PathControler.workingDirectory = path.abspath(path.dirname(path.expanduser(file)) )+os.sep

    @staticmethod
    def GetCurrentDirectory():
        return PathControler.currentDirectory

    @staticmethod
    def GetWorkingDirectory():
        return PathControler.workingDirectory


    @staticmethod
    def GetFullFilenameCurrentDirectory(filename,onpath=None):
        if os.path.isabs(filename):
            return filename
        else:
            if onpath is None:
                return path.abspath(PathControler.currentDirectory +filename)
            else:
                return path.abspath(onpath +filename)

    @staticmethod
    def GetFullPathInCurrentDirectoryUsingFile(filename):
        file = PathControler.GetFullFilenameCurrentDirectory(filename)
        return path.abspath(path.dirname(path.expanduser(file)) )+os.sep

    @staticmethod
    def GetFullPathnameOnWorkingDirectory(pathname):
        if os.path.isabs(pathname):
            return pathname
        else:
            return path.abspath(PathControler.workingDirectory +pathname)+os.sep

    @staticmethod
    def GetFullFilenameOnWorkingDirectory(filename):
        if os.path.isabs(filename):
            return filename
        else:
            return path.abspath(PathControler.workingDirectory +filename)

    @staticmethod
    def GetFullFilenameOnTempDirectory(filename):
        if os.path.isabs(filename):
            return filename
        else:
            return path.abspath(TestTempDir().GetTempPath() + filename)

def CheckIntegrity():
    print("C: " +PathControler.currentDirectory)
    print("W: " +PathControler.workingDirectory)

    print("C: " + PathControler.GetFullFilenameCurrentDirectory("tata"))
    print("W: " +PathControler.GetFullFilenameOnWorkingDirectory("tete"))
    print("T: " +PathControler.GetFullFilenameOnTempDirectory("titi"))
    print("*************")

    PathControler.SetCurrentDirectory("~")
    PathControler.SetWorkingDirectory("/tmp")

    print("C: " +PathControler.GetFullFilenameCurrentDirectory("tata"))
    print("W: " +PathControler.GetFullFilenameOnWorkingDirectory("tete"))
    print("T: " +PathControler.GetFullFilenameOnTempDirectory("titi"))
    print("*************")
    PathControler.SetCurrentDirectory("..")
    PathControler.SetWorkingDirectory("/tmp/")
    print("C: " +PathControler.GetFullFilenameCurrentDirectory("tata"))
    print("W: " +PathControler.GetFullFilenameOnWorkingDirectory("../tete"))
    print("T: " +PathControler.GetFullFilenameOnTempDirectory("../../titi"))


    print(PathControler.GetCurrentDirectory())
    print(PathControler.GetWorkingDirectory())


    return("ok")
if __name__ == '__main__':
   print(CheckIntegrity())# pragma: no cover
