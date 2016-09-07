# -*- coding: utf-8 -*-
from OTTools.Helpers.BaseOutputObject import BaseOutputObject
from OTTools.Helpers.TextFormatHelper import TFormat as TFormat

class WriterBase(BaseOutputObject):
    def __init__(self, fileName = None):
        self.fileName = None;
        self.SetFileName(fileName)
        self.__isOpen = False

    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Open(self, filename = None):
        if self.__isOpen :
            print(TFormat.InRed("The file is already open !!!!!"))
            raise Exception


        if filename is not None:
            self.SetFileName(filename)

        ## we use unbuffered so we can repaire broken files easily
        try :
            self.filePointer = open(self.fileName, 'w',0)
        except:
            print(TFormat.InRed("Error File Not Open"))
            raise

        self.__isOpen = True

    def Close(self):
        if self.__isOpen:
            self.filePointer.close()
            self.__isOpen = False
        else :
            self.PrintVerbose(TFormat.InRed("File Not Open"))

