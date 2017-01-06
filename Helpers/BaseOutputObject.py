# -*- coding: utf-8 -*-

from OTTools.Helpers.TextFormatHelper import TFormat

""" Base object to Help input and output """


class BaseOutputObject(object):
    __globalDebugMode = False
    __verboseLevel = 1

    def __init__(self, other = None):
        super(BaseOutputObject,self).__init__()
        if other is not None:
            self.__classDebugMode = other.__classDebugMode
        else:
            self.__classDebugMode = False

    @classmethod
    def SetVerboseLevel(cls,level):
        BaseOutputObject.__verboseLevel = level

    @classmethod
    def SetGlobalDebugMode(cls, mode = True):
        BaseOutputObject.__globalDebugMode = mode

    def SetInstanceDebugMode(self, mode = True):
        self.__classDebugMode = mode

    def PrintProgress(self,val,  maxx = 100,minn= 0):
        per = (float(val) -minn)*100/(maxx-minn)
        if per == round(per):
            self.Print(str(round(per)) +"%")

    def PrintDebug(self, mess):
        self.PrintInternal( mess, 3)

    def PrintVerbose(self, mess):
        self.PrintInternal(mess, 2)

    def Print(self, mess):
        self.PrintInternal(mess, 1)

    def PrintError(self, mess):
        self.PrintInternal(TFormat.InRed("Error: ") + str(mess), 1)

    def InDebugMode(self):
        return (BaseOutputObject.__globalDebugMode or self.__classDebugMode)

    @classmethod
    def Print2(cls,mess):
        a = BaseOutputObject()
        a.PrintInternal(mess,1)

    def PrintInternal(self, mess, level=1):
        if BaseOutputObject.__globalDebugMode or self.__classDebugMode :
            res = ""
            ## we only load modules in debug mode
            import inspect
            import time

            stnumber = 1
            stack = inspect.stack()[stnumber]

            # check if the we call PrintInternal directly of using the PrintDebug proxys
            if 'PrintInternal' in stack[4][0] :
                stnumber += 1
            # the real stack
            stack = inspect.stack()[stnumber]

            # nice date
            m, s = divmod(time.time(), 60.)
            h, m = divmod(m, 60)
            d, h = divmod(h, 24)
            res += ("%d:%02d:%02d" % (h+1, m, s))

            res += (": "+str(stack[1]) + ":" +str(stack[2]) )
            if level == 1 :
                res += (TFormat.InBlue(" --> "))
            elif level == 2:
                res += (TFormat.InGreen(" V-> "))
            else:
                res += (TFormat.InRed(" D-> "))
            # we recover the code of the line to print it
            # some cleaning
            if type(mess) is not str:
                stack = inspect.stack()[stnumber]
                res += (")".join("(".join(str(stack[4][0]).split("(")[1:]).split(")")[0:-1]) )
                res += (" -> ")
            #print(" [" + str(memory()) + "]"),
            res += (str(mess))
            res += "\n"
            print(res),
            return
        elif level <= BaseOutputObject.__verboseLevel :
            print(mess)

        import sys
        sys.stdout.flush()


def CheckIntegrity():

    myObj = BaseOutputObject()

    myObj.Print("Print 1")
    myObj.PrintVerbose("PrintVerbose 1")
    myObj.PrintDebug("PrintDebug 1")
    myObj.SetVerboseLevel(2)
    myObj.Print("Print 2")
    myObj.PrintVerbose("PrintVerbose 2")
    myObj.PrintDebug("PrintDebug 2")
    myObj.SetInstanceDebugMode()
    myObj.Print("Print 3")
    myObj.PrintVerbose("PrintVerbose 3")
    myObj.PrintDebug("PrintDebug 3")
    myObj.SetInstanceDebugMode(False)
    myObj.Print("Print 4")
    myObj.PrintVerbose("PrintVerbose 4")
    myObj.PrintDebug("PrintDebug 4")
    myObj.SetGlobalDebugMode(True)
    myObj.Print("Print 5")
    myObj.PrintVerbose("PrintVerbose 5")
    myObj.PrintDebug("PrintDebug 5")

    myObj.Print([3+1,8])
    myObj.PrintVerbose([3+2,8])
    myObj.PrintDebug([3+3,8])

    myObj.Print(range(4))

    var = 4
    myObj.Print(var)

    myObj.Print2("Print 5")
    myObj.PrintError("Print 6")
    myObj.InDebugMode()

    myObj.PrintProgress(50)

    myObj2 = BaseOutputObject(myObj)

    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
