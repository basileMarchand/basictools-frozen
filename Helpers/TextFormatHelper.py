# -*- coding: utf-8 -*-

# Principal developper:
#   Felipe Bordeu

import sys
import os

class TFormat(object):
    """ Format Helper class (indentation and color)"""

    __indentation = 0
    __size = 2
    __extra = 0;

    @staticmethod
    def II( extra=0 ):
        """Increase indentation level """
        TFormat.__indentation += 1;
        TFormat.AddToExtra(extra);

    @staticmethod
    def DI(extra=0):
        """Decrease indentation level """
        TFormat.__indentation -= 1;
        TFormat.AddToExtra(-extra);

    @staticmethod
    def AddToExtra(extra=0):
        TFormat.__extra += extra;

    @staticmethod
    def GetIndent():
        """Generate a string with the correct number on spaces """
        return  ' '*(TFormat.__indentation*TFormat.__size+TFormat.__extra)

    @staticmethod
    def Reset():
        """Reset the indentation state"""
        TFormat.__indentation = 0;
        TFormat.__size = 2;
        TFormat.__extra = 0;

    @staticmethod
    def GoodBad(text,test):
        if test:
            return TFormat.InGreen(str(text))
        else:
            return TFormat.InRed(str(text))

    @staticmethod
    def InRed(text=''):
        return TFormat.WithColor(text, '31')
        #return '\x1b[31m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InGreen(text=''):
        return TFormat.WithColor(text, '32')
        #return '\x1b[32m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InYellow(text=''):
        return TFormat.WithColor(text, '33')
        #return '\x1b[33m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InBlue(text=''):
        return TFormat.WithColor(text, '34')
        #return '\x1b[34m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InPurple(text=''):
        return TFormat.WithColor(text, '35')
        #return '\x1b[35m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InGrey(text=''):
        return TFormat.WithColor(text, '39')
        #return '\x1b[39m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InRedBackGround(text=''):
        return TFormat.WithColor(text, '41')
        #return '\x1b[41m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InGreenBackGround(text=''):
        return TFormat.WithColor(text, '42')
        #return '\x1b[42m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InYellowBackGround(text=''):
        return TFormat.WithColor(text, '43')
        #return '\x1b[43m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InBlueBackGround(text=''):
        return TFormat.WithColor(text, '44')
        #return '\x1b[44m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InPurpleBackGround(text=''):
        return TFormat.WithColor(text, '45')
        #return '\x1b[45m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InGreyBackGround(text=''):
        return TFormat.WithColor(text, '40')
        #return '\x1b[40m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def WithColor(text, color):
        if TFormat.SupportsColor():
            # take the \x1b[1m out to make darker colors
            if len(text) :
                return '\x1b['+color+';1m' + text + TFormat.CleanFormat()
            else:
                return '\x1b['+color+';1m'
        else:
            return text# pragma: no cover

    @staticmethod
    def WithUnderline(text=''):
        if TFormat.SupportsColor():
            if len(text) :
                return '\033[4m' + text + TFormat.CleanFormat()
            else:
                return '\033[4m'
        else:
            return text# pragma: no cover

    @staticmethod
    def WithItalic (text=''):
        if TFormat.SupportsColor():
            if len(text) :
                return '\033[3m' + text + TFormat.CleanFormat()
            else:
                return '\033[3m'
        else:
            return text# pragma: no cover

    @staticmethod
    def CleanFormat():
        if TFormat.SupportsColor():
            return '\x1b[0m'
        else:
            return ''# pragma: no cover


    @staticmethod
    def Center(text,fill= "*",width=60):

        if fill != " ":
            if (width - TFormat.__nonAnsiStringLen(text)) > 0 :
                text = ' ' + text

            if (width - TFormat.__nonAnsiStringLen(text)) > 0 :
                text = text + ' '

        #if width < len(text)  :
        #    text = text[0:width-1] + '.'

        lt = TFormat.__nonAnsiStringLen(text)
        fwidth = width - lt
        rfw = max(fwidth//2,0)
        lfw = max(width- lt - rfw ,0)

        return  fill*(lfw)  + text +  fill*(rfw)

    @staticmethod
    def __nonAnsiStringLen(text):
        import pyparsing
        ESC = pyparsing.Literal('\x1b')
        integer = pyparsing.Word(pyparsing.nums)
        escapeSeq = pyparsing.Combine(ESC + '[' + pyparsing.Optional(pyparsing.delimitedList(integer,';')) + pyparsing.oneOf(list(pyparsing.alphas)))

        nonAnsiString = lambda s : pyparsing.Suppress(escapeSeq).transformString(s)

        return len(nonAnsiString(text))
#
    @staticmethod
    def InIPython():
        try:
            __IPYTHON__
            return True# pragma: no cover
        except NameError:
            return False# pragma: no cover

    @staticmethod
    def SupportsColor():

        # from https://github.com/django/django/blob/master/django/core/management/color.py
        """
        Returns True if the running system's terminal supports color,
        and False otherwise.
        """
        plat = sys.platform
        supported_platform = plat != 'Pocket PC' and (plat != 'win32' or 'ANSICON' in os.environ)

        # isatty is not always implemented, #6223.
        is_a_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

        #hack for spyder in IPython
        if TFormat.InIPython():
            return True

        #hack for spyder
        if 'SPYDER_SHELL_ID' in os.environ:
            return True
        if hasattr(sys.stdout, 'color') :# pragma: no cover
            return sys.stdout.color


        if not supported_platform or not is_a_tty:# pragma: no cover
            return False
        return True# pragma: no cover


def CheckIntegrity():
    TFormat.Reset()
    TFormat.II(10);
    TFormat.II(5);
    if(TFormat._TFormat__indentation != 2): raise Exception()
    if(TFormat._TFormat__extra != 15): raise Exception()
    TFormat.DI(10);
    if(TFormat.GetIndent() != "       "): raise Exception()
    if(TFormat._TFormat__indentation != 1): raise Exception()
    if(TFormat._TFormat__extra != 5): raise Exception()
    TFormat.Reset()

    if(TFormat._TFormat__indentation != 0): raise Exception()
    if(TFormat._TFormat__extra != 0): raise Exception()

    print(TFormat.InGreenBackGround(TFormat.InRed("InRed with green background")))
    r  = TFormat.InRedBackGround(TFormat.InGreen("InGreen with red background")) +  '\n'
    r += TFormat.InBlueBackGround(TFormat.InYellow("InYellow with blue backfround")) + '\n'
    r += TFormat.InYellowBackGround(TFormat.InBlue('InBlue with yellow backfround')) + '\n'
    r += TFormat.InGreyBackGround(TFormat.InPurple('InPurple with grey backfround')) + '\n'
    r += TFormat.InPurpleBackGround(TFormat.InGrey('InGrey with purple backfround')) + '\n'
    r += '\n'
    r += TFormat.InRed(TFormat.Center("toto"))+'-\n'
    r += '-'+TFormat.InBlue(TFormat.Center("toto",width=8)) + '-\n'
    r += ' '+TFormat.WithUnderline('^   8  ^') + '\n\n'

    r += '-'+TFormat.InBlue(TFormat.Center("toto",width=7)) + '-\n'
    r += ' '+TFormat.WithUnderline('^  7  ^') + '\n\n'

    r += '-'+TFormat.InBlue(TFormat.Center("toto",width=6)) + '-\n'
    r += ' '+TFormat.WithUnderline('^  6 ^') + '\n\n'

    r += '-'+TFormat.InBlue(TFormat.Center("toto",width=5)) + '-\n'
    r += ' '+TFormat.WithUnderline('^ 5 ^') + '\n\n'

    r += '-'+TFormat.InBlue(TFormat.Center("toto",width=4)) + '-\n'
    r += ' '+TFormat.WithUnderline('^ 4^') + '\n\n'

    r += '-'+TFormat.InBlue(TFormat.Center("toto",width=3)) + '-\n'
    r += ' '+TFormat.WithUnderline('^3^') + '\n\n'

    TFormat.II()
    r += TFormat.GetIndent() + 'One Level Of Indentations \n'
    TFormat.II()
    r += TFormat.GetIndent() + 'Two Level Of Indentations \n'
    TFormat.DI()
    r += TFormat.GetIndent() + 'One Level Of Indentations \n'
    TFormat.DI()
    r += TFormat.GetIndent() + 'Root Level Of Indentations \n'
    print(r)
    print( TFormat.GoodBad("This is a good result : "+str(1),1<10))
    print("This is a bad result : "+ TFormat.GoodBad(10,10<1))
    print(TFormat.WithUnderline("This is a message usign WithUnderline"))
    print(TFormat.WithItalic("This is a message usign WithItalic"))

    print("Warning about the formating :")
    print(TFormat.WithUnderline()+ 'WithUnderline '  +TFormat.InRed()+"Red/WithUnderline"  + TFormat.WithItalic() )
    print("you must clean the format at the end " + TFormat.CleanFormat() + " Good")

    print(TFormat.InIPython())
    return 'ok'

if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover
