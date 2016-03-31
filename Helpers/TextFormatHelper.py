# -*- coding: utf-8 -*-

# Principal developper:
#   Felipe Bordeu    
    

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
        TFormat.__indentation = 0
        TFormat.__size = 2
        TFormat.__extra = 0;        
        
    @staticmethod
    def GoodBad(text,test):
        if test:
            return TFormat.InGreen(str(text))
        else:
            return TFormat.InRed(str(text))
    
    @staticmethod
    def InRed(text):
        return '\x1b[31m\x1b[1m' + text + '\x1b[0m'
        
    @staticmethod
    def InGreen(text):
        return '\x1b[32m\x1b[1m' + text + '\x1b[0m'

    @staticmethod
    def InYellow(text):
        return '\x1b[33m\x1b[1m' + text + '\x1b[0m'
        
    @staticmethod
    def InBlue(text):
        return '\x1b[34m\x1b[1m' + text + '\x1b[0m'        

    @staticmethod
    def InPurple(text):
        return '\x1b[35m\x1b[1m' + text + '\x1b[0m'        

    @staticmethod
    def InGrey(text):
        return '\x1b[39m\x1b[1m' + text + '\x1b[0m'        
        
    @staticmethod
    def InRedBackGround(text):
        return '\x1b[41m\x1b[1m' + text + '\x1b[0m'        
        
    @staticmethod
    def InGreenBackGround(text):
        return '\x1b[42m\x1b[1m' + text + '\x1b[0m'        
        
    @staticmethod
    def InYellowBackGround(text):
        return '\x1b[43m\x1b[1m' + text + '\x1b[0m'        
        
    @staticmethod
    def InBlueBackGround(text):
        return '\x1b[44m\x1b[1m' + text + '\x1b[0m'   
        
    @staticmethod
    def InPurpleBackGround(text):
        return '\x1b[45m\x1b[1m' + text + '\x1b[0m'   
        
    @staticmethod
    def InGreyBackGround(text):
        return '\x1b[40m\x1b[1m' + text + '\x1b[0m'   

        
    @staticmethod    
    def Center(text,fill= "*",width=60):
        width=width//2-1
        lo = len(text)
        l = (lo/2) if (lo/2)<width else width 
        return ( fill*(width-l)+ ' '  + text+ ' ' + fill*(width-l) +fill*((lo/2) == ((lo+1)/2) ) )
        

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
    r += TFormat.InRed(TFormat.Center("toto"))+'\n'
    r += TFormat.InBlue(TFormat.Center("toto",width=15)) + '\n'
    r += '^______15_____^' + '\n\n'
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
    return 'ok'
    
if __name__ == '__main__':
    CheckIntegrity() # pragma: no cover 