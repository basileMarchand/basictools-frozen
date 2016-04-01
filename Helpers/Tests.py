# -*- coding: utf-8 -*-

import traceback
""" 
Help function to access the testDataPath of the library
"""
def GetTestDataPath():
    import os
    return os.path.dirname(os.path.abspath(__file__))+os.sep + '..' + os.sep +'TestData' + os.sep
    
""" 
Function to generate and destroy a temporary directory
"""
class TestTempDir():
    path = None
    
    @classmethod 
    def GetTempPath(cls):
        if cls.path is not None:
            return cls.path
        import tempfile
        import os
        cls.path = tempfile.mkdtemp() + os.sep
        return TestTempDir.path
    
    @classmethod 
    def DeleteTempPath(cls):
        import shutil
        if cls.path is not None:
            shutil.rmtree(cls.path)
        cls.path = None
            
def __RunAndCheck(lis,bp,stopAtFirstError, cover):
    import sys
    import time
    
    if cover :
       import coverage
       ss = [ k for k in lis ]
       cov = coverage.coverage(source=ss, omit=['pyexpat','*__init__.py'])
       cov.start()
       import importlib

       for name  in lis:
           ll = importlib.import_module(name)
           reload(ll)
           
    from OTTools.Helpers.TextFormatHelper import TFormat

    res = {}
    for name in lis:
        bp.Print(TFormat.InBlue(TFormat.Center( "Running Test " +name ,width=60 ) ))
        if lis[name] is None:
            bp.Print(TFormat().InRed(TFormat().GetIndent() + 'Sub Module "' + name + '" does not have the CheckIntegrity function'))
            continue

        try:
            start_time = time.time()
            stop_time = time.time()
            print(lis[name])
            r = lis[name]()
            stop_time = time.time()
            res[name] = r
            if not isinstance(r,str): 
                bp.Print(TFormat.InRed( TFormat().GetIndent() + "Please add a correct return statement in the CheckIntegrity of the module" + name)) 
                #raise Exception()
                r = 'Not OK'
        except :
            bp.Print( "Unexpected error:" + str(sys.exc_info()[0]) )
            res[name] = "error"
            traceback.print_exc(file=bp.stdout_)
            
            r = 'Not OK'
            if stopAtFirstError : 
                bp.Restore()
                raise
            
        if r.lower() == 'ok':
            bp.Print( "OK " +name + " : %.3f seconds " %  (stop_time -start_time ))
        else: 
            bp.Print(TFormat.InRed( "NOT OK !!!! " + name )  )
    
    bp.Restore()
    
    if cover:
        cov.stop()
        cov.save()
        
        # create a temp file 
        tempdir = TestTempDir.GetTempPath()
        # tempdir = 'c:/users/d584808/appdata/local/temp/tmp4ipmul/'
        cov.html_report(directory = tempdir)
        import webbrowser
        import os
        print('Coverage Report in ')
        print(tempdir +"index.html")
        webbrowser.open(tempdir +os.sep+"index.html")
        
    return res

    
def __tryImport(noduleName):
    
    try :   
        import importlib
        tocheck = {}
        m = importlib.import_module(noduleName)
        for submod in [ noduleName+ '.'+x for x in  m.__all__ ]:
            sm = importlib.import_module(submod)
            cif = getattr( sm, "CheckIntegrity", None)            
            if cif is not None:
                tocheck[submod ] = cif
            else :
                try:
                    for subsubmod  in [ submod +'.' + x for x in sm.__all__]:
                        ssm = importlib.import_module(subsubmod)
                        cif  = getattr( ssm, "CheckIntegrity", None)   
                        tocheck[subsubmod ] = cif
                except:
                    tocheck[submod ] = None
                    

        return tocheck
    except:
        print("Error loading module '" + noduleName +"'")
        print("This module will not be tested ")
        return {}
    
def TestAll(modulestotreat=['ALL'], fulloutput=False, stopAtFirstError= False, coverage= False, extraToolsBoxs= None) :

    # calls to print, ie import module1
    from OTTools.Helpers.PrintBypass import PrintBypass
    
    print("Runnig Tests : ")  
    print("--- Begin Test ---")

    tocheck = {}
    
    tocheck.update(__tryImport('OTTools'))
    tocheck.update(__tryImport('rmtools'))
    tocheck.update(__tryImport('TopoTools'))
    if extraToolsBoxs is not None:
        for tool in extraToolsBoxs:
            tocheck.update(__tryImport(tool))        
        
    
    bp = PrintBypass()
    if not fulloutput:
        bp.ToSink()

    if modulestotreat[0] == 'ALL':
        __RunAndCheck(tocheck,bp,stopAtFirstError,coverage)
    else:
        filtered =  dict((k, v) for k, v in tocheck.items() if any(s in k for s in modulestotreat ) )
        __RunAndCheck(filtered,bp,stopAtFirstError,coverage);
        
    
    print("--- End Test ---")
    
def CheckIntegrity():
    return "Ok"
    
if __name__ == '__main__':
    #TestAll() # pragma: no cover 
    #TestAllWithCoverage(  fulloutput=False,stopAtFirstError= True) # pragma: no cover 
    #TestAllWithCoverage( modulestotreat=['IO.XdmfWriter' ], fulloutput=False,stopAtFirstError= True) # pragma: no cover 
    TestAll(modulestotreat=['OTTools'], fulloutput=False,coverage=False,extraToolsBoxs=["MyOtherToolBox"])# pragma: no cover 
    #TestAll(modulestotreat=['ALL'], fulloutput=False,coverage= True,extraToolsBoxs=["MyOtherToolBox"])# pragma: no cover 
    #TestAllWithCoverage(modulestotreat=['T.Formats' ], fulloutput=True)