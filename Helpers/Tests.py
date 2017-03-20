# -*- coding: utf-8 -*-
from __future__ import print_function

import traceback

"""
Function to generate and destroy a temporary directory
"""

class TestTempDir(object):
    path = None

    @classmethod
    def GetTempPath(cls):
        if cls.path is not None:
            return cls.path
        import tempfile
        import os
        cls.path = tempfile.mkdtemp(prefix="BasicTools_Test_Directory_",suffix="_safe_to_delete") + os.sep
        return TestTempDir.path

    #  we cant test this funciotn, because the temp path will be delete
    @classmethod
    def DeleteTempPath(cls):# pragma: no cover
        import shutil
        if cls.path is not None:
            shutil.rmtree(cls.path)
        cls.path = None

    @classmethod
    def OpenTempFolder(cls):# pragma: no cover
        import subprocess
        import os
        if os.name == "nt":
            subprocess.Popen('explorer "' + cls.GetTempPath() +'"')

    @classmethod
    def SetTempPath(cls,path):# pragma: no cover
        import os
        cls.path = os.path.abspath(path+os.sep) + os.sep

def __RunAndCheck(lis,bp,stopAtFirstError):# pragma: no cover

    from BasicTools.Helpers.TextFormatHelper import TFormat
    import sys
    import time

    res = {}
    for name in lis:
        bp.Print(TFormat.InBlue(TFormat.Center( "Running Test " +name ,width=60 ) ))
        if lis[name] is None:
            bp.Print(TFormat().InRed(TFormat().GetIndent() + 'Sub Module "' + name + '" does not have the CheckIntegrity function'))
            continue

        try:
            start_time = time.time()
            stop_time = time.time()
            #print(lis[name])
            r = lis[name]()
            sys.stdout.flush()
            sys.stderr.flush()
            stop_time = time.time()
            res[name] = r
            if not isinstance(r,str):
                bp.Print(TFormat.InRed( TFormat().GetIndent() + "Please add a correct return statement in the CheckIntegrity of the module" + name))
                #raise Exception()
                r = 'Not OK'
        except :
            sys.stdout.flush()
            sys.stderr.flush()
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

    return res

def __tryImportRecursive(submod,tocheck,stopAtFirstError):


  import importlib
  sm = importlib.import_module(submod)

  cif = getattr( sm, "CheckIntegrity", None)
  if cif is not None:
     tocheck[submod ] = cif
  else :
     try:
         for subsubmod  in [ submod +'.' + x for x in sm.__all__]:
             try:
                 __tryImportRecursive(subsubmod,tocheck,stopAtFirstError)
             except:
                 print('Error Loading File : ' + subsubmod + '.py'  )
                 if(stopAtFirstError): raise
     except:
        tocheck[submod ] = None
        if(stopAtFirstError): raise


def __tryImport(noduleName,stopAtFirstError):# pragma: no cover

    tocheck = {}
    try :
        __tryImportRecursive(noduleName,tocheck,stopAtFirstError)
    except:
        print("Error loading module '" + noduleName +"'")
        print("This module will not be tested ")
        if(stopAtFirstError): raise
    return tocheck


def TestAll(modulestotreat=['ALL'], fulloutput=False, stopAtFirstError= False, coverage= False, extraToolsBoxs= None) :# pragma: no cover

    print("")
    print("modulestotreat   : ",end="")
    print(modulestotreat)
    print("fulloutput       : ",end="")
    print(fulloutput)
    print("stopAtFirstError : ",end="")
    print(stopAtFirstError)
    print("coverage         : ",end="")
    print(coverage)
    print("extraToolsBoxs: ",end="")
    print(extraToolsBoxs)

    cov = None
    if coverage :
       import coverage
       #ss = [ k for k in lis ]
       #cov = coverage.coverage(source=ss, omit=['pyexpat','*__init__.py'])
       cov = coverage.coverage(omit=['pyexpat','*__init__.py'])
       cov.start()

    # calls to print, ie import module1
    from BasicTools.Helpers.PrintBypass import PrintBypass

    print("Runnig Tests : ")
    print("--- Begin Test ---")

    tocheck = {}


    if extraToolsBoxs is not None:
        for tool in extraToolsBoxs:
            tocheck.update(__tryImport(tool,stopAtFirstError))


    with PrintBypass() as bp:
        if not fulloutput:
            bp.ToSink()
            bp.Print("Sending all the output of the tests to sink")


        if modulestotreat[0] is not  'ALL':
            #bp.Print(str(modulestotreat))
            filtered =  dict((k, v) for k, v in tocheck.items() if all(s in k for s in modulestotreat ) )
            #bp.Print(str(filtered))
            tocheck = filtered

        __RunAndCheck(tocheck,bp,stopAtFirstError);

        #bp.Restore()
        #now the restore is done automaticaly


    if coverage :
        cov.stop()
        cov.save()

        # create a temp file
        tempdir = TestTempDir.GetTempPath()
        # tempdir = 'c:/users/d584808/appdata/local/temp/tmp4ipmul/'
        ss = [ ("*"+k.split(".")[-1]+"*") for k in tocheck ]
        cov.html_report(directory = tempdir, include=ss )
        import webbrowser
        import os
        print('Coverage Report in ')
        print(tempdir +"index.html")
        print(cov.report(show_missing=False))
        webbrowser.open(tempdir+"index.html")


    print("--- End Test ---")

def CheckIntegrity():
    TestTempDir().GetTempPath()
    TestTempDir().GetTempPath()
    return "Ok"

if __name__ == '__main__':# pragma: no cover
    import sys, getopt
    if len(sys.argv) == 1:
        TestAll(modulestotreat=['ALL'],extraToolsBoxs= ["BasicTools"], fulloutput=False,coverage=False)# pragma: no cover
    else:
      try:
          opts, args = getopt.getopt(sys.argv[1:],"hcfse:m:")
      except getopt.GetoptError:
          print( 'python  Tests.py -c -f -e <extraModules> -m <moduleFilter>')
          print( 'options :')
          print( '       -c    To activate coverage and generate a html report')
          print( '       -f    Full output for all the test')
          print( '       -s    Stop at first error')
          print( '       -e    To test extra Modules (-e can be repeated)')
          print( '       -m    To filter the output by this string (-m can be repeated)')
          sys.exit(2)

      coverage = False
      fulloutput = False
      stopAtFirstError = False
      extraToolsBoxs = []
      modulestotreat=[]


      for opt, arg in opts:
         if opt == '-h':
             print( 'python  Tests.py -c -f -s -e <extraModules> -m <moduleFilter>')
             print( 'options :')
             print( '       -c    To activate coverage and generate a html report')
             print( '       -f    Full output for all the test')
             print( '       -s    Stop at first error')
             print( '       -e    To test extra Modules (-e can be repeated)')
             print( '       -m    To filter the output by this string (-m can be repeated)')
             sys.exit()
         elif opt in ("-c"):
            coverage = True
         elif opt in ("-f"):
            fulloutput = True
         elif opt in ("-s"):
            stopAtFirstError = True
         elif opt in ("-e"):
            extraToolsBoxs.append(arg)
         elif opt in ("-m"):
            modulestotreat.append(arg)

      if len(modulestotreat) == 0:
         modulestotreat.append("ALL")

      if len(extraToolsBoxs) == 0:
         extraToolsBoxs.append("BasicTools")


      TestAll(  modulestotreat=modulestotreat,
                coverage=coverage,
                fulloutput=fulloutput,
                stopAtFirstError= stopAtFirstError,
                extraToolsBoxs=extraToolsBoxs )



