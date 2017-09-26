# -*- coding: utf-8 -*-
"""Testing infrastructure for BasicTools extra modules

"""

#for python 2.6+ compatibility
from __future__ import print_function
__author__ = "Felipe Bordeu"

import traceback

import time

"""
python  Tests.py -c -f -s -e <extraModules> -m <moduleFilter>
options :
    -c    To activate coverage and generate a html report
    -f    Full output for all the test
    -s    Stop at first error
    -e    To test extra Modules (-e can be repeated)
    -m    To filter the output by this string (-m can be repeated)
    -d    Dry run do not execute only show what will be executed

"""

from BasicTools.Helpers.which import which

def WriteTempFile(filename,content=None,mode="w" ):
    pfile = TestTempDir.GetTempPath() + filename
    with open(pfile, mode) as f:
        if content is not None:
            f.write( content)
        return pfile
    raise(Exception("Unable ot create file :" + pfile))

class TestTempDir(object):
    """Class to generate and to destroy a temporary directory
    """

    path = None

    @classmethod
    def GetTempPath(cls):
        if cls.path is not None:
            return cls.path
        import tempfile
        import os
        cls.path = tempfile.mkdtemp(prefix="BasicTools_Test_Directory_",suffix="_safe_to_delete") + os.sep
        cls.__saveTempPath()
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
        elif which("nautilus"):
            print(cls.GetTempPath())
            subprocess.Popen(['nautilus',  cls.GetTempPath() ])

    @classmethod
    def SetTempPath(cls,path,create=True):# pragma: no cover
        import os
        cls.path = os.path.abspath(path+os.sep) + os.sep
        if create and not os.path.exists(cls.path):
            os.makedirs(cls.path)
        cls.__saveTempPath()

    #very useful in conbination of a alias
    #alias cdtmp='source ~/.BasicToolsTempPath'
    @classmethod
    def __saveTempPath(cls):
        from os.path import expanduser
        home = expanduser("~")
        import os
        with open(home + os.sep+".BasicToolsTempPath","w") as f:
            f.write("cd " + TestTempDir.path + "\n")
            import stat
            os.chmod(home + os.sep+".BasicToolsTempPath", stat.S_IWUSR | stat.S_IRUSR |stat.S_IXUSR)

def __RunAndCheck(lis,bp,stopAtFirstError,dryrun):# pragma: no cover

    from BasicTools.Helpers.TextFormatHelper import TFormat
    import sys

    res = {}
    for name in lis:
        bp.Print(TFormat.InBlue(TFormat.Center( "Running Test " +name ,width=80 ) ))
        if lis[name] is None:
            bp.Print(TFormat().InRed(TFormat().GetIndent() + 'Sub Module "' + name + '" does not have the CheckIntegrity function'))
            continue

        try:
            start_time = time.time()
            stop_time = time.time()
            #print(lis[name])
            if dryrun:
                r = "Dry Run "
            else:
                r = lis[name]()
            sys.stdout.flush()
            sys.stderr.flush()
            stop_time = time.time()
            res[name] = r
            if not isinstance(r,str):
                bp.Print(TFormat.InRed( TFormat().GetIndent() + "Please add a correct return statement in the CheckIntegrity of the module" + name))
                #raise Exception()
                r = 'Not OK'
        except UserWarning as e :
            sys.stdout.flush()
            sys.stderr.flush()
            bp.Print( "Unexpected Warning:" + str(sys.exc_info()[0]) )
            res[name] = "error"
            traceback.print_exc(file=bp.stdout_)

            r = 'Not OK'
        except:
            sys.stdout.flush()
            sys.stderr.flush()
            bp.Print( "Unexpected Error:" + str(sys.exc_info()[0]) )
            res[name] = "error"
            traceback.print_exc(file=bp.stdout_)

            r = 'Not OK'
            if stopAtFirstError :
                bp.Restore()
                raise


        if r.lower() == 'ok':
            bp.Print( "OK " +name + " : %.3f seconds " %  (stop_time -start_time ))
        else:
            bp.Print(TFormat.InRed( str(r) + " !!!! " + name )  )

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


def TestAll(modulestotreat=['ALL'], fulloutput=False, stopAtFirstError= False, coverage= False, extraToolsBoxs= None,dryrun=False) :# pragma: no cover

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
    print("dryrun: ",end="")
    print(dryrun)

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
    start_time = time.time()
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

        __RunAndCheck(tocheck,bp,stopAtFirstError,dryrun);

        #bp.Restore()
        #now the restore is done automaticaly


    if coverage :
        cov.stop()
        cov.save()

        # create a temp file
        tempdir = TestTempDir.GetTempPath()
        ss = [ ("*"+k.split(".")[-1]+"*") for k in tocheck ]
        cov.html_report(directory = tempdir, include=ss )
        import webbrowser
        import os
        print('Coverage Report in ')
        print(tempdir +"index.html")
        print(cov.report(show_missing=False))
        webbrowser.open(tempdir+"index.html")

    stop_time = time.time()
    bp.Print( "Total Time : %.3f seconds " %  (stop_time -start_time ))

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
          opts, args = getopt.getopt(sys.argv[1:],"hcfsde:m:")
      except getopt.GetoptError:
          print( 'python  Tests.py -c -f -e <extraModules> -m <moduleFilter>')
          print( 'options :')
          print( '       -c    To activate coverage and generate a html report')
          print( '       -f    Full output for all the test')
          print( '       -s    Stop at first error')
          print( '       -e    To test extra Modules (-e can be repeated)')
          print( '       -m    To filter the output by this string (-m can be repeated)')
          print( '       -d    Dry run do not execute only show what will be executed')
          sys.exit(2)

      coverage = False
      fulloutput = False
      stopAtFirstError = False
      extraToolsBoxs = []
      modulestotreat=[]
      dryrun = False



      for opt, arg in opts:
         if opt == '-h':
             print( 'python  Tests.py -c -f -s -e <extraModules> -m <moduleFilter>')
             print( 'options :')
             print( '       -c    To activate coverage and generate a html report')
             print( '       -f    Full output for all the test')
             print( '       -s    Stop at first error')
             print( '       -e    To test extra Modules (-e can be repeated)')
             print( '       -m    To filter the output by this string (-m can be repeated)')
             print( '       -d    Dry run do not execute only show what will be executed')
             sys.exit()
         elif opt in ("-c"):
            coverage = True
         elif opt in ("-f"):
            fulloutput = True
         elif opt in ("-s"):
            stopAtFirstError = True
         elif opt in ("-d"):
            dryrun = True
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
                extraToolsBoxs=extraToolsBoxs,
                dryrun= dryrun )



