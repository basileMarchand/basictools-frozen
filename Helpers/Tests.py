# -*- coding: utf-8 -*-
"""Testing infrastructure for BasicTools extra modules

"""

#for python 2.6+ compatibility
from __future__ import print_function
__author__ = "Felipe Bordeu"

import traceback
import time


Test_Help_String = """
python  Tests.py -c -f -s -e <extraModules> -m <moduleFilter>
options :
    -c    To activate coverage and generate a html report
          -cb (coverage with browser and local file index.html generated)
    -f    Full output for all the test
    -s    Stop at first error
    -e    To test extra Modules (-e can be repeated)
    -m    To filter the output by this string (-m can be repeated)
    -d    Dry run do not execute anything, only show what will be executed
    -p    Activate profiling
"""

from BasicTools.Helpers.which import which


def GetUniqueTempFile(suffix="",prefix='tmp'):
    #solution from
    # https://stackoverflow.com/questions/2961509/python-how-to-create-a-unique-file-name

    import tempfile
    (fd, filename) = tempfile.mkstemp(suffix=suffix, prefix=prefix,dir=TestTempDir.GetTempPath())
    return (fd,filename)


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

def __RunAndCheck(lis,bp,stopAtFirstError,dryrun,profiling):# pragma: no cover

    res = {}
    from BasicTools.Helpers.TextFormatHelper import TFormat
    import sys


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

                if profiling :
                    import cProfile, pstats
                    #for python2
                    try:
                        from StringIO import StringIO
                    except:
                        #python3
                        from io import StringIO
                    pr = cProfile.Profile()
                    pr.enable()
                    r = lis[name]()
                    pr.disable()
                    s = StringIO()
                    sortby = 'cumulative'
                    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
                    ps.print_stats()
                    print(s.getvalue())
#                    glob = {"lis":lis,"name":name}
#                    loc = {}
#                    cProfile.runctx("r = lis[name]()",glob,loc)
#                    r = loc["r"]
                else:
                    r = lis[name]()
            sys.stdout.flush()
            sys.stderr.flush()
            stop_time = time.time()
            res[name] = r
            #Python 2
            if sys.version_info >= (3,0):
                if not isinstance(r,str):
                    bp.Print(TFormat.InRed( TFormat().GetIndent() + "Please add a correct return statement in the CheckIntegrity of the module" + name))
                    #raise Exception()
                    r = 'Not OK'
            else:
                if not (isinstance(r,str) or isinstance(r,unicode) ):
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
                 print('-*-*-*-*-*-*> missing CheckIntegrity()??? <*-*-*-*-*-*--'  )
                 if(stopAtFirstError): raise
     except:
        tocheck[submod ] = None
        if(stopAtFirstError): raise


def __tryImport(noduleName,bp,stopAtFirstError):# pragma: no cover

    tocheck = {}
    try :
        __tryImportRecursive(noduleName,tocheck,stopAtFirstError)
    except:
        print("Error loading module '" + noduleName +"'")
        print("This module will not be tested ")

        sys.stdout.flush()
        sys.stderr.flush()
        bp.Print( "Unexpected Error:" + str(sys.exc_info()[0]) )
        traceback.print_exc(file=bp.stdout_)

        if(stopAtFirstError): raise
    return tocheck


def TestAll(modulestotreat=['ALL'], fulloutput=False, stopAtFirstError= False, coverage= False, extraToolsBoxs= None,dryrun=False,profiling=False,browser = False) :# pragma: no cover

    print("")
    print("modulestotreat   : " + str(modulestotreat))
    print("fulloutput       : " + str(fulloutput) )
    print("stopAtFirstError : " + str(stopAtFirstError))
    print("coverage         : " + str(coverage))
    print("browser          : " + str(browser))
    print("profiling        : " + str(profiling))
    print("extraToolsBoxs   : " + str(extraToolsBoxs))
    print("dryrun           : " + str(dryrun))


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

    with PrintBypass() as bp:
        if extraToolsBoxs is not None:
            for tool in extraToolsBoxs:
                tocheck.update(__tryImport(tool,bp,stopAtFirstError))



        if not fulloutput:
            bp.ToSink()
            bp.Print("Sending all the output of the tests to sink")


        if modulestotreat[0] is not  'ALL':
            #bp.Print(str(modulestotreat))
            filtered =  dict((k, v) for k, v in tocheck.items() if all(s in k for s in modulestotreat ) )
            #bp.Print(str(filtered))
            tocheck = filtered

        res = __RunAndCheck(tocheck,bp,stopAtFirstError,dryrun,profiling);

        #bp.Restore()
        #now the restore is done automaticaly


    if coverage :
        cov.stop()
        cov.save()

        # create a temp file
        if browser:
            tempdir = TestTempDir.GetTempPath()
        else:
            tempdir = "./"

        ss = [ ("*"+k.split(".")[-1]+"*") for k in tocheck ]
        cov.html_report(directory = tempdir, include=ss  ,title="Coverage report of "+ " and ".join(extraToolsBoxs) )
        print('Coverage Report in : ' + tempdir +"index.html")
        print(cov.report(show_missing=False))
        if browser:
          import webbrowser
          webbrowser.open(tempdir+"index.html")


    stop_time = time.time()
    bp.Print( "Total Time : %.3f seconds " %  (stop_time -start_time ))

    print("--- End Test ---")
    return res

def CheckIntegrity():
    TestTempDir().GetTempPath()
    TestTempDir().GetTempPath()
    return "Ok"

if __name__ == '__main__':# pragma: no cover
    import sys, getopt
    if len(sys.argv) == 1:
        res = TestAll(modulestotreat=['ALL'],extraToolsBoxs= ["BasicTools"], fulloutput=False,coverage=False)# pragma: no cover
    else:
      try:
          opts, args = getopt.getopt(sys.argv[1:],"hcfsdpe:m:")
      except getopt.GetoptError as e:
          print(e)
          print(Test_Help_String)
          sys.exit(2)

      coverage = False
      fulloutput = False
      stopAtFirstError = False
      extraToolsBoxs = []
      modulestotreat=[]
      dryrun = False
      profiling = False
      browser = False



      for opt, arg in opts:
         if opt == '-h':
             print(Test_Help_String)
             sys.exit()
         elif opt in ("-c"):
            coverage = True
            browser = False
         elif opt in ("-cb"):
            coverage = True
            browser = True
         elif opt in ("-f"):
            fulloutput = True
         elif opt in ("-s"):
            stopAtFirstError = True
         elif opt in ("-d"):
            dryrun = True
         elif opt in ("-p"):
            profiling = True
         elif opt in ("-e"):
            extraToolsBoxs.append(arg)
         elif opt in ("-m"):
            modulestotreat.append(arg)

      if len(modulestotreat) == 0:
         modulestotreat.append("ALL")

      if len(extraToolsBoxs) == 0:
         extraToolsBoxs.append("BasicTools")

      res = TestAll(  modulestotreat=modulestotreat,
                coverage=coverage,
                fulloutput=fulloutput,
                stopAtFirstError= stopAtFirstError,
                extraToolsBoxs=extraToolsBoxs,
                dryrun = dryrun,
                profiling  =profiling,
                browser = browser)
    errors = { x:y  for x,y in res.items() if y.lower() != "ok"}
    oks = { x:y  for x,y in res.items() if y.lower() == "ok"}
    #print(errors)
    #print(len(errors))
    sys.exit(len(errors))



