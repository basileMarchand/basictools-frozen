# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Testing infrastructure for BasicTools extra modules

"""

import sys
import traceback
import time
import getopt
import os

from BasicTools.Helpers.which import which


Test_Help_String = """
python  Tests.py -c -f -s -e <extraModules> -m <moduleFilter>
options :
    -f    Full output for all the test
    -s    Stop at first error
    -e    To test extra Modules (-e can be repeated)
    -m    To filter the output by this string (-m can be repeated)
    -k    Skip test base on string
    -d    Dry run do not execute anything, only show what will be executed
    -c    To activate coverage and generate a html report
          -cb (coverage with browser and local file index.html generated)
    -l    Generation of the coverage (html) report localy (current path)
    -b    Launch browser after the covarage report generation
    -p    Activate profiling
    -v    Activate maximal level of verbosity
    -y    Generate .pyc when inporting modules (default False)
    -L    Output Localy, Use the current path for all the outputs
    -P <dir> Settemp output directory to
"""


def SkipTest(evironementVaribleName):

    if evironementVaribleName in os.environ:
        print("Warning skiping test (environement variable "+str(evironementVaribleName)+" set)")
        return True
    return False


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
    prefix = "BasicTools_Test_Directory_"
    createdOnRam = False
    @classmethod
    def GetTempPath(cls, onRam=None):
        if cls.path is not None and (onRam is None or onRam == cls.createdOnRam):
            return cls.path
        import tempfile

        if onRam is None:
            onRam = False

        if onRam:
            cls.path = tempfile.mkdtemp(prefix=cls.prefix,suffix="_safe_to_delete",dir="/dev/shm/") + os.sep
            cls.createdOnRam = True
        else:
            cls.path = tempfile.mkdtemp(prefix=cls.prefix,suffix="_safe_to_delete") + os.sep
            cls.createdOnRam = False

        cls.__saveTempPath()
        return TestTempDir.path

    #  we cant test this function, because the temp path will be delete
    @classmethod
    def DeleteTempPath(cls):# pragma: no cover
        import shutil
        if cls.path is not None:
            shutil.rmtree(cls.path)
        cls.path = None

    @classmethod
    def OpenTempFolder(cls):# pragma: no cover
        import subprocess
        if os.name == "nt":
            subprocess.Popen('explorer "' + cls.GetTempPath() +'"')
        elif which("nautilus"):
            print(cls.GetTempPath())
            subprocess.Popen(['nautilus',  cls.GetTempPath() ])

    @classmethod
    def SetTempPath(cls,path,create=True):# pragma: no cover
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
        with open(home + os.sep+".BasicToolsTempPath","w") as f:
            f.write("cd " + TestTempDir.path + "\n")
            import stat
            os.chmod(home + os.sep+".BasicToolsTempPath", stat.S_IWUSR | stat.S_IRUSR |stat.S_IXUSR)

def __RunAndCheck(lis,bp,stopAtFirstError,dryrun,profiling):# pragma: no cover

    res = {}
    from BasicTools.Helpers.TextFormatHelper import TFormat

    for name in lis:
        bp.Print(TFormat.InBlue(TFormat.Center( "Running Test " +name ,width=80 ) ))
        if lis[name] is None:
            bp.Print(TFormat().InRed(TFormat().GetIndent() + 'Sub Module "' + name + '" does not have the CheckIntegrity function'))
            r = 'Not OK'
            if stopAtFirstError :
                raise(Exception(TFormat().GetIndent() + 'Sub Module "' + name + '" does not have the CheckIntegrity function'))
        try:
            start_time = time.time()
            stop_time = time.time()
            if dryrun:
                r = "Dry Run "
            else:
                teststr = """
from {} import CheckIntegrity
res = CheckIntegrity()""".format(name)
                local = {}

                if profiling :
                    import cProfile, pstats
                    from io import StringIO
                    pr = cProfile.Profile()
                    pr.enable()
                    exec(teststr,{},local)
                    r = local["res"]
                    pr.disable()
                    s = StringIO()
                    sortby = 'cumulative'
                    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
                    ps.print_stats()
                    print(s.getvalue())
                else:
                    exec(teststr,{},local)
                    r = local["res"]
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
            bp.Print( TFormat.InGreen( f"OK {name} : {stop_time -start_time:.3f} seconds " ) )
        elif r.lower() == 'skip':
            bp.Print(TFormat.InYellow( str(r) + " !!!! " + name )  )
        else:
            bp.Print(TFormat.InRed( str(r) + " !!!! " + name )  )

    return res

def __tryImportRecursive(submod, tocheck, stopAtFirstError, modulestotreat, modulestoskip):


    import importlib
    try:
        sm = importlib.import_module(submod)
    except ImportError as e:
        print(e)
        if stopAtFirstError:
            raise
        sm = None
    except Exception:
        sm = None

    cif = getattr( sm, "CheckIntegrity", None)
    if cif is not None:
        tocheck[submod ] = cif
    else :
        try:
            subsubmods = []
            try:
                listToTest = getattr(sm,"__all__",getattr(sm,"_test",None))
                subsubmods =  [ submod +'.' + x for x in listToTest]
            except Exception as e:
                print(e)
                print('Error Loading module : ' + submod + ' (Current folder'+os.getcwd()+')'  )
                print('-*-*-*-*-*-*> missing _test variable ??? <*-*-*-*-*-*--'  )
                raise

            for subsubmod  in subsubmods:
                if any([x.lower() in subsubmod.lower() for x in modulestoskip]):
                    print("skip module :" +str(subsubmod))
                    continue
                try:
                    __tryImportRecursive(subsubmod,tocheck,stopAtFirstError,modulestotreat,modulestoskip)
                except Exception as e:
                    print(e)
                    print('Error Loading module : ' + subsubmod + ' (Current folder'+os.getcwd()+')'  )
                    print('-*-*-*-*-*-*> missing CheckIntegrity()??? <*-*-*-*-*-*--'  )
                    raise
        except:
            tocheck[submod] = None
            if stopAtFirstError:
                raise


def __tryImport(noduleName, bp, stopAtFirstError, modulestotreat, modulestoskip):# pragma: no cover

    tocheck = {}
    try :
        __tryImportRecursive(noduleName, tocheck, stopAtFirstError, modulestotreat, modulestoskip)
    except:
        print("Error loading module '" + noduleName +"'")
        print("This module will not be tested ")

        sys.stdout.flush()
        sys.stderr.flush()
        bp.Print( "Unexpected Error:" + str(sys.exc_info()[0]) )
        traceback.print_exc(file=bp.stdout_)

        if stopAtFirstError:
            raise
    return tocheck


def TestAll(modulestotreat=['ALL'],modulestoskip=[], fulloutput=False, stopAtFirstError= False, extraToolsBoxs= None,dryrun=False,profiling=False,coverage=None) :# pragma: no cover

    print("")
    print("modulestotreat   : " + str(modulestotreat))
    print("modulestoskip    : " + str(modulestoskip))
    print("fulloutput       : " + str(fulloutput) )
    print("stopAtFirstError : " + str(stopAtFirstError))
    print("coverage         : " + str(coverage))
    print("profiling        : " + str(profiling))
    print("extraToolsBoxs   : " + str(extraToolsBoxs))
    print("dryrun           : " + str(dryrun))


    cov = None
    if coverage["active"]:
        import coverage as modcoverage
        #ss = [ k for k in lis ]
        #cov = coverage.coverage(source=ss, omit=['pyexpat','*__init__.py'])
        cov = modcoverage.coverage(omit=['pyexpat','*__init__.py'])
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
                tocheck.update(__tryImport(tool,bp,stopAtFirstError,modulestotreat,modulestoskip))

        if not fulloutput:
            bp.ToSink()
            bp.Print("Sending all the output of the tests to sink")

        if modulestotreat[0] != 'ALL':
            filtered =  dict((k, v) for k, v in tocheck.items() if any(s in k for s in modulestotreat ) )
            tocheck = filtered

        if len(modulestoskip):
            filtered =  dict((k, v) for k, v in tocheck.items() if not any(s in k for s in modulestoskip ) )
            tocheck = filtered

        res = __RunAndCheck(tocheck,bp,stopAtFirstError,dryrun,profiling)

    if coverage["active"] :
        cov.stop()
        cov.save()

        # create a temp file
        if coverage["localhtml"]:
            #tempdir = "./"
            tempdir = os.getcwd() + os.sep
        else:
            tempdir = TestTempDir.GetTempPath()

        ss = [ ("*"+os.sep.join(k.split("."))+"*") for k in tocheck ]
        cov.html_report(directory = tempdir, include=ss  ,title="Coverage report of "+ " and ".join(extraToolsBoxs) )
        print('Coverage Report in : ' + tempdir +"index.html")
        cov.report( include=ss )
        if coverage["lauchBrowser"]:
            import webbrowser
            webbrowser.open(tempdir+"index.html")


    stop_time = time.time()
    bp.Print( "Total Time : %.3f seconds " %  (stop_time -start_time ))

    print("--- End Test ---")
    #print(tocheck)
    if len(modulestotreat) == 1 and modulestotreat[0] == "ALL" :
        #we verified that all thepython files in the repository are tested
        CheckIfAllThePresentFilesAreTested(extraToolsBoxs,res)
    return res

def CheckIfAllThePresentFilesAreTested(extraToolsBoxs,testedFiles):

    pythonPaths = os.environ.get('PYTHONPATH', os.getcwd()).split(":")

    toignore = ['.git', # git file
                "__pycache__", # python 3
                "__init__.py",# infra
                "setup.py", # compilation
                "BasicTools/docs/conf.py", #documentation
                ]
    testedFiles = testedFiles.keys()
    presentFiles = []
    for pythonPath in pythonPaths:
        for eTB in extraToolsBoxs:
            path = pythonPath + os.sep + eTB

            for dirname, dirnames, filenames in os.walk(path):
                cleandirname = dirname.replace(pythonPath+os.sep,"")#.replace(os.sep,".")


                # print path to all filenames.
                for filename in filenames:
                    if filename in toignore:
                        continue

                    if len(filename) > 3 and filename[-3:] == ".py":
                        if len(cleandirname) == 0:
                            presentFiles.append(filename)
                        else:
                            presentFiles.append(cleandirname+os.sep+filename)

                # Advanced usage:
                # editing the 'dirnames' list will stop os.walk() from recursing into there.
                for dti in toignore:
                    if dti in dirnames:
                        dirnames.remove(dti)

    for tf in testedFiles:
        if tf.replace(".",os.sep)+".py" in presentFiles:
            presentFiles.remove(tf.replace(".",os.sep) +".py")

    if len(presentFiles) > 0 :
        print("files Present in the repository but not tested")
        for i in presentFiles:
            print(i)

def CheckIntegrity():
    TestTempDir().GetTempPath()
    TestTempDir().GetTempPath()
    return "Ok"

def RunTests():
    if len(sys.argv) == 1:
        res = TestAll(modulestotreat=['ALL'],extraToolsBoxs= ["BasicTools"], fulloutput=False,coverage={"active":False})# pragma: no cover
    else:
        try:
            opts, args = getopt.getopt(sys.argv[1:],"hcblfsdpvyLe:m:k:")
        except getopt.GetoptError as e:
            print(e)
            print(Test_Help_String)
            sys.exit(2)

        coverage = False
        fulloutput = False
        stopAtFirstError = False
        extraToolsBoxs = []
        modulestotreat=[]
        modulestoskip=[]
        dryrun = False
        profiling = False
        browser = False
        localhtml = False

        sys.dont_write_bytecode = True

        for opt, arg in opts:
            if opt == '-h':
                print(Test_Help_String)
                sys.exit()
            elif opt in ("-c"):
                coverage = True
                browser = False
            elif opt in ("-l"):
                localhtml = True
                browser = False
            elif opt in ("-b"):
                browser = True
            elif opt in ("-f"):
                fulloutput = True
            elif opt in ("-s"):
                stopAtFirstError = True
            elif opt in ("-d"):
                dryrun = True
            elif opt in ("-v"):
                from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
                myObj = BaseOutputObject()
                myObj.SetGlobalDebugMode(True)
            elif opt in ("-p"):
                profiling = True
            elif opt in ("-e"):
                extraToolsBoxs.append(arg)
            elif opt in ("-m"):
                modulestotreat.append(arg)
            elif opt in ("-k"):
                modulestoskip.append(arg)
            elif opt in ("-P"):
                print('Setting temp output directory to ' + arg)
                TestTempDir.SetTempPath(arg)
            elif opt in ("-y"):
                sys.dont_write_bytecode = False
            elif opt in ("-L"):
                TestTempDir().SetTempPath(os.getcwd())
                from BasicTools.Helpers.Tests import TestTempDir as TestTempDir2
                TestTempDir2().SetTempPath(os.getcwd())

        if len(modulestotreat) == 0:
            modulestotreat.append("ALL")

        if len(extraToolsBoxs) == 0:
            extraToolsBoxs.append("BasicTools")

        res = TestAll(  modulestotreat=modulestotreat,
                        modulestoskip=modulestoskip,
                    coverage={"active":coverage,"localhtml":localhtml,"lauchBrowser": browser},
                    fulloutput=fulloutput,
                    stopAtFirstError= stopAtFirstError,
                    extraToolsBoxs=extraToolsBoxs,
                    dryrun = dryrun,
                    profiling  = profiling
                    )
    errors = {}
    oks = {}
    skipped = {}

    for x,y in res.items():
        if str(y).lower() == "ok":
            oks[x] = y
        elif str(y).lower() == "skip":
            skipped[x] = y
        else:
            errors[x] = y
    print("Number of test OK       : " + str(len(oks)))
    print("Number of skipped test  : " + str(len(skipped)))
    print("Number of test KO       : " + str(len(errors)))

    noks = len(oks)
    nerrors = len(errors)
    if noks ==  0 and nerrors == 0:
        noks = 1
    print(f"Percentage of OK test : {(noks*100.0)/(noks+nerrors):.2f} %")
    return errors

if __name__ == '__main__':# pragma: no cover
    errors = RunTests()
    sys.exit(len(errors))
