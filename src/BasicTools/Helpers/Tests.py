# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Testing infrastructure for BasicTools extra modules

"""
from typing import Callable, List, Tuple
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
    -l    Generation of the coverage (html) report locally (current path)
    -b    Launch browser after the coverage report generation
    -p    Activate profiling
    -t    use mypy to check the typing of every module
    -v    Activate maximal level of verbosity
    -y    Generate .pyc when importing modules (default False)
    -L    Output Locally, Use the current path for all the outputs
    -P <dir> Set temporary output directory to
    -g    to use the Froze decorator during testing
"""


def SkipTest(environnementVariableName:str) -> bool:
    """chef if a environnement variable is present to help the user to skip test

    Parameters
    ----------
    environnementVariableName : str
        the name of a environnement variable

    Returns
    -------
    bool
        true if the variable is present
        false if absent
    """

    if environnementVariableName in os.environ:
        print("Warning skipping test (environnement variable "+str(environnementVariableName)+" set)")
        return True
    return False

def RunListOfCheckIntegrities(  toTest: List[Callable[[bool],str] ], GUI:bool = False ) -> str:
    """ Execute all the functions in the list. Stop at the first error.
        the functions must return "ok", "skip" to be treated as successful function

    Parameters
    ----------
    toTest : List[Callable[[bool],str] ]
        a list of functions
    GUI : bool, optional
        bool argument passed to every function , by default False

    Returns
    -------
    str
        "ok" if all functions are executed correctly (the return of every function is "ok")
        "error..."  at the first error encountered
    """
    for func in toTest:
        print("running test : " + str(func))
        res = func(GUI)
        if any([str(res).lower().startswith(x) for x  in ["ok","skip"]]) :
            continue
        return f"error in {func} res"
    return "ok"


def GetUniqueTempFile(suffix:str="", prefix:str='tmp') -> Tuple[int,str]:
    """Create a unique temporary file on the temp directory

    Parameters
    ----------
    suffix : str, optional
        If 'suffix' is not None, the file name will end with that suffix, otherwise there will be no suffix., by default ""
    prefix : str, optional
        If 'prefix' is not None, the file name will begin with that prefix, otherwise a default prefix is used., by default 'tmp'

    Returns
    -------
    Tuple[int,str]
        please read the doc of module tempfile.mkstemp
    """

    #solution from
    # https://stackoverflow.com/questions/2961509/python-how-to-create-a-unique-file-name

    import tempfile
    (fd, filename) = tempfile.mkstemp(suffix=suffix, prefix=prefix,dir=TestTempDir.GetTempPath())
    return (fd,filename)


def WriteTempFile(filename:str, content:str=None, mode:str="w") ->  str:
    """Create and fill a temporary file on the temporary folder

    Parameters
    ----------
    filename : str
        the name of the file to be created
    content : str, optional
        the content of the file, by default None
    mode : str, optional
        read python "open" function documentation for more detail, by default "w"

    Returns
    -------
    str
        the name of the filename (with the path attached)

    Raises
    ------
    Exception
        if error during the file creation
    """

    pFile = TestTempDir.GetTempPath() + filename
    with open(pFile, mode) as f:
        if content is not None:
            f.write( content)
        return pFile

    raise Exception(f"Unable ot create file :{pFile}" )

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

    #very useful in combination of a alias
    #alias cdtmp='source ~/.BasicToolsTempPath'
    @classmethod
    def __saveTempPath(cls):
        from os.path import expanduser
        home = expanduser("~")
        with open(home + os.sep+".BasicToolsTempPath","w") as f:
            f.write("cd " + TestTempDir.path.replace("\\","/") + "\n")
            import stat
            os.chmod(home + os.sep+".BasicToolsTempPath", stat.S_IWUSR | stat.S_IRUSR |stat.S_IXUSR)

def __RunAndCheck(lis,bp,stopAtFirstError,dryRun,profiling,typing):# pragma: no cover

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
            if dryRun:
                r = "Dry Run "
            else:
                testStr = """
from {} import CheckIntegrity
res = CheckIntegrity()""".format(name)
                local = {}

                if profiling :
                    import cProfile, pstats
                    from io import StringIO
                    pr = cProfile.Profile()
                    pr.enable()
                    exec(testStr,{},local)
                    r = local["res"]
                    pr.disable()
                    s = StringIO()
                    sortby = 'cumulative'
                    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
                    ps.print_stats()
                    print(s.getvalue())
                else:
                    exec(testStr,{},local)
                    r = local["res"]
            sys.stdout.flush()
            sys.stderr.flush()
            stop_time = time.time()
            res[name] = r
            if not isinstance(r,str):
                bp.Print(TFormat.InRed( TFormat().GetIndent() + "Please add a correct return statement in the CheckIntegrity of the module" + name))
                #raise Exception()
                r = 'Not OK'

            if typing:
                import subprocess
                subprocess.call(["mypy","--ignore-missing-imports","-m",name])

        except UserWarning as e :
            sys.stdout.flush()
            sys.stderr.flush()
            bp.Print( "Unexpected Warning:" + str(sys.exc_info()[0]) )
            res[name] = "error"
            traceback.print_exc(file=bp.stdout_)
            r = 'Not OK'
        except KeyboardInterrupt as e:
            sys.stdout.flush()
            sys.stderr.flush()
            bp.Print( "User interruption " )
            res[name] = "interrupted"
            r = 'Interrupted'
            bp.Restore()
            raise
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


        if r.lower().startswith('ok'):
            bp.Print( TFormat.InGreen( f"OK {name} : {stop_time -start_time:.3f} seconds " ) )
        elif r.lower().startswith('skip'):
            bp.Print(TFormat.InYellow( str(r) + " !!!! " + name )  )
        else:
            bp.Print(TFormat.InRed( str(r) + " !!!! " + name )  )

    return res

def __tryImportRecursive(subModule, toCheck, stopAtFirstError, modulesToTreat, modulesToSkip):


    import importlib
    try:
        sm = importlib.import_module(subModule)
    except ImportError as e:
        print(e)
        if stopAtFirstError:
            raise
        sm = None
    except Exception:
        sm = None

    cif = getattr( sm, "CheckIntegrity", None)
    if cif is not None:
        toCheck[subModule ] = cif
    else :
        try:
            subSubModules = []
            try:
                listToTest = getattr(sm,"__all__",getattr(sm,"_test",None))
                subSubModules =  [ subModule +'.' + x for x in listToTest]
            except Exception as e:
                print(e)
                print('Error Loading module : ' + subModule + ' (Current folder'+os.getcwd()+')'  )
                print('-*-*-*-*-*-*> missing _test variable ??? <*-*-*-*-*-*--'  )
                raise

            for subSubModule  in subSubModules:
                if any([x.lower() in ("."+subSubModule.lower()+".") for x in modulesToSkip]):
                    print("skip module :" +str(subSubModule))
                    continue
                try:
                    __tryImportRecursive(subSubModule,toCheck,stopAtFirstError,modulesToTreat,modulesToSkip)
                except Exception as e:
                    print(e)
                    print('Error Loading module : ' + subSubModule + ' (Current folder'+os.getcwd()+')'  )
                    print('-*-*-*-*-*-*> missing CheckIntegrity()??? <*-*-*-*-*-*--'  )
                    raise
        except:
            toCheck[subModule] = None
            if stopAtFirstError:
                raise


def __tryImport(noduleName, bp, stopAtFirstError, modulesToTreat, modulesToSkip):# pragma: no cover

    toCheck = {}
    try :
        __tryImportRecursive(noduleName, toCheck, stopAtFirstError, modulesToTreat, modulesToSkip)
    except:
        print("Error loading module '" + noduleName +"'")
        print("This module will not be tested ")

        sys.stdout.flush()
        sys.stderr.flush()
        bp.Print( "Unexpected Error:" + str(sys.exc_info()[0]) )
        traceback.print_exc(file=bp.stdout_)

        if stopAtFirstError:
            raise
    return toCheck


def TestAll(modulesToTreat=['ALL'],modulesToSkip=[], fullOutput=False, stopAtFirstError= False, extraToolsBoxes= None, dryRun=False, profiling=False, coverage=None, typing=False) :# pragma: no cover

    print("")
    print("modulesToTreat   : " + str(modulesToTreat))
    print("modulesToSkip    : " + str(modulesToSkip))
    print("fullOutput       : " + str(fullOutput) )
    print("stopAtFirstError : " + str(stopAtFirstError))
    print("coverage         : " + str(coverage))
    print("typing           : " + str(typing))
    print("profiling        : " + str(profiling))
    print("extraToolsBoxes   : " + str(extraToolsBoxes))
    print("dryRun           : " + str(dryRun))


    cov = None
    if coverage["active"]:
        import coverage as modcoverage
        #ss = [ k for k in lis ]
        cov = modcoverage.Coverage(omit=['pyexpat','*__init__.py'])
        cov.start()

    # calls to print, ie import module1
    from BasicTools.Helpers.PrintBypass import PrintBypass

    print("Running Tests : ")
    start_time = time.time()
    print("--- Begin Test ---")

    toCheck = {}

    with PrintBypass() as bp:
        if extraToolsBoxes is not None:
            for tool in extraToolsBoxes:
                toCheck.update(__tryImport(tool,bp,stopAtFirstError,modulesToTreat,modulesToSkip))

        if not fullOutput:
            bp.ToSink()
            bp.Print("Sending all the output of the tests to sink")

        if modulesToTreat[0] != 'ALL':
            filtered =  dict((k, v) for k, v in toCheck.items() if any(s in k for s in modulesToTreat ) )
            toCheck = filtered

        if len(modulesToSkip):
            filtered =  dict((k, v) for k, v in toCheck.items() if not any(s in k for s in modulesToSkip ) )
            toCheck = filtered

        res = __RunAndCheck(toCheck,bp,stopAtFirstError,dryRun,profiling,typing)

    if coverage["active"] :
        cov.stop()
        cov.save()

        # create a temp file
        if coverage["localHtml"]:
            #tempdir = "./"
            tempdir = os.getcwd() + os.sep
        else:
            tempdir = TestTempDir.GetTempPath()

        ss = [ ("*/"+os.sep.join(k.split("."))+"*") for k in toCheck ]
        cov.html_report(directory = tempdir, include=ss  ,title="Coverage report of "+ " and ".join(extraToolsBoxes) )
        print('Coverage Report in : ' + tempdir +"index.html")
        cov.report( include=ss )
        if coverage["launchBrowser"]:
            import webbrowser
            webbrowser.open(tempdir+"index.html")


    stop_time = time.time()
    bp.Print( "Total Time : %.3f seconds " %  (stop_time -start_time ))

    print("--- End Test ---")
    if len(modulesToTreat) == 1 and modulesToTreat[0] == "ALL" :
        #we verified that all the python files in the repository are tested
        CheckIfAllThePresentFilesAreTested(extraToolsBoxes,res)
    return res

def CheckIfAllThePresentFilesAreTested(extraToolsBoxes,testedFiles):

    pythonPaths = os.environ.get('PYTHONPATH', os.getcwd()).split(":")

    toIgnore = ['.git', # git file
                "__pycache__", # python 3
                "__init__.py",# infra
                "setup.py", # compilation
                "BasicTools/docs/conf.py", #documentation
                ]
    testedFiles = testedFiles.keys()
    presentFiles = []
    for pythonPath in pythonPaths:
        for eTB in extraToolsBoxes:
            path = pythonPath + os.sep + eTB

            for dirname, dirnames, filenames in os.walk(path):
                cleanDirName = dirname.replace(pythonPath+os.sep,"")#.replace(os.sep,".")


                # print path to all filenames.
                for filename in filenames:
                    if filename in toIgnore:
                        continue

                    if len(filename) > 3 and filename[-3:] == ".py":
                        if len(cleanDirName) == 0:
                            presentFiles.append(filename)
                        else:
                            presentFiles.append(cleanDirName+os.sep+filename)

                # Advanced usage:
                # editing the 'dirnames' list will stop os.walk() from recursing into there.
                for dti in toIgnore:
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

def RunTests() -> int:
    """Base function to run all tests, all the option are captured from the command line arguments.

    Returns
    -------
    int
        return the number of failed tests
    """
    if len(sys.argv) == 1:
        res = TestAll(modulesToTreat=['ALL'],extraToolsBoxes= ["BasicTools"], fullOutput=False,coverage={"active":False},typing=False)# pragma: no cover
    else:
        try:
            opts, args = getopt.getopt(sys.argv[1:],"thcblfsdpvyLe:m:k:P:g")
        except getopt.GetoptError as e:
            print(e)
            print(Test_Help_String)
            sys.exit(2)

        coverage = False
        typing = False
        fullOutput = False
        stopAtFirstError = False
        extraToolsBoxes = []
        modulesToTreat=[]
        modulesToSkip=[]
        dryRun = False
        profiling = False
        browser = False
        localHtml = False

        sys.dont_write_bytecode = True

        for opt, arg in opts:
            if opt == '-h':
                print(Test_Help_String)
                sys.exit()
            elif opt in ("-c"):
                coverage = True
                browser = False
            elif opt in ("-t"):
                typing = True
                os.environ["MYPYPATH"] = os.environ["PYTHONPATH"]
            elif opt in ("-l"):
                localHtml = True
                browser = False
            elif opt in ("-b"):
                browser = True
            elif opt in ("-f"):
                fullOutput = True
            elif opt in ("-s"):
                stopAtFirstError = True
            elif opt in ("-d"):
                dryRun = True
            elif opt in ("-v"):
                from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
                myObj = BaseOutputObject()
                myObj.SetGlobalDebugMode(True)
                myObj.SetVerboseLevel(2)
            elif opt in ("-p"):
                profiling = True
            elif opt in ("-e"):
                extraToolsBoxes.append(arg)
            elif opt in ("-m"):
                modulesToTreat.append(arg)
            elif opt in ("-k"):
                modulesToSkip.append(arg)
            elif opt in ("-P"):
                print('Setting temp output directory to ' + arg)
                TestTempDir.SetTempPath(arg)
            elif opt in ("-y"):
                sys.dont_write_bytecode = False
            elif opt in ("-L"):
                TestTempDir().SetTempPath(os.getcwd())
                from BasicTools.Helpers.Tests import TestTempDir as TestTempDir2
                TestTempDir2().SetTempPath(os.getcwd())
            elif  opt in ("-g"):
                import BasicTools.Helpers.BaseOutputObject as BaseOutputObject
                BaseOutputObject.useFroze_itDecorator = True
                print("BaseOutputObject.useFroze_itDecorator =", BaseOutputObject.useFroze_itDecorator)



        if len(modulesToTreat) == 0:
            modulesToTreat.append("ALL")

        if len(extraToolsBoxes) == 0:
            extraToolsBoxes.append("BasicTools")

        res = TestAll(  modulesToTreat=modulesToTreat,
                        modulesToSkip=modulesToSkip,
                    coverage={"active":coverage,"localHtml":localHtml,"launchBrowser": browser},
                    fullOutput=fullOutput,
                    stopAtFirstError= stopAtFirstError,
                    extraToolsBoxes=extraToolsBoxes,
                    dryRun = dryRun,
                    profiling  = profiling,
                    typing=typing
                    )
    errors = {}
    oks = {}
    skipped = {}

    for x,y in res.items():
        if str(y).lower().startswith("ok"):
            oks[x] = y
        elif str(y).lower().startswith("skip")  :
            skipped[x] = y
        else:
            errors[x] = y
    print("Number of test OK       : " + str(len(oks)))
    print("Number of skipped test  : " + str(len(skipped)))
    print("Number of test KO       : " + str(len(errors)))

    nbOks = len(oks)
    nbErrors = len(errors)
    if nbOks ==  0 and nbErrors == 0:
        nbOks = 1
    print(f"Percentage of OK test : {(nbOks*100.0)/(nbOks+nbErrors):.2f} %")
    return errors

if __name__ == '__main__':# pragma: no cover
    errors = RunTests()
    sys.exit(len(errors))
