# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import subprocess
import os

import numpy as np
import shlex

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject as BOO
import BasicTools.Helpers.ParserHelper as PH

samresExec = "samres"

class SamcefOuputReader(BOO):
    def __init__(self):
        super(SamcefOuputReader,self).__init__()
        self.canHandleTemporal = True
        self.subprocess = None
        self.filename = ""
        self.LCP = "me"
        self.timeToRead = 0.
        self.fieldsToRead = ["U","Element_Stress" ]

    def SetFileName(self, fileName):
        self.PrintDebug(fileName)
        if len(fileName) > 4 and fileName[-4:] == ".fac":
            fileName = fileName[0:-4]
        if len(fileName) > 3 and fileName[-3] == "_":
            self.LCP = fileName[-2:]
            fileName = fileName[0:-3]
        self.filename = fileName

    def SetTimeToRead(self,time= None,timeIndex=None):
        if timeIndex is not None:
            self.timeToRead = self.times[timeIndex]
        elif time is not None:
            self.timeToRead = PH.ReadFloat(time)

    def Open(self):
        if self.subprocess is not None:
            raise(Exception("File already open"))
        casename = self.filename.split(".fac")[0]
        args = [samresExec,"NOM="+casename, "LCP="+self.LCP]
        self.PrintDebug("Opening "+ str(casename) )
        self.subprocess = subprocess.Popen(args, bufsize=1,  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                                           #executable=None,stderr=None, preexec_fn=None, close_fds=True, shell=False, cwd=None, env=None, startupinfo=None, creationflags=0, restore_signals=True, start_new_session=False, pass_fds=(), *, encoding=None, errors=None)
    def OpenIfNeeded(self):
        if self.subprocess == None:
            self.Open()
            return True
        else:
            return False

    def Close(self,needToClose=True):
        if needToClose :
            self.subprocess.terminate()
            self.subprocess.communicate()
            self.subprocess = None

    def ReadMetaData(self):
        needToClose = self.OpenIfNeeded()
        self.times = []
        self._writeCommand('$$GET_VALUE "code 26" " " " "')
        ntimes = PH.ReadInt(self._readResponse())
        for i in range(ntimes):
            data = self._readResponse()
            key,val = shlex.split(data)[0].split()
            fval = PH.ReadFloat(val)
            self.times.append(fval)

        self._skipResponse(3)
        self._skipResponse(ntimes)
        self.PrintDebug(self.times)
        self.Close(needToClose)

    def GetAvailableTimes(self):
        return self.times

    def Read(self):
        self.PrintDebug("start reading")
        from BasicTools.IO.SamcefReader import DatReader
        dreader = DatReader()
        mesh = dreader.Read(fileName=self.filename+".dat")
        self.PrintDebug("reading mesh")

        needClose = self.OpenIfNeeded()

        name3 = ["X","Y","Z"]
        name6 = ["XX","YY","ZZ","XY", "XZ", "YZ"]
        names = {3:name3, 6:name6}

        offsets= {}
        cpt = 0
        for elemtype, data in mesh.elements.items():
            offsets[elemtype] =  cpt
            cpt += data.GetNumberOfElements()

        for fname in self.fieldsToRead:
            self.PrintDebug("Reading "+ fname )
            oid, ofield = self.ReadField(fname)
            if InternalCodes[fname].on == "Nodes":
                field = np.zeros( (mesh.GetNumberOfNodes(),ofield.shape[1]) )
                foid = dreader.filetointernalid
                for i in range(ofield.shape[0]):
                    for j in range(ofield.shape[1]):
                        field[foid[oid[i,j]],j ] = ofield[i,j]
                storage =mesh.nodeFields

            elif InternalCodes[fname].on == "Elements":
                field = np.zeros( (mesh.GetNumberOfElements(),ofield.shape[1]) )
                foid = dreader.filetointernalidElement
                for i in range(ofield.shape[0]):
                    for j in range(ofield.shape[1]):
                        elem,cpt = foid[oid[i,j]]
                        enumber = offsets[elem.elementType]+cpt
                        field[enumber,j ] = ofield[i,j]
                storage =mesh.elemFields

            if field.shape[1] > 3:
                for i in range(field.shape[1]):
                    storage[fname+names[field.shape[1]][i]] = field[:,i]
            else:
                    storage[fname] = field

        self.Close(needClose)

        return mesh

    def ReadField(self,fieldname=None,time=None,timeIndex=None):
        needToClose = self.OpenIfNeeded()
        # file:///softs/samcef/16.1-03/doc/m024/samres-m024.pdf
        #page 25
        self.SetTimeToRead(time=time,timeIndex=timeIndex)

        arg1,arg2,arg3,arg4,arg5,arg6 = InternalCodes[fieldname].GetArgs()
        nbcomp= InternalCodes[fieldname].nbcomp
        arg5 = 'Time '+ str(self.timeToRead)

        cmd = "$$GET_VALUE " + " ".join(['"'+str(x)+'"' for x in [arg1,arg2,arg3,arg4,arg5,arg6]])

        self._writeCommand(cmd)
        NRef = PH.ReadInt(self._readResponse())
        refNames = []
        for i in range(NRef):
            refNames.append(self._readResponse())

        SEL_SIZE = PH.ReadInt(self._readResponse())//nbcomp
        entityname = self._readResponse()

        oid = self._readNumpyArray((SEL_SIZE,nbcomp), dtype=int)
        data = self._readNumpyArray((SEL_SIZE,nbcomp), dtype=float)

        self.Close(needToClose)

        if entityname not in ['"Node %"','"Element %"']:
            raise(Exception("Dont know how to treat->"+ entityname+"<-"))
        return oid,data

    def _writeCommand(self,cmd):
        self.PrintDebug("Command:"+cmd)
        self.subprocess.stdin.write(str(cmd)+os.linesep)
        self.subprocess.stdin.flush()
        self._skipResponse()

    def _skipResponse(self,nblines=1):
        self.PrintDebug(f"Skiping {nblines} lines ")
        for i in range(nblines):
            self._readResponse()

    def _readResponse(self):
        data = self.subprocess.stdout.readline().rstrip()
        if data.find("ERROR") >= 0:
            print(data)
            raise(Exception(data))
        self.PrintDebug("responce: "+ data)
        return data

    def _readNumpyArray(self, shape, dtype=float):
        data = np.loadtxt(self.subprocess.stdout, dtype=dtype, max_rows=int(np.prod(shape)))
        data.shape = shape
        return data

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".fac",SamcefOuputReader)

class SamcefData():
    def __init__(self,arg1=" ",arg2=" ",arg3=" ",arg4=" ",name=" ",nbcomp=1):
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3
        self.arg4 = arg4
        self.nbcomp = nbcomp
    def GetArgs(self):
        return (self.arg1,self.arg2,self.arg3,self.arg4," "," ")

class SamcefDataN(SamcefData):
    def __init__(self,*k,**nk):
        super(SamcefDataN,self).__init__(*k,**nk)
        self.on = "Nodes"

class SamcefDataE(SamcefData):
    def __init__(self,*k,**nk):
        super(SamcefDataE,self).__init__(*k,**nk)
        self.on = "Elements"

class SamcefDataG(SamcefData):
    def __init__(self,*k,**nk):
        super(SamcefDataE,self).__init__(*k,**nk)
        self.on = "Global"

InternalCodes = {}
InternalCodes["Time"] =  SamcefData(arg1='Code 26',name="Time")
InternalCodes["U"] =  SamcefDataN(arg1='Code 163',name="U",nbcomp=3)
InternalCodes["U0"] = SamcefDataN(arg1='Code 163',arg2="Component 1", name="U0")
InternalCodes["U1"] = SamcefDataN(arg1='Code 163',arg2="Component 2", name="U1")
InternalCodes["U2"] = SamcefDataN(arg1='Code 163',arg2="Component 3", name="U2")
InternalCodes["Reaction"] =  SamcefDataN(arg1='Code 221',name="Reaction",nbcomp=3)
InternalCodes["Nodal_Stress"] =  SamcefDataN(arg1='Code 1411',name="Nodal_Stress",nbcomp=6)

InternalCodes["Element_Stress"] =  SamcefDataE(arg1='Code 3411',name="Element_Stress",nbcomp=6)
InternalCodes["Element_Stress_VM"] =  SamcefDataE(arg1='Code 3411',arg2="Von Mises",name="Element_Stress_VM",nbcomp=1)

"""
 "Sdb 4" "[Sdb] Element Mass" " " ".SUPPORT_0102"
 "Sdb 3" "[Sdb] Element Volume" " " ".SUPPORT_0102"
 "Sdb 2" "[Sdb] Element Surface" " " ".SUPPORT_0102"
 "Sdb 1" "[Sdb] Element Length" " " ".SUPPORT_0102"
 "Code 26" "[Code   26] Time" " " ".SUPPORT_0050"
 "Code 163" "[Code  163] Nodal displacements (DX,DY,DZ)" ".TOSCALAR_0003" ".SUPPORT_0003"
 "Code 191" "[Code  191] Nodal rotation (RX,RY,RZ)" ".TOSCALAR_0009" ".SUPPORT_0009"
 "Code 221" "[Code  221] Nodal reaction" ".TOSCALAR_0010" ".SUPPORT_0010"
 "Code 901" "[Code  901] Time" ".TOSCALAR_0901" ".SUPPORT_0901"
 "Code 981" "[Code  981] Kinetic energy" ".TOSCALAR_0981" ".SUPPORT_0981"
 "Code 982" "[Code  982] Potential energy" ".TOSCALAR_0982" ".SUPPORT_0982"
 "Code 984" "[Code  984] External forces work" ".TOSCALAR_0984" ".SUPPORT_0984"
 "Code 985" "[Code  985] Time steps" ".TOSCALAR_0985" ".SUPPORT_0985"
 "Code 986" "[Code  986] Time integration errors" ".TOSCALAR_0986" ".SUPPORT_0986"
 "Code 1305" "[Code 1305] Contact pressure" " " ".SUPPORT_0026"
 "Code 1306" "[Code 1306] Friction stress" " " ".SUPPORT_0026"
 "Code 1307" "[Code 1307] Normal distance" " " ".SUPPORT_0026"
 "Code 1342" "[Code 1342] Sliding velocity direction 1" " " ".SUPPORT_0026"
 "Code 1343" "[Code 1343] Sliding velocity direction 2" " " ".SUPPORT_0026"
 "Code 1344" "[Code 1344] Sliding" " " ".SUPPORT_0026"
 "Code 1411" "[Code 1411] Stress tensor" ".TOSCALAR_0013" ".SUPPORT_0013"
 "Code 1421" "[Code 1421] Strain tensor" ".TOSCALAR_0013" ".SUPPORT_0013"
 "Code 1431" "[Code 1431] 2D linear stress tensor" ".TOSCALAR_0015" ".SUPPORT_0015"
 "Code 1437" "[Code 1437] Normal force" ".TOSCALAR_0015" ".SUPPORT_0015"
 "Code 1445" "[Code 1445] 2D strain tensor (upper skin)" ".TOSCALAR_0015" ".SUPPORT_0015"
 "Code 1446" "[Code 1446] 2D strain tensor (lower skin)" ".TOSCALAR_0015" ".SUPPORT_0015"
 "Code 1524" "[Code 1524] Force vector" ".TOSCALAR_0025" ".SUPPORT_0025"
 "Code 1525" "[Code 1525] Moment vector" ".TOSCALAR_0025" ".SUPPORT_0025"
 "Code 3408" "[Code 3408] 2D Cauchy stress tensor (low. GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3409" "[Code 3409] 2D Biot stress tensor (low. GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3411" "[Code 3411] Stress tensor" ".TOSCALAR_0014" ".SUPPORT_0014"
 "Code 3418" "[Code 3418] 2D Cauchy strain tensor (low. GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3419" "[Code 3419] 2D Biot strain tensor (low. GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3431" "[Code 3431] 2D linear stress tensor" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3448" "[Code 3448] 2D Cauchy stress tensor (middle GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3449" "[Code 3449] 2D Biot stress tensor (middle GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3452" "[Code 3452] Cauchy strain tensor" ".TOSCALAR_0014" ".SUPPORT_0014"
 "Code 3478" "[Code 3478] 2D Cauchy stress tensor (upper GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3479" "[Code 3479] 2D Biot stress tensor (upper GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3488" "[Code 3488] 2D Cauchy strain tensor (upper GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 3489" "[Code 3489] 2D Biot strain tensor (upper GP)" ".TOSCALAR_0006" ".SUPPORT_0006"
 "Code 9821" "[Code 9821] Resultants - 3 components" ".TOSCALAR_9821" ".SUPPORT_9821"
 "Code 9822" "[Code 9822] Number of plastic elements" ".TOSCALAR_9822" ".SUPPORT_9822"
 "Code 9823" "[Code 9823] Number of damaged elements" ".TOSCALAR_9823" ".SUPPORT_9823"
 "Code 9824" "[Code 9824] Number of nodes in contact" ".TOSCALAR_9824" ".SUPPORT_9824"
 "Code 9850" "[Code 9850] Number of Cycles" ".TOSCALAR_9850" ".SUPPORT_9850"
"""

def ReadFieldFromDataBase(fieldname,field, time=None):
    reader = SamcefOuputReader()
    reader.SetFileName(fieldname)
    reader.Open()
    return  reader.ReadField(field,time=time)

def CheckIntegrity():
    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("SAMCEF_NO_FAIL"): return "ok"

    # no .fac in the database for test for the moment
    return "ok"
    tempfileName = "XXX.fac"

    reader = SamcefOuputReader()
    reader.SetGlobalDebugMode(False)
    reader.SetFileName(tempfileName)
    reader.Open()
    data = reader.ReadMetaData()
    data = reader.ReadField("U0")
    reader.Close()

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
