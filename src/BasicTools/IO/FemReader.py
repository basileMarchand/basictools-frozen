# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Fem file reader
"""

import numpy as np

import BasicTools.Containers.UnstructuredMesh as UM
import BasicTools.Containers.ElementNames as ElementNames

from BasicTools.IO.ReaderBase import ReaderBase


def ReadFem(fileName=None,string=None):
    """Function API for reading a Fem result file

    Parameters
    ----------
    fileName : str, optional
        name of the file to be read, by default None
    string : str, optional
        data to be read as a string instead of a file, by default None

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    obj = FemReader()
    if fileName is not None:
        obj.SetFileName(fileName)

    if string is not None:
        obj.SetStringToRead(string)

    return obj.Read()

keysToIgnore = [ "CORD1R",
                "DESGLB",
                "MAT1",
                "PSOLID",
                "MAT1",
                "LABEL",
                "FORMAT",
                "OUTPUT",
                "ELFORCE",
                "ESE(H3D)",
                "GPSTRESS(H3D,",
                "SPCFORCE(SPARSE)",
                "STRESS",
                "OUTPUT,",
                "",
                ]

keysToIgnoreSubcase = [ "",
                    ]

keysToIgnoreBulk = [ "PSOLID",
                "DOPTPRM",
                "DOPTPRM,",
                "GAPPRM",
                "GAPPRM,",
                "DCONADD",
                "PBUSH",
                "DCONSTR",
                "LOAD",
                "SPCADD",
                "LOADADD",
                "MAT1",
                "OUTPUT,",
                ]

class ignored():
    def __init__(self,ignored):
        pass

class NastranLineParcer(object):
    def __init__(self):
        self.fields = []
        self.data = {}
        self.structured = False

    def Parse(self,line,file=None):
        #k = GetField(line,0)
        #if k != self.fields[0][0]:
        #    raise(Exception("Keyword does not match : " + str(k) + " != "+str(self.fields[0][0])))

        for i in range(1,len(self.fields)):
            f = self.GetField(line,i)
            name = self.fields[i][0]
            types = self.fields[i][1]
            self.data[name] = self.ParseTypes(f,types)
        if self.GetField(line,9) == "+" :
            line = file.getLine()
            self.ContinueParsing(line)

    def GetFields(self,line,fn,ffn):
        return (self.GetField(line,x) for x in range(fn,ffn) )

    def GetField(self,line,fn):
        if self.structured:
            if (fn+1)*8 > len(line):
                line += (8*(fn+1) - len(line))*" "
            return line[8*fn:8*(fn+1)].lstrip().rstrip()
        else:
            return line.split()[fn]

    def ParseType(self,field,typ):
        if typ == float:
            return self.ParseFloat(field)
        else:
            return typ(field)

    def ParseTypes(self,fields,ts):
        res = []
        for i,t in enumerate(ts):
            f = self.GetField(fields,i)
            res.append(self.ParseType(f, t))
        return res

    def ParseFloat(self,line):
            line += (8-len(line))*' '
            if line[6] == "-"  and line[5] != " ":
                    line = line[0:6] + "e" + line[6:8]
            elif line[5] == "-"  and line[4] != " ":
                    line = line[0:5] + "e" + line[5:8]
            elif line[4] == "-"  and line[3] != " ":
                    line = line[0:4] + "e" + line[4:8]
            elif line[3] == "-"  and line[2] != " ":
                    line = line[0:3] + "e" + line[3:8]
            elif line[2] == "-"  and line[1] != " ":
                    line = line[0:2] + "e" + line[2:8]
            try:
                    res = float(line)
            except:# pragma: no cover
                    print("string is ->",line)
                    raise
            return res
class FORCE(NastranLineParcer):
    def __init__(self):
        super(FORCE,self).__init__()
        self.structured = True
        self.fields = [("KeyWord",str), ("SIG",int),("G",int),("CID",(None,int)),("F",float),("N1",float),("N2",float),("N3",float)]

class LOAD(NastranLineParcer):
    def __init__(self):
        super(LOAD,self).__init__()
        self.structured = True
        self.fields = [("KeyWord",str), ("S",float),("S1",float),("L1",float)]

class IgnoredOneLine(NastranLineParcer):
    pass
class IgnoredOneMultiline(NastranLineParcer):
    pass

class FemReader(ReaderBase):
    """Fem Reader class
    """
    def __init__(self,fileName = None):
        super(FemReader,self).__init__()
        self.commentChar = "$"
        self.SetFileName(fileName)
        self.ignoreNotTreated = True

    def GetField(self,line,number):
        return line[8*number:8*(number+1)]

    def ReadCleanLine(self,withError=False):
        res = super(FemReader,self).ReadCleanLine(withError=withError)
        import re
        ## add space after each comma
        return re.sub(r'([,/\\]+)(?![ ])', r'\1 ', res)

    def Read(self):
        """Function that performs the reading of a Fem file

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        self.StartReading()

        NLP = NastranLineParcer()

        res = UM.UnstructuredMesh()
        resdic = {}
        resdic["CordinateSystems"] = {}
        resdic["Objective"] = None
        resdic["Optim Functions"] = {}
        resdic["Cases"] = {}
        ids = []
        xs =[]
        ys = []
        zs = []
        filetointernalid = {}
        filetointernalidElement = {}

        cpt =0

        elements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
        elements.Reserve(1000000)
        need_to_read = True
        while(True):
            if need_to_read :
                line = self.ReadCleanLine()
            need_to_read = True
            if line is None:
                break

            l = line.split()
            key = l[0]
            #print(key)
            if key in keysToIgnore and self.ignoreNotTreated :
                while(True):
                    line = self.ReadCleanLine()
                    if line[0] != "+":
                        break
                need_to_read = False
                continue

            if key == 'SUBCASE' :
            #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-38938C32-4019-4DDA-AE94-2987C8355669-htm.html
                #subcase = int(self.GetField(line,1))
                subcase = int(l[1])
                data = {}

                while(True):
                    line = self.ReadCleanLine()
                    key = line.split()[0]

                    if key in keysToIgnoreSubcase and self.ignoreNotTreated :
                        while(True):
                            line = self.ReadCleanLine()
                            if line[0] != "+":
                                break
                        need_to_read = False
                        continue


                    #print("tt" +line)
                    if line.find("LABEL") > -1:
                    #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-99CAEA1D-352E-4CDA-B74B-C58E4D12D020-htm.html?st=label
                        data['label'] = line.split("LABEL")[1].strip()
                        continue
                    if key == 'ANALYSIS':
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-E7589B2F-C046-4B34-A668-B4E796172645-htm.html
                        #print(line)
                        data['analisys'] = line.split()[1].strip()
                        continue
                    if key == 'SPC' :
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-A456B121-CEB2-48B6-9EBA-7F19EF250AA5-htm.html
                        data['SPC'] = int(line.split("=")[1].strip())
                        continue
                    if key == 'LOAD':
                        # External Static Load Set Selection
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-02DF840A-A5EB-4F04-AB79-DC82A5691105-htm.html
                        data['LOAD'] = int(line.split("=")[1].strip())
                        continue

                    if key == 'MPCFORCE'  :
                        print("MPCFORCE intensionally ignored")
                        continue
                        # OUTPOUT
                        #  Requests multipoint constraint force vector output
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-97B55735-D7B3-4396-A161-BC6B8B694B98-htm.html?st=MPCFORCES

                    if key == 'STRAIN'  :
                        print("STRAIN intensionally ignored")
                        continue
                        # OUTPOUT
                        #  Element Strain Output Request
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-5DE89858-947A-4B35-B244-03BEE2CA779C-htm.html
                    if key == 'DESSUB'  :
                        # Subcase Level Design Constraints Selection
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-093E2B66-D5C4-4C2E-8494-7E6AD36F8B01-htm.html
                        data['DESSUB'] = int(line.split("=")[1].strip())
                        continue
                    break
                resdic['Cases'][subcase] = data
                need_to_read = False
                continue

            if line[0:6] == 'DESOBJ' :
                vec = line.split("=")
                resdic["Objective"] = {"id":int(vec[1]),"type":vec[0].split("(")[1].split(")")[0]}
                continue

            if key == "SET":
                need_to_read = False
                l = line.split()
                key = l[0]
                self.PrintDebug("Reading Sets ")
                setNumber = l[1]
                data = []

                if l[2] != "=":
                    raise Exception("Error in the file format")

                line = " ".join(l[3:])

                while(True):
                    if need_to_read :
                        line = self.ReadCleanLine().strip()
                        need_to_read = True
                    if line is None:
                        break

                    units = [x.strip() for x in  line.split(",")]
                    for u in units:
                        if len(u) == 0:
                            continue
                        if "THRU" in u:
                            partition = u.split()
                            begin = int(partition[0])
                            end = int(partition[-1])
                            data.extend(range(begin,end+1) )
                        else:
                            idd = int(u)
                            data.append(idd)
                    need_to_read = True
                    if line[-1] != ",":
                        break
                continue


            if line == "BEGIN BULK" :
                need_to_read = True
                while(True):
                    if need_to_read :
                        line = self.ReadCleanLine()
                    need_to_read = True
                    if line is None:
                        break

                    l = line.split()
                    key = l[0]

                    if key in keysToIgnoreBulk and self.ignoreNotTreated :
                        while(True):
                            line = self.ReadCleanLine()
                            if line[0] != "+":
                                break
                        need_to_read = False
                        continue

                    if key == 'DOPTPRM'  :
                        resdic["DOPTPRM"] = {self.GetField(line,1):NLP.ParseFloat(self.GetField(line,2))}
                        continue
                    if key == 'DTPL' :
                        print("DTPL intensionally ignored")
                        while(True):
                            line = self.ReadCleanLine()
                            if line[0] != " ":
                                break
                        need_to_read = True
                        continue


                    if key == "PBUSH":
                        #PBUSH is in the keysToIgnoreBulk list
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2018/ENU/NSTRN-Reference/files/GUID-47AAEE71-02F5-4DB2-9817-39F5EE1B0CE8-htm.html
                        NLP.structured = True
                        PID,K  = NLP.GetFields(line,1,3)
                        if K != "K" :
                            print("--",K,"--")
                            raise Exception("Error in the file format")

                        (K1,K2,K3) = NLP.GetFields(line,3,6)
                        print( (K1,K2,K3) )
                        need_to_read = True

                        continue

                    if key == "CBUSH":
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2018/ENU/NSTRN-Reference/files/GUID-86B41C9D-41DB-4664-AE0D-B4B55981F183-htm.html
                        NLP.structured = True

                        EID,PID,GA,GB  = NLP.GetFields(line,1,5)
                        GA = int(GA)
                        GB = int(GB)

                        #for the moment the CBUSH are converted to ntags
                        #elements = res.GetElementsOfType(ElementNames.Bar_2)
                        #conn = [filetointernalid[xx] for xx in  [GA, GB] ]
                        #localid = elements.AddNewElement(conn,EID)
                        #filetointernalidElement[EID] = (elements,localid-1)
                        need_to_read = True

                        tagGA = res.GetNodalTag("CBUSH_"+PID+"_GA")
                        tagGB = res.GetNodalTag("CBUSH_"+PID+"_GB")
                        tagGA.AddToTag(filetointernalid[GA])
                        tagGB.AddToTag(filetointernalid[GB])

                        continue

                    if key == 'DRESP1' :
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-4C44D4BF-E70D-4548-8BED-D4C0497E5479-htm.html?st=DRESP1
                        idd = int(self.GetField(line,1))
                        userlabel = self.GetField(line,2).strip()
                        RTYPE = self.GetField(line,3).strip()
                        REGION = self.GetField(line,5).strip()
                        ATTA = self.GetField(line,6).strip()
                        ATTB = self.GetField(line,7).strip()
                        ATTC = self.GetField(line,8).strip()
                        bulkdata = []
                        line = self.ReadCleanLine()
                        while line[0] == "+":
                            try:
                                #print("--" + line)
                                for x in [1,2,3,4,5,6,7,8]:
                                    bulkdata.append(int(self.GetField(line,x)))
                                line = self.ReadCleanLine()
                            except Exception as inst:
                                #print(inst)
                                line = self.ReadCleanLine()
                                break
                        need_to_read    = False
                        #print(bulkdata)
                        continue

                    if line[0:5] == 'PARAM' and self.ignoreNotTreated : continue
                    if key == 'CORD2R' :
                        # from
                        # https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2018/ENU/NSTRN-Reference/files/GUID-F048B85B-79B9-49AD-94B4-87CC69122C3B-htm.html
                        idd = int(self.GetField(line,1))
                        floats = [NLP.ParseFloat(self.GetField(line,x)) for x in [3,4,5,6,7,8]]
                        line = self.ReadCleanLine()
                        floats.extend([NLP.ParseFloat(self.GetField(line,x)) for x in [1,2,3]] )
                        from BasicTools.Linalg.Transform import Transform
                        orient = Transform()
                        orient.system = "CORD2R"

                        A = np.array(floats[0:3])
                        B = np.array(floats[3:6])
                        C = np.array(floats[6:9])
                        z = (B-A)/np.linalg.norm(B-A)
                        x = (C-A)- z*np.dot(z,C-A)
                        y = np.cross(z,x)

                        orient.SetOffset(A)
                        orient.SetFirst(x)
                        orient.SetSecond(y)
                        resdic["CordinateSystems"][idd] = orient

                        continue
                    if key == 'CORD1R' and self.ignoreNotTreated : continue
                    if key == 'RBE2' :
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2017/ENU/NSTRN-Reference/files/GUID-9C8F3EE5-D900-487A-88E7-9502838E0A3F-htm.html?st=RBE2
                        NLP.structured = True
                        EID, GN,CM = NLP.GetFields(line,1,4)
                        EID = int(EID)
                        GN = int(GN)

                        elements = res.GetElementsOfType(ElementNames.Bar_2)
                        tag = elements.tags.CreateTag("RBE2_"+str(EID),False)

                        oidd = filetointernalid[GN]

                        start =  4
                        while(True):
                            line = self.ReadCleanLine()
                            if line[0] != "+":
                                break
                            for i in range(start,9):
                                data = line[8*i:8*(i+1)].strip()
                                if len(data) ==0:
                                    break
                                GMi = int(data)
                                localid = elements.AddNewElement([oidd, filetointernalid[GMi] ],-1)
                                tag.AddToTag(localid-1)


                            start =1
                        need_to_read = False
                        continue
                    if key == 'RBE3' :
                        # from https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2017/ENU/NSTRN-Reference/files/GUID-69695F7D-3B6E-4183-8BB0-838E8A014CD2-htm.html?st=RBE3

                        EID, _, REFGRID,REFC,WTi,C1,G11,G12  = NLP.GetFields(line,1,9)
                        EID =int(EID)
                        REFGRID = int(REFGRID)
                        REFC = int(REFC)
                        ##Reading weights
                        WTi = float(WTi)
                        C1 = int(C1)
                        G11 = int(G11)
                        G12 = int(G12)

                        oidd = filetointernalid[REFGRID]
                        elements = res.GetElementsOfType(ElementNames.Bar_2)
                        tag = elements.tags.CreateTag("RBE3_"+str(EID)+"_"+str(REFC))
                        #localid = elements.AddNewElement([oidd, filetointernalid[P1] ],-1)
                        #tag.AddToTag(localid-1)
                        P2 = int(line[8*7:8*8])
                        localid = elements.AddNewElement([oidd, filetointernalid[P2]],-1)
                        tag.AddToTag(localid-1)
                        P3 = int(line[8*8:8*9])
                        localid = elements.AddNewElement([oidd, filetointernalid[P3]],-1)
                        tag.AddToTag(localid-1)

                        while(True):
                            line = self.ReadCleanLine()
                            if line[0] != "+":
                                break
                            for i in range(1,9):
                                data = line[8*i:8*(i+1)].strip()
                                if len(data) == 0:
                                    break
                                P1 = int(data)
                                localid = elements.AddNewElement([oidd ,filetointernalid[P1]],-1)
                                tag.AddToTag(localid-1)


                        need_to_read = False
                        continue



                    if key == '+'  :
                        print(line )
                        print("+ intensionally ignored")
                        continue

                    if l[0] == 'SPC' :
                        ##https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2018/ENU/NSTRN-Reference/files/GUID-942678EC-4197-4D05-8CEC-6D61961444B0-htm.html
                        NLP.structured = True
                        SID = int(line[8:8*2])
                        tag = res.GetNodalTag("SPC"+str(SID))

                        G1 = int(line[8*2:8*3])
                        C1 = int(line[8*3:8*4])
                        D1 = float(line[8*4:8*5])
                        tag.AddToTag(filetointernalid[G1])
                        if len(NLP.GetField(line,5)) :
                            G2 = int(line[8*5:8*6])
                            C2 = int(line[8*6:8*7])
                            D2 = float(line[8*7:8*8])
                            tag.AddToTag(filetointernalid[G2 ])

                        continue
                    if key == 'PLOAD4' :
                        name = "PLOAD4_"+str(int(line[8:8*2]) )
                        tag = res.GetNodalTag(name)
                        idd = int(line[8*2:8*3])
                        node = filetointernalid[idd ]
                        ElementsAndNumber = filetointernalidElement[idd]
                        ElementsAndNumber[0].AddElementToTag(ElementsAndNumber[1],name)
                        continue

                    if key == 'FORCE' :
                        NLP.structured = True
                        key,name,oid,_,factor,fx,fy,fz = NLP.ParseTypes(line,[str,str,int,ignored,float,float,float,float])
                        name = "FORCE_"+name
                        node = filetointernalid[oid ]
                        tag = res.GetNodalTag(name).AddToTag(node)
                        val = [fx,fy,fz]
                        resdic[name] = {'N':node, "factor":factor, "val":val }
                        continue

                    if key == 'GRID' :
                        NLP.structured = True
                        key,oid,_,x,y,z = NLP.ParseTypes(line,[str,int,str,float,float,float])
                        ids.append(oid)
                        csystem= self.GetField(line,2)
                        if csystem != " "*8:
                            x,y,z = resdic["CordinateSystems"][int(csystem)].GetPointInOriginalSystem([x,y,z])
                        xs.append(x)
                        ys.append(y)
                        zs.append(z)
                        filetointernalid[oid] = cpt
                        cpt +=1
                        continue

                    if key ==  'CTRIA3':
                        oid = int(line[8:8*2])
                        idd2 = str(int(line[8*2:8*3]))
                        P1 = int(line[8*3:8*4])
                        P2 = int(line[8*4:8*5])
                        P3 = int(line[8*5:8*6])

                        conn = [filetointernalid[xx] for xx in  [P1, P2,P3 ] ]
                        elements = res.GetElementsOfType(ElementNames.Triangle_3)

                        localid = elements.AddNewElement(conn,oid)
                        filetointernalidElement[oid] = (elements,localid-1)

                        elements.tags.CreateTag("ET"+idd2,False).AddToTag(localid-1)
                        #res.AddElementToTagUsingOriginalId(oid,)
                        continue
                    if key ==  'CTETRA':
                        oid = int(line[8:8*2])
                        idd2 = str(int(line[8*2:8*3]))
                        points = []
                        for i in range(3,9):
                            try:
                                points.append(int(line[8*i:8*(i+1)]))
                            except:
                                pass

                        if len(points) > 4:
                            line = self.ReadCleanLine()
                            for i in range(1,5):
                                points.append(int(line[8*i:8*(i+1)]))

                        conn = [filetointernalid[xx] for xx in  points ]

                        if len(conn) == 4:
                            elements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
                        elif len(conn) == 10:
                            elements = res.GetElementsOfType(ElementNames.Tetrahedron_10)

                        localid = elements.AddNewElement(conn,oid)
                        filetointernalidElement[oid] = (elements,localid-1)

                        elements.tags.CreateTag("ET"+idd2,False).AddToTag(localid-1)
                        #res.AddElementToTagUsingOriginalId(oid,)
                        continue


                    if key == 'ENDDATA' : break
                    if line is None: break

                    print("key " + str(key) )# pragma: no cover
                    raise(ValueError("string  '" + str(line) + "' not treated"))# pragma: no cover

            if key == 'ENDDATA' : break

            print("key " + str(key) ) # pragma: no cover
            raise(ValueError("string '" + str(line) + "' not treated"))# pragma: no cover

        res.SetNodes(np.array([xs,ys,zs],dtype=float).T, ids )
        self.EndReading()
        res.PrepareForOutput()
        self.output = res
        self.data = resdic
        return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".fem",FemReader)


def CheckIntegrity(GUI = False):


    from BasicTools.Helpers.Tests import TestTempDir, WriteTempFile

    testData=u"""
DESOBJ(MIN) = 1
BEGIN BULK
DRESP1         1   wcomp   WCOMP
DRESP1         2 volfrac VOLFRAC
DTPL           1  PSOLID      72
        MEMBSIZ .001
DCONSTR        1       2        .1
DCONADD        2       1
DOPTPRM,OPTMETH,DUAL
DOPTPRM,OBJTOL,0.005
GAPPRM,GAPCMPL,      NO
PARAM,CHECKEL,      NO
PARAM,RENUMOK,     YES
PARAM,PRINFACC,       1
PARAM,CONTFEL,     YES
$$----------------------
$$  System Data
$$----------------------
$HMNAME SYSTCOL       80 "For_Support 1"
$HWCOLOR  SYSTCOL        80       9
$
CORD2R         1       0-.08441 -.0764260.      -.08441 -.076426-1.     +
+       .9155901-.0764260.
GRID      162969        53.00389234.456633.29904
GRID      156839        50.53256233.782536.31633
GRID      156554        50.67131231.871634.13537
GRID      146810        54.21149233.586335.82001
CTETRA    581279       2  162969  156839  156554  146810
CTRIA3    581280       3  162969  156839  156554
SPC            1  162969  1234560.0
FORCE          2  156839       01.0     2949.8  -5792.120.0
RBE2     1649872 162969  123456    156839    156554
CBUSH    581279      74  162969  156839                               1
PBUSH        74       K   RIGID   RIGID   RIGID
$$0000001111111122222222333333334444444455555555666666667777777788888888
RBE3     1649879          146810  1234561.0          123  162969  156554
+         156839
ENDDATA
+
"""
    filename = "CheckIntegrity_FemReader.fem"
    WriteTempFile(filename,testData)
    res = ReadFem(TestTempDir.GetTempPath() + filename)

    res = ReadFem(string=testData)

    FR = FemReader()
    FR.SetStringToRead(testData)
    res = FR.Read()

    print(res.nodes)
    print(res)
    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    WriteMeshToXdmf(TestTempDir().GetTempPath()+"FemReaderTest.xdmf",res)
    print(TestTempDir().GetTempPath())

    if GUI:# pragma: no cover
        import os
        os.system('vglrun paraview ' + TestTempDir().GetTempPath()+"FemReaderTest.xdmf")


    nlp=NastranLineParcer()

    assert nlp.ParseFloat("1") == 1
    assert nlp.ParseFloat("1.-3") == 1e-3
    assert nlp.ParseFloat("1.2-4") == 1.2e-4
    assert nlp.ParseFloat("1.23-5") == 1.23e-5
    assert nlp.ParseFloat("1.234-6") == 1.234e-6
    assert nlp.ParseFloat("1.2345-7") == 1.2345e-7

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(GUI = True))# pragma: no cover
