# -*- coding: utf-8 -*-

#for python 2.6+ compatibility
from __future__ import print_function

import numpy as np



import BasicTools.FE.ElementNames as EN

OdbName = {}
OdbName[EN.Tetrahedron_4] = 'C3D4'
OdbName[EN.Triangle_3] = 'S3'
#OdbName[EN.Triangle_3] = 'S3R'
OdbName[EN.Bar_2] = 'CONN3D2'

def WriteMaterial(odb, material):
    import abaqusConstants as AC

    if material is None:
        #Creating material for odb
        pMat = odb.Material(name='Elastic Material')
        pMat.Elastic(type=AC.ISOTROPIC,
            temperatureDependency=AC.OFF, dependencies=0,
            noCompression=AC.OFF, noTension=AC.OFF,
            moduli=AC.LONG_TERM, table=((12000,0.3),))
    else:
        raise
    print("print Material")

def WriteSection(odb, section):

    ##Creating section for odb
    if  section is None:
         sectionName = 'Homogeneous Solid Section'
         mySection = odb.HomogeneousSolidSection( name = sectionName,
         material = 'Elastic Material',
         thickness = 1.0)
         pCat = odb.SectionCategory(name='odbSection',description = 'Section for odb')

    else:
        raise

def WriteOdb(filename,mesh,PointFields=None,CellFields=None,PointFieldsNames=None,CellFieldsNames=None,  __insubprocess= False, abaqusExec="abaqus"):


    if PointFields is None:
        PointFields = [];

    if CellFields  is None:
        CellFields   = [];

#    if GridFields is None:
#        GridFields  = [];

    if PointFieldsNames is None:
        PointFieldsNames  = [];

    if CellFieldsNames is None:
        CellFieldsNames  = [];

#    if GridFieldsNames is None:
#        GridFieldsNames  = [];


    try :
        import abaqusConstants as AC
        import odbAccess
    except :
        if __insubprocess :
            print("Error Loading libraries in the subprocess")
            return
        # it was not possible to load the libraries, we tried to launch a
        # writing service
        import BasicTools.IO.CodeInterface as CodeInterface
        import os
        from BasicTools.IO.Wormhole import WormholeClient


        path = os.sep.join(filename.split("/")[0:-1])

        interface = CodeInterface.Interface(path)

        from BasicTools.Helpers.Tests import TestTempDir

        interface.SetWorkingDirectory(TestTempDir.GetTempPath())
        interface.processDirectory = TestTempDir.GetTempPath()
        absfilename   = os.path.abspath(filename)

        port = 12346
        interface.tpl = """
from BasicTools.IO.Wormhole import WormholeServer
WormholeServer(""" +str(port) +""",dry=False)
"""
        interface.inputFilename = "ServerCode"
        interface.inputFileExtension = ".py"

        interface.WriteFile(0)


        #interface.SetCodeCommand('"C:\\Program Files (x86)\\Notepad++\\notepad++.exe"')
        from BasicTools.Helpers.which import which
        if which(abaqusExec) is None:
            raise Exception("Abaqus not available in your system")
        interface.SetCodeCommand(abaqusExec + " python")
        import sys
        proc = interface.SingleRunComputation(0)
        import time
        time.sleep(2)
        print(interface.lastCommandExecuted)

        client = WormholeClient()
        client.Connect(port)
        client.SendData("filename",absfilename)
        print("Writing file :")
        print(absfilename)
        client.SendData("mesh",mesh)
        client.SendData("PointFields",PointFields)
        client.SendData("CellFields",CellFields)
        client.SendData("PointFieldsNames",PointFieldsNames)
        client.SendData("CellFieldsNames",CellFieldsNames)
        client.RemoteExec("from BasicTools.IO.OdbWriter import WriteOdb")
        client.RemoteExec("WriteOdb(filename,mesh,PointFields,CellFields,PointFieldsNames,CellFieldsNames)")
        client.RemoteExec("print('Done')")
        client.Exit()

        proc.wait()
        return

    mesh.PrepareForOutput()

    odbName = filename.split("/")[-1]
    pOdb = odbAccess.Odb(name=odbName,analysisTitle='MyFirstAnalisys',path=filename,description='1D beam')
    pOdb.save()

    WriteMaterial(pOdb,None)
    WriteSection(pOdb,None)


    ##Creating the 3D solid part
    pPart = pOdb.Part(name='beamTaylor',embeddedSpace = AC.THREE_D,type= AC.DEFORMABLE_BODY)
    nodeLabels = list(range(1,1+mesh.GetNumberOfNodes()))
    pPart.addNodes(labels = nodeLabels,coordinates = mesh.GetPosOfNodes())


    ##Create a node and element set
    for tag in mesh.nodesTags:
        #print(pPart.nodes)
        #print(pPart.nodes[1])
        #print(tag.GetIds().astype(np.int32))
        #print(pPart.nodes[tag.GetIds().astype(np.int32)])
        #print("----")
        pPartNodes= pPart.NodeSetFromNodeLabels(tag.name,(tag.GetIds()+1).astype(np.int32))
#

    cpt =0;
    for ntype, data in mesh.elements.items():
        elemtype = OdbName[ntype]

        dd = list(range(data.globaloffset+1,data.globaloffset+1+data.GetNumberOfElements()))
        dd = np.array(dd,dtype=np.int32)
        print(dd)
        print(data.connectivity.astype(np.int32))
        print(elemtype)
        print(ntype)

        pPart.addElements(labels = dd , connectivity=(data.connectivity+1).astype(np.int32) , type = elemtype, elementSetName=ntype)

#

    for name in mesh.GetNamesOfElemTags():
        ids = mesh.GetElementsInTag(name)+1
        elementSet = pPart.ElementSetFromElementLabels(    name=name,elementLabels=ids.astype(np.int32))


##Creating the instance for the solid part
    pAssembly = pOdb.rootAssembly.Instance(name = 'Principal',object = pPart)
    #pOdb.update()
    #pOdb.save()
    #pOdb.close()
    #return

#
##Creating section for odb
#
##Creating the analysis step
    pStep = pOdb.Step(name='StaticAnalysis',description='Analysis type - 101',domain=AC.TIME,timePeriod = 1.0)
    pFrame0 = pStep.Frame(incrementNumber=0,frameValue=0.0000)
    #pDisp0 = pFrame0.FieldOutput(name='U',description='Displace ment',type=AC.VECTOR,componentLabels=('1','2','3'))

    #dispValues = np.hstack((np.array(nodeLabels)[:,np.newaxis],)*3).ravel().astype(float)
    #dispValues.shape = (len(nodeLabels),3)
    #dispValues /= len(nodeLabels)*10

    #pDisp0.addData(position = AC.NODAL,instance = pAssembly,labels = nodeLabels,data=dispValues)



    cellLabels = range(1,1+mesh.GetNumberOfElements())
    for i in range(len(CellFieldsNames)):
        if mesh.GetNumberOfElements() == CellFields[i].size:
            ftype = AC.SCALAR
        elif mesh.GetNumberOfElements()*3 == CellFields[i].size:
            ftype = AC.VECTOR
        else:
            raise
        fo = pFrame0.FieldOutput(name=CellFieldsNames[i],description=CellFieldsNames[i],type=ftype)
        fo.addData(position = AC.CENTROID,instance = pAssembly,labels = cellLabels,data=CellFields[i])

    for i in range(len(PointFieldsNames)):
        if mesh.GetNumberOfNodes() == PointFields[i].size:
            ftype = AC.SCALAR
        elif mesh.GetNumberOfNodes()*3 == PointFields[i].size:
            ftype = AC.VECTOR
        else:
            raise
        fo = pFrame0.FieldOutput(name=PointFieldsNames[i],description=PointFieldsNames[i],type=ftype )
        fo.addData(position = AC.NODAL,instance = pAssembly,labels = nodeLabels,data=PointFields[i])


##Creating the frame for the step
#
#pFrame1 = pStep.Frame(incrementNumber=1,frameValue=1.0000)
#
##Reading the result file
#dispValues = ny.loadtxt('DISPL_POINTS.dat',skiprows=1,usecols = (3,4,5))
#
##Creating the Field Output - Displacement
#
#pDisp1 = pFrame1.FieldOutput(name='U',description='Displace ment',type=VECTOR,componentLabels=('1','2','3'))
#
##Adding data
#
#pDisp1.addData(position = NODAL,instance = pAssembly,labels = nodeLabels,data=dispValues)
#
##Setting default display options
#pStep.setDefaultField(pDisp0)
#
    pOdb.update()
    pOdb.save()
    pOdb.close()

def CheckIntegrity():
    import BasicTools.FE.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()
    print(tempdir)

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,1]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=np.int)

    tag = mymesh.nodesTags.CreateTag("coucou")
    tag.AddToTag(0)
    tag.AddToTag(1)

    tet = mymesh.GetElementsOfType(EN.Tetrahedron_4)
    tet.AddNewElement([0,1,2,3],0)
    tet.tags.CreateTag("TheOnlyTet").AddToTag(0)

    tris = mymesh.GetElementsOfType(EN.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)
    tris.tags.CreateTag("OneTri").AddToTag(0)

    #bars = mymesh.GetElementsOfType(EN.Bar_2)
    #bars.AddNewElement([0,1],0)
    #bars.AddNewElement([1,3],1)
    #bars.tags.CreateTag("firstBar").AddToTag(0)

    #point = mymesh.GetElementsOfType(EN.Point_1)
    #point.AddNewElement([0],0)
    #point.tags.CreateTag("onlyPoint").AddToTag(0)


    #mymesh.AddElementToTagUsingOriginalId(3,"Tag1")
    #mymesh.AddElementToTagUsingOriginalId(5,"Tag3")

    tempdir = "./"
    WriteOdb(tempdir+"Test_OdbWriter.odb",mymesh)

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
