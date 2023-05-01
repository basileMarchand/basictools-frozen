# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Ut file writer (Zset mesh files)
"""

import numpy as np
import os

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.IO.GeofWriter as GW
from BasicTools.IO.GeofWriter import GeofName as GeofName
from BasicTools.IO.GeofReader import nbIntegrationsPoints as nbIntegrationsPoints
from BasicTools.Containers.Filters import ElementFilter


class UtWriter(WriterBase):
    """Class to write an ut file, and the files of the data defined in this ut file.
    This class can generate a .ut, .geof, .ctnod, .node, .integ to vizualise finite element computational results
    """
    def __init__(self):
        super(UtWriter,self).__init__()
        self.canHandleTemporal = True
        self.canHandleAppend = True

    def __str__(self):
        res  = 'UtWriter : \n'
        res += '   Name : '+str(self.fileName)+'\n'
        return res

    def Open(self):
        pass

    def SetFileName(self, fileName):
        """Sets the name of ut file to write

        Parameters
        ----------
        fileName : str
            file name to set
        """
        if fileName is None:
            self.fileName = None
            self.folder = None
            return
        from os import path
        self.fileName = path.basename(fileName);
        self.folder = path.dirname(fileName)

    def GetFolder(self):
        """Returns folder where the file is written

        Returns
        -------
        str
            folder where the file is written
        """
        return self.folder + os.sep

    def GetBaseName(self):
        """Returns the file name to write without extension

        Returns
        -------
        str
            file name to write without extension
        """
        from os import path
        return path.splitext(self.fileName)[0]

    def AttachMesh(self, mesh):
        """Attach mesh to writer

        Parameters
        ----------
        mesh : UnstructuredMesh
            mesh to attach to writer
        """
        self.mesh = mesh

    def AttachData(self, data_node, data_ctnod = {}, data_integ = {}, Nnode = None, Nint = None):
        """Attach fields to the writer

        Parameters
        ----------
        data_node : dict
            primal fields, defined at the vertices of the mesh
        data_ctnod : dict, optional
            dual fields, defined at the vertices of the mesh, by default {}
        data_integ : dict, optional
            dual fields, defined at the integration points of the mesh, by default {}
        Nnode : int, optional
            number of nodes, by default None
        Nint : int, optional
            number of integration points, by default None
        """
        self.data_node        = data_node
        self.data_ctnod       = data_ctnod
        self.data_integ       = data_integ
        self.data_node_names  = list(data_node.keys())
        self.data_integ_names = list(data_integ.keys())
        self.NnodeVar         = len(data_node)
        self.NintVar          = len(data_integ)
        if Nnode is not None:
            self.Nnode = Nnode
        else:
            if len(self.data_node_names):
                self.Nnode            = data_node[self.data_node_names[0]].shape[0]
                #print("warning empty self.dat_node_names")

        if Nint is not None:
            self.Nint = Nint
        else:
            try:
                self.Nint = data_integ[self.data_integ_names[0]].shape[0]
            except IndexError:
                True

    def AttachDataFromProblemData(self, problemData, tag, Nnode = None, Nint = None):
        """Attach the fields defined in a problemData to the writer

        Parameters
        ----------
        problemData : ProblemData
            description of a physical problem
        tag : str
            solution name to attache to the writer
        Nnode : int, optional
            number of nodes, by default None
        Nint : int, optional
            number of integration points, by default None
        """
        self.AttachData(problemData.solutions[tag].data_node, problemData.solutions[tag].data_ctnod, problemData.solutions[tag].data_integ, Nnode = Nnode, Nint = Nint)

    def AttachSequence(self, timeSequence):
        """Attach the time sequence to the writer

        Parameters
        ----------
        timeSequence : np.ndarray
            time sequence to attach to the writer
        """
        self.timeSequence = timeSequence
        self.Ntime = timeSequence.shape[0]


    def WriteMesh(self):
        """Write mesh to file in Geof format
        """
        if self.mesh == None:
            print("please attach a mesh to the UtWriter object to be able to write a mesh; script terminated")
            exit()
        else:
            OW = GW.GeofWriter()
            OW.Open(self.GetFolder()+self.fileName+".geof")
            OW.Write(self.mesh, useOriginalId=True, lowerDimElementsAsElsets=False)
            OW.Close()


    def InitWrite(self, writeGeof, geofName = None, skipCtnod = False):
        """Initialize the written files

        Parameters
        ----------
        writeGeof : bool
            if True, the mesh is written in Geof format
        geofName : str, optional
            name of the mesh file, by default None
        skipCtnod : bool, optional
            if True, the dual quantities will not be exported at the verticies, by default False

        Returns
        -------
        _io.TextIOWrapper, _io.TextIOWrapper, _io.TextIOWrapper
            nodeFile, ctnodFile and integFile initialized
        """
        if geofName == None:
            __string = u"**meshfile "+self.fileName+".geof\n"
        else:
            __string = u"**meshfile "+geofName+"\n"

        if writeGeof==True:
            self.WriteMesh()

        if self.NnodeVar>0:
            __string += "**node "
            for field in self.data_node_names:
                __string += field+" "
            __string += "\n"

        if self.NintVar>0:
            __string += "**integ "
            for field in self.data_integ_names:
                __string += field+" "
            __string += "\n"

        __string += "**element\n"

        with open(self.GetFolder()+self.GetBaseName()+".ut", "w") as f:
            f.write(__string)
        f.close()

        for ff in [".node",".ctnod",".integ"]:
            try:
                os.remove(self.GetFolder()+self.GetBaseName()+ff)
            except OSError:
                pass

        nodeFile = open(self.GetFolder()+self.GetBaseName()+".node", "a")
        if skipCtnod == False:
            ctnodFile = open(self.GetFolder()+self.GetBaseName()+".ctnod", "a")
        else:
            ctnodFile = None
        if self.NintVar > 0:
            integFile = open(self.GetFolder()+self.GetBaseName()+".integ", "a")
        else:
            integFile = None

        return nodeFile, ctnodFile, integFile


    def WriteTimeStep(self, nodeFile, ctnodFile, integFile, timeSequenceStep, skipCtnod = False):
        """Write one time step into initialized .ut, .ctnod, .node and .integ files

        Parameters
        ----------
        nodeFile : _io.TextIOWrapper
            .node file containing primal fields at verticies
        ctnodFile : _io.TextIOWrapper
            .ctnod file containing dual fields at verticies
        integFile : _io.TextIOWrapper
            .ctnod file containing dual fields at integration points
        timeSequenceStep : np.ndarray
            time step information
        skipCtnod : bool, optional
            if True, the dual quantities will not be exported at the verticies, by default False
        """
        if self.NnodeVar>0:
            data_node = np.empty(self.NnodeVar*self.Nnode)
            for k in range(self.NnodeVar):
                data_node[k*self.Nnode:(k+1)*self.Nnode] = self.data_node[self.data_node_names[k]]
            data_node.astype(np.float32).byteswap().tofile(nodeFile)
            del data_node

        if self.NintVar>0:
            if skipCtnod == False:
                data_ctnod = np.empty(self.NintVar*self.Nnode)
                for k in range(self.NintVar):
                    data_ctnod[k*self.Nnode:(k+1)*self.Nnode] = self.data_ctnod[self.data_integ_names[k]]
                data_ctnod.astype(np.float32).byteswap().tofile(ctnodFile)
                del data_ctnod

            numberElements = []
            nbPtIntPerElement = []

            for name, data, ids in ElementFilter(self.mesh, dimensionality=3):
                numberElements.append(data.GetNumberOfElements())
                nbPtIntPerElement.append(nbIntegrationsPoints[GeofName[name]])
            nbTypeEl = len(numberElements)

            data_integ = np.empty(self.NintVar*self.Nint)
            count0 = 0
            field = np.empty((self.NintVar,self.Nint))
            for k in range(self.NintVar):
                field[k,:] = self.data_integ[self.data_integ_names[k]]
            for l in range(nbTypeEl):
                for m in range(numberElements[l]):
                    for k in range(self.NintVar):
                        data_integ[count0:count0+nbPtIntPerElement[l]] = field[k,nbPtIntPerElement[l]*m:nbPtIntPerElement[l]*m+nbPtIntPerElement[l]]
                        count0 += nbPtIntPerElement[l]
            data_integ.astype(np.float32).byteswap().tofile(integFile)
            del field; del data_integ

        __string = str(int(timeSequenceStep[0]))+" "+str(int(timeSequenceStep[1]))+" "+str(int(timeSequenceStep[2]))+" "+str(int(timeSequenceStep[3]))+" "+str(timeSequenceStep[4])+"\n"

        with open(self.GetFolder()+self.GetBaseName()+".ut", "a") as f:
            f.write(__string)
        f.close()

    def Write(self, outmesh, PointFieldsNames, PointFields, CellFieldsNames = None, CellFields = None):
        """Initialize writer and write .ut, .geof, .ctnod, .node and .integ files

        Parameters
        ----------
        outmesh : UnstructuredMesh
            the mesh to be written
        PointFieldsNames : list[str]
            name of the fields to write defined at the vertices of the mesh
        PointFields : list[np.ndarray]
            fields to write defined at the vertices of the mesh
        CellFieldsNames : None
            Not Used, by default None
        CellFields : None
            Not Used, by default None
        """

        self.AttachMesh(outmesh)
        self.AttachData( {x:y for x,y in zip(PointFieldsNames,PointFields) }, {}, {})
        self.AttachSequence(np.array([[0,0,0,0,0]]))
        self.WriteFiles(writeGeof=True)


    def WriteFiles(self, writeGeof, geofName = None, skipCtnod = False):
        """Write .ut, .geof, .ctnod, .node and .integ files

        Parameters
        ----------
        writeGeof : bool
            if True, the mesh is written in Geof format
        geofName : str, optional
            name of the mesh file, by default None
        skipCtnod : bool, optional
            if True, the dual quantities will not be exported at the verticies, by default False
        """
        nodeFile, ctnodFile, integFile = self.InitWrite(writeGeof, geofName, skipCtnod)

        if self.NnodeVar > 0:
            data_node = np.empty(self.NnodeVar*self.Nnode*self.Ntime)
            for i in range(self.Ntime):
                for k in range(self.NnodeVar):
                    data_node[self.NnodeVar*self.Nnode*i+k*self.Nnode:self.NnodeVar*self.Nnode*i+(k+1)*self.Nnode] = self.data_node[self.data_node_names[k]][:,i]
            data_node.astype(np.float32).byteswap().tofile(nodeFile)
            del data_node

        if self.NintVar > 0:
            if skipCtnod == False:
                data_ctnod = np.empty(self.NintVar*self.Nnode*self.Ntime)
                for i in range(self.Ntime):
                    for k in range(self.NintVar):
                        data_ctnod[self.NintVar*self.Nnode*i+k*self.Nnode:self.NintVar*self.Nnode*i+(k+1)*self.Nnode] = self.data_ctnod[self.data_integ_names[k]][:,i]
                data_ctnod.astype(np.float32).byteswap().tofile(ctnodFile)
                del data_ctnod

            numberElements = []
            nbPtIntPerElement = []
            if self.mesh.GetNumberOfElements(3) != 0:
                dim = 3
            else:
                dim = 2
            for name, data, ids in ElementFilter(self.mesh, dimensionality=3):
                numberElements.append(data.GetNumberOfElements())
                nbPtIntPerElement.append(nbIntegrationsPoints[GeofName[name]])
            nbTypeEl = len(numberElements)

            data_integ = np.empty(self.NintVar*self.Nint*self.Ntime)
            count0 = 0
            for i in range(self.Ntime):
                field = np.empty((self.NintVar,self.Nint))
                for k in range(self.NintVar):
                    field[k,:] = self.data_integ[self.data_integ_names[k]][:,i]
                for l in range(nbTypeEl):
                    for m in range(numberElements[l]):
                        for k in range(self.NintVar):
                            data_integ[count0:count0+nbPtIntPerElement[l]] = field[k,nbPtIntPerElement[l]*m:nbPtIntPerElement[l]*m+nbPtIntPerElement[l]]
                            count0 += nbPtIntPerElement[l]
            data_integ.astype(np.float32).byteswap().tofile(integFile)
            del field; del data_integ

        __string = ""
        for i in range(self.Ntime):
            __string += str(int(self.timeSequence[i,0]))+" "+str(int(self.timeSequence[i,1]))+" "+str(int(self.timeSequence[i,2]))+" "+str(int(self.timeSequence[i,3]))+" "+str(self.timeSequence[i,4])+"\n"

        with open(self.GetFolder()+self.GetBaseName()+".ut", "a") as f:
            f.write(__string)
        f.close()

        nodeFile.close()
        if skipCtnod == False:
            ctnodFile.close()
        if self.NintVar > 0:
            integFile.close()

    def Close(self):
        # to override the WriterBase().Close()
        pass

def CheckIntegrity():

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    timeSequence = reader.time

    reader.atIntegrationPoints = False
    Nnode = reader.ReadField(fieldname="U1", timeIndex=0).shape[0]
    reader.atIntegrationPoints = True
    Nint = reader.ReadField(fieldname="sig11", timeIndex=0).shape[0]

    Ntime = reader.time.shape[0]
    NnodeVar = len(reader.node)
    NintVar = len(reader.integ)

    import collections
    data_node = collections.OrderedDict()
    for dnn in reader.node:
        data_node[dnn] = np.empty((Nnode,Ntime))
    data_ctnod = collections.OrderedDict()
    data_integ = collections.OrderedDict()
    for din in reader.integ:
        data_ctnod[din] = np.empty((Nnode,Ntime))
        data_integ[din] = np.empty((Nint,Ntime))

    reader.atIntegrationPoints = False
    for i in range(Ntime):
        for dn in data_node:
            data_node[dn][:,i] = reader.ReadField(fieldname=dn, timeIndex=i)
        for dc in data_ctnod:
            data_ctnod[dc][:,i] = reader.ReadField(fieldname=dc, timeIndex=i)
    reader.atIntegrationPoints = True
    for i in range(Ntime):
        for di in data_integ:
            data_integ[di][:,i] = reader.ReadField(fieldname=di, timeIndex=i)

    import BasicTools.IO.GeofReader as GR
    mymesh = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")

    ##################################
    # EXEMPLE SYNTAXE DU WRITER
    import BasicTools.IO.UtWriter as UW
    UtW = UW.UtWriter()
    UtW.SetFileName(tempdir+os.sep+"toto.ut")
    UtW.AttachMesh(mymesh)
    UtW.AttachData(data_node, data_ctnod, data_integ)
    UtW.AttachSequence(timeSequence)
    UtW.WriteFiles(writeGeof=True)
    ##################################

    print(UtW)

    print("Temp directory =", tempdir)

    import filecmp
    res = filecmp.cmp(tempdir + "toto.node",  BasicToolsTestData.GetTestDataPath() + "UtExample/cube.node", shallow=False)
    print("node files equals  ?", res )
    if res == False:
        return 'Problem comparing toto.node'
    res = filecmp.cmp(tempdir + "toto.ctnod", BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ctnod", shallow=False)
    print("ctnod files equals ?", res)
    if res == False:
        return 'Problem comparing toto.ctnode'

    res = filecmp.cmp(tempdir + "toto.integ", BasicToolsTestData.GetTestDataPath() + "UtExample/cube.integ", shallow=False)
    print("integ files equals ?", res)
    if res == False:
        return 'Problem comparing toto.integ'
    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
