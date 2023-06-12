

cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM
import BasicTools.Containers.NativeUnstructuredMesh as NUM

cimport BasicTools.FE.Numberings.NativeDofNumbering as cNDN
import BasicTools.FE.Numberings.NativeDofNumbering as NDN

cimport BasicTools.CythonDefs as CD
import BasicTools.NumpyDefs as ND
cimport BasicTools.Containers.NativeTransfer as cNT
import BasicTools.Containers.ElementNames as ElementNames
import numpy as np
cimport numpy as np
from libcpp cimport bool


cdef class NativeTransfer:
#    def __init__(self):
#       self.mesh = None

    cdef cNT.TransferClass* GetCppPointer(self) nogil:
        return &(self.cpp_object)

    def ToStr(self):
        return str(self.cpp_object.ToStr());

    def SetVerbose(self, verbose:bool):
        self.cpp_object.SetVerbose(verbose)

    def SetTransferMethod(self, method):
        self.cpp_object.SetTransferMethod(method.encode())

    def GetTransferMethod(self):
        return str(self.cpp_object.GetTransferMethod())

    def SetSourceFEField(self, sourceField, elementFilter = None ):
        from BasicTools.Containers.Filters import ElementFilter
        from BasicTools.Containers.UnstructuredMeshModificationTools import CleanLonelyNodes
        from BasicTools.Containers.UnstructuredMeshInspectionTools import ExtractElementsByElementFilter

        self.sourceMesh = NUM.CUnstructuredMesh()
        iMeshDim = sourceField.mesh.GetElementsDimensionality()

        if elementFilter == None:
            elementFilter = ElementFilter(mesh = sourceField.mesh, dimensionality= iMeshDim)


        iMesh = ExtractElementsByElementFilter(sourceField.mesh, elementFilter, copy=False )
        CleanLonelyNodes(iMesh)

        if iMesh.nodes.shape[1] < 3 :
            newNodes = np.zeros((iMesh.nodes.shape[0],3), dtype= ND.PBasicFloatType, order="c")
            newNodes[0:,:iMesh.nodes.shape[1]] = iMesh.nodes
            iMesh.nodes = newNodes

        self.pythonmesh = iMesh

        self.sourceMesh.SetDataFromPython(iMesh)
        cdef cNUM.CUnstructuredMesh sourceMesh = self.sourceMesh
        self.cpp_object.SetSourceMesh(sourceMesh.GetCppPointer())

        # hack to find the correct name of the space used
        from BasicTools.FE.Spaces.FESpaces import AllSpaces
        for spaceName, obj in AllSpaces.items():
            if sourceField.space is obj:
                self.cpp_object.SetSourceSpace(spaceName.encode())
                break
        else:
            raise Exception("Space not available on the cpp side.")

        cdef cNDN.NativeDofNumbering ndn = sourceField.numbering
        cdef cNDN.DofNumbering* dn = ndn.GetCppPointer()
        self.cpp_object.SetSourceNumbering(dn)

    def SetTargetPoints(self, cnp.ndarray[CBasicFloatType, ndim=2,mode="c"] targetPoints not None ):

        if targetPoints.shape[1] == 3 :
            self.cpp_object.SetTargetPoints(FlattenedMapWithOrder[Matrix, CBasicFloatType, Dynamic, Dynamic, RowMajor](targetPoints))
            return

        cdef cnp.ndarray[CBasicFloatType, ndim=2,mode="c"]  newTargetPoints = np.zeros((targetPoints.shape[0],3), dtype= ND.PBasicFloatType, order="c")
        newTargetPoints[0:,:targetPoints.shape[1]] = targetPoints
        self.pythonNewTargetPoints = newTargetPoints
        self.cpp_object.SetTargetPoints(FlattenedMapWithOrder[Matrix, CBasicFloatType, Dynamic, Dynamic, RowMajor](newTargetPoints))

    def Compute(self):
        #print("Start Compute cython ")
        self.cpp_object.Compute()
        #print("Done Compute cython")

    def GetOperator(self):
        #print("Start GetOperator")
        from scipy.sparse import coo_matrix
        cdef CBasicIndexType nb_source_Dofs
        nb_source_Dofs  = self.cpp_object.nb_source_Dofs
        cdef CBasicIndexType targetPoints
        targetPoints =   self.cpp_object.nb_targetPoints

        if self.cpp_object.data.size() == 0:
            #print( (targetPoints,nb_source_Dofs) )
            res = coo_matrix(([],([],[])), shape= (targetPoints,nb_source_Dofs), copy=True)
            #print("Done GetOperator empty")
            return res

        cdef CBasicFloatType[::1] d = <CBasicFloatType[:self.cpp_object.data.size()]>self.cpp_object.data.data()
        cdef CBasicIndexType[::1] r = <CBasicIndexType[:self.cpp_object.rows.size()]>self.cpp_object.rows.data()
        cdef CBasicIndexType[::1] c = <CBasicIndexType[:self.cpp_object.cols.size()]>self.cpp_object.cols.data()

        #print( (targetPoints,nb_source_Dofs) )
        res = coo_matrix((d,(r,c)), shape= (targetPoints,nb_source_Dofs), copy=True)
        #print("Done GetOperator")
        return res

    def GetStatus(self):
        return ndarray_copy(self.cpp_object.GetStatus())

    def SetUseEdges(self, val:bool):
            self.cpp_object.useEdges = val;
