from BasicTools.Containers import UnstructuredMesh
from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
from BasicTools.Containers import UnstructuredMeshCreationTools as UMCT
from scipy.spatial import Delaunay
import numpy as np

def create_mesh_example() -> UnstructuredMesh:
    """This function creates a dummy mesh using Delaunay triangulation.
    """
    # Create a grid of points
    x, y = np.meshgrid(np.linspace(0,1,10), np.linspace(0,1,10))
    points = np.stack([x.ravel(), y.ravel()], axis=1)
    # Generate the triangles with Delaunay
    tri = Delaunay(points)
    triangles = tri.simplices
    # Create a BasicTools UnstructuredMesh using the CreateMeshOfTriangles utility
    mesh = UMCT.CreateMeshOfTriangles(points, triangles)
    return mesh

if __name__ == "__main__":

    # Create a simple mesh
    mesh = create_mesh_example()
    
    # Make six dummy nodal fields
    fields = np.random.randn(mesh.GetNumberOfNodes(), 6)
    
    # Dump the mesh and nodal fields into a XDMF file
    WriteMeshToXdmf(filename="WriteXdmf.xdmf",                          # path where the file will be stored
                    baseMeshObject=mesh,                                # UnstructuredMesh object 
                    PointFields=[fields[:,i] for i in range(6)],        # list of scalar fields
                    PointFieldsNames=[f"field_{i}" for i in range(6)],  # list of names for each scalar field
    )