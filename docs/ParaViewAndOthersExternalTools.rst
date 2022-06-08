******************************
ParaViewAndOthersExternalTools
******************************

BasicTools offer some capabilities to communicate with external mesh dedicated tools.

Vtk [#vtk]_
###########

If you have a working VTK  installation, you can set the module :py:mod:`BasicTools.Containers.vtkBridge` to convert meshes back and forth BasicTools meshes.
::

    from BasicTools.Containers.vtkBridge import MeshToVtk, VtkToMesh
    BTMesh =  #<-- this is my BasicTools Mesh
    vtkMesh = MeshToVtk(BTMesh,TagsAsFields=True)
    # do some work with VTK
    vtk MeshII = #<- create a new vtk mesh from the output of a filter
    BTMeshII = VtkToMesh(vtkMeshII,FieldsAsTags=True)

ParaView [#paraview]_
#####################

Conda/Mamba Users Installation
******************************

For the moment this functionality is not available on conda environments installation even if the plugin is installed in the path ``/conda_env_path/ParaViewPlugins/BasicToolsParaViewBridge.py``
The reason is the incompatibility of the ParaView python with the conda python.
We are working on a solution for this problem

Developer Installation
**********************

Some functionalities (like readers, writers) can be added to ParaView as a plugin.
The configuration consists in setting the ``PYTHONPATH`` environment variable to your BasicTools installation

    ``PYTHONPATH=/path/to/BasicTools/src``

Then you can load the plugin ``/path/to/BasicTools/extras/BasicToolsParaViewBridge.py`` using the **Tools->Manage Plugins...** menu.
Also, you can set the ``PV_PLUGIN_PATH`` environment variable to indicate ParaView to load automatically the plugin at start up.

    ``PV_PLUGIN_PATH=/path/to/BasicTools/extras``

Three type of object are added to ParaView by the plugin:

* Readers: The BasicTools capabilities to reading data from different file formats.
* Writers: The BasicTools capabilities to export data to different file formats.
* Filters: Some of the mesh treatment functionalities of BasicTools are exposed as vtk filters.

Be aware that the use of this functionalities involve a format conversion between the vtk and the BasicTools internal format.
Be aware that your Python installation version may not be compatible with Python version of ParaView.


MeshIO [#meshio]_
###################
MeshIO is a library capable of reading and writing to various mesh file formats.

If you have a working MeshIO installation, you can set the module :py:mod:`BasicTools.Containers.MeshIOBridge` to convert meshes back and forth BasicTools meshes.
MeshIO offer some reading and writing capabilities.
More information in :py:mod:`BasicTools.Containers.MeshIOBridge`.

PyVista [#pyvista]_
###################
If you have a working PyVista installation, you can set the module :py:mod:`BasicTools.Containers.PyVistaBridge` to convert meshes back and forth BasicTools meshes.

PyVista offer a very simple interface for the visualisation of 3D meshes.
More information in BasicTools.Containers.PyVistaBridge


.. rubric:: Footnotes
.. [#vtk] https://vtk.org/
.. [#paraview] https://www.paraview.org/
.. [#meshio] https://github.com/nschloe/meshio
.. [#pyvista] https://www.pyvista.org/