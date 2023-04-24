*******************
Mesh File Converter
*******************

Located in the ``src/BasicTools/IO/`` folder, this application uses the readers and writers of BasicTools to convert mesh files to a new file format.
The tools uses the extension of the file name to select the correct reader/writer automatically.
If no output file is provided, a summary of the information present in the input file is printed.
Some Readers/Writer can handle fields (solution files) and/or time.

The command line option to show the help is :

.. code-block::

   MeshFileConverter -h

And the corresponding output

.. code-block::

   python  MeshFileConverter -i <inputfile> -o <outputfile>
   options:
        -i    Input file name
        -o    output file name
        -h    this help
        -t    time to read (if the input file can handle time)
             (default last time step is written)
        -p    print available times to read
        -v    more verbose output
        -m    Activate MeshIO Readers and Writers
        -s    Plot mesh before continue (press key "q" exit)

   Available Readers :
        ['.inp', '.asc', '.ansys', '.geof', '.geo',
        '.msh', '.mesh', '.meshb', '.sol', '.solb', '.gcode', '.fem',
        '.stl', '.xdmf', '.pxdmf', '.PIPE', '.odb', '.ut', '.utp', '.vtu',
        '.dat', '.datt']

   Available Writers :
        ['.PIPE', '.geof', '.msh', '.mesh', '.meshb',
        '.odb', '.stl', '.xdmf', '.xmf', '.csv']


The available Readers/Writers may vary dependent on the BasicTools version installed.
