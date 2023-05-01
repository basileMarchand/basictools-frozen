# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import os
import numpy as np
from BasicTools.NumpyDefs import PBasicFloatType

"""OpenFoam file reader.
Partial reader; only field are read at given time step, no mesh reader
"""


def ReadField(timeStep:int, fieldName:str, folderName:str = '.') -> np.ndarray:
    """Return a numpy.ndarray containing a field at a given time step
        read from an openfoam computation

    Parameters
    ----------
    timeStep : int
        time step at which the field is read
    fieldName : str
        name of the field to read
    folderName : str, optional
        name of the folder containing the openfoam computation., by default '.'

    Returns
    -------
    np.ndarray
        field being read
    """


    fileName = folderName+os.sep+str(timeStep)+os.sep+fieldName
    f = open(fileName,"r")
    allLines = f.readlines()
    f.close()

    data = []

    lines = iter(allLines)

    while True:
        line = next(lines)
        if "internalField" in line:
            line = next(lines)
            nbeOfCells = int(line.strip('\n()').split()[0])
            line = next(lines)
            break

    for i in range(nbeOfCells):
        line = next(lines)
        data.append([float(x) for x in line.strip('\n()').split()])

    line = next(lines)
    assert line.strip('\n').split()[0] == ")", "problem with end of data line, should be ')'"

    data = np.array(data, dtype = PBasicFloatType)

    return data

def CheckIntegrity():

    string = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2006                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0.0005";
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];


internalField   nonuniform List<vector>
9
(
(0.0728798 -0.581068 -0.00350291)
(-0.00403898 -1.21489 -0.948814)
(-0.107096 -0.742119 -0.255411)
(-0.0140631 0.474135 -0.830478)
(-0.0130636 -0.441406 -1.83613)
(-0.128549 -0.590474 -0.14046)
(0.113899 -1.04733 0.194137)
(-0.0719712 -0.886893 0.114526)
(-0.0117159 -1.4759 -0.410963)
)
;

boundaryField
{
    wall
    {
        type            slip;
    }
    inflow
    {
        type            fixedValue;
        value           uniform (0 -1 0);
    }
    outflow
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
"""

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()

    os.makedirs(tempdir+"0.1")

    f = open(tempdir+"0.1"+os.sep+"U", "w")
    f.write(string)
    f.close()

    field = ReadField(0.1, "U", tempdir)
    assert field[0,0] == 0.0728798
    assert field[-1,-1] == -0.410963

    return 'OK'

if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover
