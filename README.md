gmx2qmmm
======
`gmx2qmmm` is a python interface for QM/MM computation.

## Overview

gmx2qmmm is a python package to bridge Gaussian and Gromacs. The test runs were performed using Gaussian16 and Gromacs 5.0.2, but the code should be able to read earlier Gaussian and other Gromacs versions. The only limits are the formats of the human-readable input and output files of each program, as such, conversion scripts can be written to make the interface work with any version, if the current code does not support it.
Conceptually, gmx2qmmm creates a QM/MM potential and performs either single point calculations (i.e., the current energy of your system) or optimizations (via steepest decent). Other functions are easily implemented in the future.

## System requirments
 - [python 2.7] (numpy, sqlite3)
 - [Gaussian16]
 - [Gromacs 5.0.2]
 
## Input files
The required input files in the package are

|Input files|Command|Default input name|
| ------ | ------ | ------ |
|Coordinate file (.g96)|-c|conf.g96|
|Topology (.top)|-p|topol.top|
|QM atoms file(.ndx)|-n|qmatoms.ndx|
|QM parameters (.dat)|-qm|qm.dat|
|MM parameters (.dat)|-mm|mm.dat|
|QM/MM parameters (.dat)|-qmmm|qmmm.dat|
|Active atoms (.ndx)|-act|act.ndx|
|Path file (.dat)|-path|path.dat|
|Logfile (.log)|-g|logfile|

If the parameters are not entered in the input files, the program will run the job with the default values.

Input files guide please enter

**python gmx2qmmm.py -h**

The executed example with default parameters is

**python gmx2qmmm.py -c conf.gro -p topol.top -n qmatoms -qm qm -mm mm -qmmm qmmm -act active_atoms -path path -g logfile**

For advance information please read [gmx2qmmm reference].
## Support and development
For bug reports/suggestions/complaints please raise an issue on [GitHub].

Or contact us directly: [gmx2qmmm@gmail.com]


[python 2.7]:<https://www.python.org/download/releases/2.7>
[Gaussian16]:<https://gaussian.com/gaussian16/>
[Gromacs 5.0.2]:<http://www.gromacs.org>
[GitHub]:<https://github.com/gmx2qmmm/gmx2qmmm_portable>
[gmx2qmmm@gmail.com]:<mailto:gmx2qmmm@gmail.com>
[gmx2qmmm reference]:<https://drive.google.com/file/d/1ynMEmijAfRfvyLQYSU30lsO8aqFmTrgx/view?usp=sharing>
