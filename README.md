gmx2qmmm v.1.0.2
======
**gmx2qmmm** is a python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculation.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents
- [Overview](#Overview)
- [System requirments](#system-requirments)
- [**gmx2qmmm** job](#job-type)
- [Input files](#input-files)
- [Examples](#example)
- [References](#references)
- [Support and development](#support-and-development)
<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Overview

**gmx2qmmm** is a python package to bridge [Gaussian] and [Gromacs]. The test runs were performed using [Gaussian16] and [Gromacs 5.0.2], but the code should be able to read earlier Gaussian and other Gromacs versions. The only limits are the formats of the human-readable input and output files of each program, as such, conversion scripts can be written to make the interface work with any version, if the current code does not support it.
Conceptually, **gmx2qmmm** creates a QM/MM potential and performs either single point calculations (i.e., the current energy of your system) or geometry optimizations. Other functions are easily implemented in the future.

## System requirments
 - [python 2.7] (numpy, sqlite3)
 - [Gaussian16] (and earlier version)
 - [Gromacs 5.0.2] (and earlier version)
## **gmx2qmmm** job
|Job type|Default input name|
| ------ | ------ |
|Single point calcuation (SP)|Calculate single point energy and forces (.xyz) |
|Geometry optimizations (OPT)|Optimize the system energy via optimizer ([Steepest descent], [Conjugate gradient] or [BFGS])|

[Steepest descent]:<https://en.wikipedia.org/wiki/Gradient_descent>
[Conjugate gradient]:<https://en.wikipedia.org/wiki/Conjugate_gradient_method>
[BFGS]:<https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm>
## Input files
The required input files in the package are

|Input files|Command|Default input name|
| ------ | ------ | ------ |
|Coordinate file (.g96 or .gro)|-c|conf.g96|
|Topology (.top)|-p|topol.top|
|QM atoms file(.ndx)|-n|qmatoms.ndx|
|QM parameters (.dat)|-qm|qm.dat|
|MM parameters (.dat)|-mm|mm.dat|
|QM/MM parameters (.dat)|-qmmm|qmmm.dat|
|Active atoms (.ndx)|-act|act.ndx|
|Path file (.dat)|-path|path.dat|
|Logfile (.log)|-g|logfile|

If the parameters are not entered in the input files, the program will run the job with the default values.

Input files guide please execute  `python gmx2qmmm.py -h`

The executed example with default parameters is  `python gmx2qmmm.py`

For advance information please read [gmx2qmmm reference].

## Examples
The directory example contains SP and OPT calculation for glycine serine (GLYSER). 
![alt text](https://github.com/gmx2qmmm/gmx2qmmm_portable/blob/master/example/glyser.png?raw=true)

The names of the input files are default. Go to the sp/opt directory on the command line and run:

```
python gmx2qmmm.py
```

## References

> A user‐friendly, Python‐based quantum mechanics/Gromacs interface: gmx2qmmm
> Jan P. Götze, Yuan‐Wei Pi, Simon Petry, Fabian Langkabel,  Jan Felix Witte, Oliver Lemke
> https://doi.org/10.1002/qua.26486

## Support and development
For bug reports/suggestions/complaints please raise an issue on [GitHub].

Or contact us directly: [gmx2qmmm@gmail.com]


[python 2.7]:<https://www.python.org/download/releases/2.7>
[Gaussian16]:<https://gaussian.com/gaussian16/>
[Gromacs 5.0.2]:<http://www.gromacs.org>
[Gaussian]:<https://gaussian.com/gaussian16/>
[Gromacs]:<http://www.gromacs.org>
[GitHub]:<https://github.com/gmx2qmmm/gmx2qmmm_portable>
[gmx2qmmm@gmail.com]:<mailto:gmx2qmmm@gmail.com>
[gmx2qmmm reference]:<https://drive.google.com/file/d/1B6YNfCFRB4jqweVABamPQWlgziFlNIDK/view?usp=sharing>
