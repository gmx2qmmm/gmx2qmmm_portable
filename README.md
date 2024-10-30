gmx2qmmm v2 (in preparation)
============================

**gmx2qmmm** is a python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculation.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents
- [gmx2qmmm v2 (in preparation)](#gmx2qmmm-v2-in-preparation)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [System requirements](#system-requirements)
  - [**gmx2qmmm** job](#gmx2qmmm-job)
  - [Input files](#input-files)
  - [Installation](#installation)
  - [Examples](#examples)
  - [References](#references)
  - [Support and development](#support-and-development)
  - [Links](#links)
<!-- END doctoc generated TOC please keep comment here to allow auto update -->

---

## Overview

**gmx2qmmm** is a Python package to bridge [Gaussian] and [Gromacs]. The test runs were performed using [Gaussian16] and [Gromacs 5.0.2], but the code should be able to read earlier Gaussian and other Gromacs versions. The only limits are the formats of the human-readable input and output files of each program, as such, conversion scripts can be written to make the interface work with any version, if the current code does not support it.

Conceptually, **gmx2qmmm** creates a QM/MM potential and performs either single point calculations (i.e., the current energy of your system), geometry optimizations, and linear relaxed scan. (Other ultilities are ongoing)

## System requirements
 - [Python 3.6+]
 - [Gaussian16] (and earlier version)
 - [Gromacs 5.0.2] (and earlier version)

## **gmx2qmmm** job
|Job type|Calculation|
| ------ | ------ |
|Single point calcuation (SP)|Calculate single point energy and forces (.xyz) |
|Geometry optimizations (OPT)|Optimize the system energy via optimizer ([Steepest descent], [Conjugate gradient] or [BFGS])|
|Relaxed Scan (SCAN)|Relaxed linear scan (angle and dihedral angle are in development)|

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

For advance information please read the [Documentation].

---

## Installation

Clone the repository and install with `pip`. We recommend using a fresh virtual environment.

```bash
~ $ git clone https://github.com/gmx2qmmm/gmx2qmmm_portable
~ $ cd gmx2qmmm
~/gmx2qmmm_portable $ pip install .
```

For development installation, consider installing in editable mode:

```bash
~/gmx2qmmm_portable $ pip install -e .
```

---

## Examples
The directory example contains SP, OPT, linear SCAN calculation for glycine serine (GLYSER).
![alt text](https://github.com/gmx2qmmm/gmx2qmmm_portable/blob/master/examples/glyser.png?raw=true)

The names of the input files are default. Go to the sp/opt directory on the command line and run:

```
$ gmx2qmmm
```
---

## References

> A user‐friendly, Python‐based quantum mechanics/Gromacs interface: gmx2qmmm
> Jan P. Götze, Yuan‐Wei Pi, Simon Petry, Fabian Langkabel, Jan Felix Witte, Oliver Lemke
> https://doi.org/10.1002/qua.26486

## Support and development
For bug reports/suggestions/complaints please raise an issue on [GitHub].

Or contact us directly: [gmx2qmmm@gmail.com]

## Links
- [Documentation]
- [gmx2qmmm reference]


[python 3.6+]:<https://docs.python.org/3.6>
[Gaussian16]:<https://gaussian.com/gaussian16/>
[Gromacs 5.0.2]:<http://www.gromacs.org>
[Gaussian]:<https://gaussian.com/gaussian16/>
[Gromacs]:<http://www.gromacs.org>
[GitHub]:<https://github.com/gmx2qmmm/gmx2qmmm_portable>
[gmx2qmmm@gmail.com]:<mailto:gmx2qmmm@gmail.com>
[Documentation]:<https://gmx2qmmm.github.io/gmx2qmmm_io>
[gmx2qmmm reference]:<https://drive.google.com/file/d/1B6YNfCFRB4jqweVABamPQWlgziFlNIDK/view?usp=sharing>
