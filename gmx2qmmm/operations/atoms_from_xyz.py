"""
To change this license header, choose License Headers in Project
Properties. To change this template file, choose Tools | Templates and
open the template in the editor.

This will read an orca.xyz file and produce the matrix of coordinates.
Will default to angstrom. Unit can be set via A or B as second line
argument.
"""

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import re
import sys


def atoms_from_xyz(inpname):
    ifile = open(inpname)
    n_a = 0
    list_atoms = []
    for line in ifile:
        n_a = int(line)
        break
    for line in ifile:
        break
    for line in ifile:
        match = re.search(r"\s(\w+)\s+", line, flags=re.MULTILINE)
        list_atoms.append(match.group(1))
        if len(list_atoms) == n_a:
            break
    ifile.close()
    if len(list_atoms) != n_a:
        print("Did not find the same number of coordinates as indicated by number of atoms at beginning of file!")
        exit(1)
    return list_atoms


if __name__ == "__main__":
    atoms = atoms_from_xyz(sys.argv[1])
    for element in atoms:
        print(element)
