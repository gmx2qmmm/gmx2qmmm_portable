# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will read an orca.xyz file and produce the matrix of coordinates. Will default to angstrom. Unit can be set via A or B as second line argument.

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import sys


def geo_from_xyz(inpname, inppar):
    ifile = open(inpname)
    n_a = 0
    list_xyz = []
    for line in ifile:
        n_a = int(line)
        break
    for line in ifile:  # skip a line
        break
    for line in ifile:
        matchlist = re.findall("\s*([-]*\d+\.\d{14})", line)
        list_xyz += matchlist
    if (inppar) == "B":
        list_bohrs = []
        for i in range(0, n_a):
            for j in range(0, 3):
                list_bohrs.append(
                    "{:.14f}".format(float(list_xyz[i * 3 + j]) / 0.529177)
                )
        list_xyz = list_bohrs
    ifile.close()
    if len(list_xyz) != n_a * 3:
        print("geo_from_xyz: Did not find the same number of coordinates as indicated by 3*(number of atoms at beginning of file)!")
        print("geo_from_xyz: Alternatively, formatting may be off! Please provide an output from a converter script as input!")
        exit(1)
    return list_xyz


if __name__ == "__main__":
    coords = geo_from_xyz(sys.argv[1], sys.argv[2])
    for element in coords:
        print(element)
