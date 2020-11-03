# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will read an orca.hes file and produce the matrix of coordinates. Will default to angstrom. Unit can be set via A or B as second line argument.

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"


def geo_from_hes(inpname, param):
    import re

    ifile = open(inpname)
    n_a = 0
    list_xyz = []
    for line in ifile:
        match = re.search("\$atoms", line)
        if match:
            for line in ifile:
                n_a = int(line)
                break
            break
    for line in ifile:
        matchlist = re.findall("\s+([-]*\d+\.\d{12})", line)
        list_xyz += matchlist
        if len(list_xyz) == 3 * n_a:
            break
    list_ang = []
    for i in range(0, n_a):
        for j in range(0, 3):
            list_ang.append("{:.12f}".format(float(list_xyz[i * 3 + j]) * 0.529177))
    list_xyz = list_ang
    list_bohrs = []
    if (param) == "B":
        for i in range(0, n_a):
            for j in range(0, 3):
                list_bohrs.append(
                    "{:.12f}".format(float(list_xyz[i * 3 + j]) / 0.529177)
                )
        list_xyz = list_bohrs
    ifile.close()
    if len(list_xyz) != n_a * 3:
        print "Did not find the same number of coordinates as indicated by 3*(number of atoms at beginning of file)!"
        exit(1)
    return list_xyz


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 2:
        param = sys.argv[2]
    else:
        param = "A"
    geo_from_hes(sys.argv[1], param)
