# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will read an orca.hess file and report the full hessian.

__author__ = "jangoetze"
__date__ = "$02-Jan-2018 14:45:17$"

import re


def hessian_from_hes(inpname):
    n_dof = 0
    curr_col = 0
    list_hess = []
    ifile = open(inpname)
    for line in ifile:
        match = re.search("\$hessian", line)
        if match:
            for line in ifile:
                n_dof = int(line)
                break
    ifile.close()
    matchlist = []
    for i in range(0, n_dof):
        ifile = open(inpname)
        for line in ifile:
            match = re.search("\$hessian", line)
            if match:
                curr_col = 0
                for line in ifile:
                    match = re.search("\s{6}\s*\d+", line)
                    if match:
                        for line in ifile:
                            match = re.findall(
                                r"^\s*" + str(i) + "\s+", line, flags=re.MULTILINE
                            )
                            if match:
                                matchlist = re.findall(
                                    "\s+([-]*\d\.\d{10}E[+-]\d\d)", line
                                )
                                list_hess += matchlist
                                break
                        curr_col += len(matchlist)
                        if curr_col == n_dof:
                            break
                break
        ifile.close()
    matrix_hess = []
    for i in range(0, n_dof):
        line_hess = []
        for j in range(0, n_dof):
            line_hess.append(float(list_hess[i * n_dof + j]))
        matrix_hess.append(line_hess)
    if len(list_hess) != n_dof * n_dof:
        print("Did not find the same number of Hessian entries as indicated by (degrees of freedom)**2!")
        exit(1)
    return matrix_hess


if __name__ == "__main__":
    hessian_from_hes(name1)  # name1 is not defined
