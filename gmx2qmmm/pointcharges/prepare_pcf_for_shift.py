"""
To change this license header, choose License Headers in Project
Properties. To change this template file, choose Tools | Templates and
open the template in the editor.

Will read a generate_pcf_from_XXX output, will prepare a
generate_charge_shift input requires input PCF, qmatoms(format below),
charge of QM region as in the QM calculation, connectivity list (format
below), output1 (qmcoords), output2 (m1list), output3(m2list),
output4(updated chargelist)

qmatoms format: contains all atom numbers in the qmregion. Remaining
charge difference after subtracting QM charge (sysargv3) will be split
between m1atoms in the output PCF connectivity list format: per line:
ATOMNR BONDPARTNERNR1 BONDPARTNERNR2 ...
"""

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # during a rain storm

import re

import numpy as np


def get_bondpartners(connlist, target):
    partnerlist = []
    for entry in connlist:
        found = False

        for i in range(0, len(entry)):
            if int(entry[i]) == int(target):
                found = True
                break
        if found:
            if int(entry[0]) == int(target):
                for i in range(1, len(entry)):
                    if int(entry[i]) not in np.array(partnerlist).astype(int):
                        partnerlist.append(int(entry[i]))
            else:
                if int(entry[0]) not in np.array(partnerlist).astype(int):
                    partnerlist.append(int(entry[0]))
    return partnerlist


def identify_m2(qmlist, m1list, connlist):
    m2list = []
    for element in m1list:
        m2line = []
        bondlist = get_bondpartners(connlist, element)
        for entry in bondlist:
            if int(entry) not in np.array(qmlist).astype(int):
                m2line.append(entry)
        m2list.append(m2line)
    return m2list


def identify_m1(qmlist, connlist):
    m1list = []
    for element in qmlist:
        bondlist = get_bondpartners(connlist, element)
        for entry in bondlist:
            if (int(entry) not in np.array(qmlist).astype(int)) and (
                int(entry) not in np.array(m1list).astype(int)
            ):
                m1list.append(entry)
    return m1list


def read_conn_list(inp):
    connlist = []
    with open(inp) as ifile:
        for line in ifile:
            connline = re.findall("\d+", line)
            if connline:
                connlist.append(connline)
    return np.array(connlist)


def read_qmatom_list(inp):
    qmatomlist = []
    with open(inp) as ifile:
        for line in ifile:
            if "[" in line or "]" in line:
                continue
            atomlist = re.findall("\d+", line)
            if atomlist:
                for element in atomlist:
                    qmatomlist.append(element)
    sortedlist = sorted(np.array(qmatomlist).astype(int))
    return sortedlist


def get_qmcoords(qmatoms, charges):
    qmcoordlist = []
    for element in qmatoms:
        qmline = [
            float(charges[int(element) - 1][0]),
            float(charges[int(element) - 1][1]),
            float(charges[int(element) - 1][2]),
            float(charges[int(element) - 1][3]),
        ]
        qmcoordlist.append(qmline)
    return qmcoordlist


def read_charge_list(inp):
    chargelist = []
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(
                r"^\s*([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                chargeline = [
                    float(match.group(1)),
                    float(match.group(2)),
                    float(match.group(3)),
                    float(match.group(4)),
                ]
                chargelist.append(chargeline)
    return chargelist


def eliminate_and_shift_to_m1(qmatoms, charges, m1list, qmcharge, qmcoordsq):
    updated_charges = []
    curr_charge = float(0.0)
    for element in qmcoordsq:
        curr_charge += float(element[3])
    count = 0
    parentcharge = float(curr_charge) - float(qmcharge)
    for element in charges:
        count += 1
        if count in np.array(qmatoms).astype(int):
            updated_charges.append(["QM"])
        else:
            chargeline = []
            for entry in element:
                chargeline.append(float(entry))
            updated_charges.append(chargeline)
    if len(m1list) != 0:
        parentcharge /= float(len(m1list))
    count = 0
    for element in updated_charges:
        count += 1
        if count in np.array(m1list).astype(int):
            updated_charges[count - 1][3] = float(
                updated_charges[count - 1][3]
            ) + float(parentcharge)
    return updated_charges


def prepare_pcf_for_shift_fieldsonly(charges, qmatomlist, qmcharge, connlist):
    m1list = identify_m1(qmatomlist, connlist)
    m2list = identify_m2(qmatomlist, m1list, connlist)
    qmcoordlist = get_qmcoords(qmatomlist, charges)
    updated_chargelist = eliminate_and_shift_to_m1(
        qmatomlist, charges, m1list, qmcharge, qmcoordlist
    )
    return qmcoordlist, m1list, m2list, updated_chargelist


def prepare_pcf_for_shift(
    inp, qmatoms, qmcharge, connfile, qmcoords, m1file, m2file, outfile
):
    chargelist = np.array(read_charge_list(inp))
    qmatomlist = read_qmatom_list(qmatoms)
    connlist = read_conn_list(connfile)
    m1list = identify_m1(qmatomlist, connlist)
    m2list = identify_m2(qmatomlist, m1list, connlist)
    qmcoordlist = get_qmcoords(qmatomlist, chargelist)
    updated_chargelist = eliminate_and_shift_to_m1(
        qmatomlist, chargelist, m1list, qmcharge, qmcoordlist
    )
    with open(qmcoords, "w") as ofile:
        for element in np.array(qmcoordlist):
            ofile.write(
                "{:<.10f} {:<.10f} {:<.10f}\n".format(
                    float(element[0]), float(element[1]), float(element[2])
                )
            )
    with open(m1file, "w") as ofile:
        for element in np.array(m1list):
            ofile.write(str(int(element)) + "\n")
    with open(m2file, "w") as ofile:
        for element in np.array(m2list):
            for entry in element:
                ofile.write(str(int(entry)) + " ")
            ofile.write("\n")
    with open(outfile, "w") as ofile:
        for element in np.array(updated_chargelist):
            if element[0] != "QM":
                ofile.write(
                    "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                        float(element[0]),
                        float(element[1]),
                        float(element[2]),
                        float(element[3]),
                    )
                )
            else:
                ofile.write("QM\n")
        ofile.write("$end\n")


if __name__ == "__main__":
    prepare_pcf_for_shift(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        sys.argv[4],
        sys.argv[5],
        sys.argv[6],
        sys.argv[7],
        sys.argv[8],
    )
