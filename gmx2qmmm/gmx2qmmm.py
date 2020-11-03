"""
This will read a gmx data set (.g96 or .gro, .top) and a list of atoms
to be computed with QM in a QM/MM scheme. Will also read the
parameters of the QM, MM and QMMM calculations and an active region.
Will set up a QM/MM calculation with the QM and MM parameters of
choice. Will produce a log file containing details and the QM/MM
energy.
"""

__author__ = "jangoetze"
__date__ = "$06-Feb-2018 12:45:17$"

import re
import os
import subprocess
import sys

import numpy as np

from gmx2qmmm._helper import _flatten, logger, stepper
from gmx2qmmm._helper import get_linkatoms_ang, make_xyzq
from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.pointcharges import generate_charge_shift as final_pcf
from gmx2qmmm.operations import generate_top as topprep
from gmx2qmmm.operations import qmmm


def get_mollength_direct(molname, top):
    """"""
    mollength = 0
    with open(top) as ifile:
        # print str(top) + " is open"
        for line in ifile:
            match = re.search(r"^;", line, flags=re.MULTILINE)
            if match:
                continue
            match = re.search(r"^\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
            if match:
                # print "moltype was found"
                for line in ifile:
                    match = re.search(r"^;", line, flags=re.MULTILINE)
                    if match:
                        continue
                    else:
                        matchstring = r"^\s*" + re.escape(molname)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            # print str(matchstring) + " was found"
                            for line in ifile:
                                match = re.search(r"^;", line, flags=re.MULTILINE)
                                if match:
                                    continue
                                else:
                                    match = re.search(
                                        r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                    )
                                    if match:
                                        # print "atoms was found"
                                        for line in ifile:
                                            match = re.search(
                                                r"^;", line, flags=re.MULTILINE
                                            )
                                            if match:
                                                continue
                                            else:
                                                match = re.search(
                                                    r"^\s*(\d+)\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+[-]*\d+[\.]*[\d+]*",
                                                    line,
                                                    flags=re.MULTILINE,
                                                )
                                                if match:
                                                    mollength = int(match.group(1))
                                                match = re.search(
                                                    r"^\s*\n", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    break
                                        break
                            break
                break
    return mollength


def get_curr_top(molname, top, basedir):
    curr_top = top
    found = make_pcf.checkformol(molname, top)

    if not found:
        toplist = make_pcf.getincludelist(top)
        for element in toplist:
            found = make_pcf.checkformol(molname, element)
            if found:
                curr_top = element
                break
    if not found:
        print("No charges found for " + str(molname) + ". Exiting.")
        exit(1)
    return curr_top


def get_connlist(offset, molname, top):
    connlist = []
    with open(top) as ifile:
        for line in ifile:
            match = re.search(r"^;", line, flags=re.MULTILINE)
            if match:
                continue
            match = re.search(r"^\[ moleculetype \]", line, flags=re.MULTILINE)
            if match:
                for line in ifile:
                    match = re.search(r"^;", line, flags=re.MULTILINE)
                    if match:
                        continue
                    else:
                        matchstring = r"^\s*" + re.escape(molname)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            for line in ifile:
                                match = re.search(r"^;", line, flags=re.MULTILINE)
                                if match:
                                    continue
                                else:
                                    match = re.search(
                                        r"^\[ bonds \]", line, flags=re.MULTILINE
                                    )
                                    if match:
                                        for line in ifile:
                                            match = re.search(
                                                r"^;", line, flags=re.MULTILINE
                                            )
                                            if match:
                                                continue
                                            else:
                                                match = re.search(
                                                    r"^\s*(\d+)\s+(\d+)",
                                                    line,
                                                    flags=re.MULTILINE,
                                                )
                                                if match:
                                                    a = int(match.group(1)) + int(
                                                        offset
                                                    )
                                                    b = int(match.group(2)) + int(
                                                        offset
                                                    )
                                                    if a > b:
                                                        c = b
                                                        b = a
                                                        a = c
                                                    found = False
                                                    for element in connlist:
                                                        if int(element[0]) == a:
                                                            if int(b) not in np.array(
                                                                element
                                                            ).astype(int):
                                                                element.append(int(b))
                                                                found = True
                                                    if not found:
                                                        connlist.append(
                                                            [int(a), int(b)]
                                                        )
                                                match = re.search(
                                                    r"^\s*\n", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    break
                                        break
                            break
                break
    return connlist


def read_conn_list_from_top(top, mollist, basedir):
    count = 0
    connlist = []
    for element in mollist:
        curr_top = get_curr_top(element[0], top, basedir)
        for i in range(0, int(element[1])):
            mollength = get_mollength_direct(element[0], curr_top)
            connset = get_connlist(count, element[0], curr_top)
            connlist += connset
            count += int(mollength)
    return connlist


def read_pathparams(inp):
    # software path
    g16path = ""
    tmpath = ""
    orcapath = ""
    gmxpath = "gromacs-2019.1/bin/"
    # Gaussain,Turbomole,ORCA, GROMACS, gaussianCMD, TMCMD, orcaCMD, gmxCMD
    info = [
        "",
        "",
        "",
        "gromacs-2019.1/bin/",
        "rung16",
        "",
        "",
        "gmx19",
        "gromacs-2019.1/bin/",
        "gromacs-2019.1/share/gromacs/top/",
    ]
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(r"^g16path\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[0] = str(match.group(1))

            match = re.search(r"^tmpath\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[1] = str(match.group(1))

            match = re.search(r"^orcapath\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[2] = str(match.group(1))

            match = re.search(r"^gmxpath\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[3] = str(match.group(1))

            match = re.search(r"^g16cmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[4] = str(match.group(1))

            match = re.search(r"^tmcmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[5] = str(match.group(1))

            match = re.search(r"^orcacmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[6] = str(match.group(1))

            match = re.search(r"^gmxcmd\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[7] = str(match.group(1))

            match = re.search(r"^gmxtop_path\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[8] = str(match.group(1))
    return info


def read_qmparams(inp):
    info = [
        "G16",
        "BP86",
        "STO-3G",
        int(0),
        int(1),
        int(1),
        int(1000),
        "NONE",
    ]  # program method basis charge multiplicity cores memory(MB) extraopts
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(r"^program\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[0] = str(match.group(1)).upper()
            match = re.search(r"^method\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[1] = str(match.group(1)).upper()
            match = re.search(r"^basis\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[2] = str(match.group(1)).upper()
            match = re.search(r"^charge\s*=\s*([-]*\d+)", line, flags=re.MULTILINE)
            if match:
                info[3] = int(match.group(1))
            match = re.search(r"^multiplicity\s*=\s*(\d+)", line, flags=re.MULTILINE)
            if match:
                info[4] = int(match.group(1))
            match = re.search(r"^cores\s*=\s*(\d+)", line, flags=re.MULTILINE)
            if match:
                info[5] = int(match.group(1))
            match = re.search(r"^memory\s*=\s*(\d+)", line, flags=re.MULTILINE)
            if match:
                info[6] = int(match.group(1))
            match = re.search(r"^extra\s*=\s*(.+)$", line, flags=re.MULTILINE)
            if match:
                info[7] = match.group(1)
    return info


def read_mmparams(inp):
    info = ["amberGS", float(2.0)]
    flaglist = []
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(r"^ff\s*=\s*(\d*\.*\d*)", line, flags=re.MULTILINE)
            if match:
                info[0] = float(match.group(1))
                continue
            match = re.search(r"^rvdw\s*=\s*(\d*\.*\d*)", line, flags=re.MULTILINE)
            if match:
                info[1] = float(match.group(1))
                continue
            match = re.search(r"\-D(\S*)", line, flags=re.MULTILINE)
            if match:
                if match.group(1) not in flaglist:
                    flaglist.append((match.group(1)).upper())
    return info, flaglist


def read_qmmmparams(inp):
    jobname = "testjob"
    jobtype = "singlepoint"
    propa = "steep"
    info = [
        jobname,
        jobtype.upper(),
        float(0.00001),
        int(5),
        float(0.1),
        propa.upper(),
        float(0.0018897261),
        int(0),
        "MORSE",
        "YES",
    ]
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(r"^jobname\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[0] = match.group(1)
            match = re.search(r"^jobtype\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[1] = str(match.group(1)).upper()
            match = re.search(r"^f_thresh\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[2] = float(match.group(1))
            match = re.search(r"^maxcycle\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[3] = int(match.group(1))
            match = re.search(r"^initstep\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[4] = float(match.group(1))
            match = re.search(r"^propagator\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[5] = str(match.group(1)).upper()
            match = re.search(r"^disp\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[6] = float(match.group(1))
            match = re.search(r"^current_step\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[7] = int(match.group(1))
            match = re.search(r"^databasefit\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                # info[8]=match.group(1)
                info[8] = str(match.group(1)).upper()
            match = re.search(r"^optlastonly\s*=\s*(\S+)", line, flags=re.MULTILINE)
            if match:
                info[9] = str(match.group(1)).upper()
    return info


def get_pcf_self_grad(charges):
    charges_grad = []
    for i in range(0, len(charges) - 1):
        grad_line = np.array([0.0, 0.0, 0.0])
        curr_vec = np.array(
            [float(charges[i][0]), float(charges[i][1]), float(charges[i][2])]
        )
        curr_charge = float(charges[i][3])
        for j in range(i + 1, len(charges)):
            new_vec = np.arrayay(
                [float(charges[j][0]), float(charges[j][1]), float(charges[j][2])]
            )
            new_charge = float(charges[j][3])
            diff = np.array(new_vec - curr_vec)
            dist = np.linalg.norm(diff)
            for k in range(0, 3):
                grad_line[k] += (curr_charge * new_charge / (dist * dist)) * (
                    diff[k] / dist
                )
        charges_grad.append(grad_line)
    return np.array(charges_grad)


def get_pcf_self_pot(charges):
    charges_pot = 0.0
    for i in range(0, len(charges)):
        curr_vec = np.array(
            [float(charges[i][0]), float(charges[i][1]), float(charges[i][2])]
        )
        curr_charge = float(charges[i][3])
        for j in range(i + 1, len(charges)):
            new_vec = np.array(
                [float(charges[j][0]), float(charges[j][1]), float(charges[j][2])]
            )
            new_charge = float(charges[j][3])
            diff = np.array(new_vec - curr_vec)
            dist = np.linalg.norm(diff)
            charges_pot += curr_charge * new_charge / dist
    return charges_pot


def read_charges_clean(file):
    clean_charges = []
    with open(file) as ifile:
        for line in ifile:
            match = re.search(
                r"^\s*([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                clean_charges.append(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                        float(match.group(4)),
                    ]
                )
    return clean_charges


def ang2bohr(xyzqlist):
    xyzqlist_shifted = []
    a2b = 1.8897261 * (-78.3363167335) / (-78.3363157152)  # as is gaussian16
    for element in xyzqlist:
        new_element = []
        for i in range(0, len(element)):
            if i < 3:
                new_element.append(float(element[i]) * a2b)
            else:
                new_element.append(float(element[i]))
        xyzqlist_shifted.append(new_element)
    return xyzqlist_shifted


def write_charge_field(charges, outname):
    with open(outname, "w") as ofile:
        for element in charges:
            ofile.write(
                "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                    float(element[0]),
                    float(element[1]),
                    float(element[2]),
                    float(element[3]),
                )
            )
        ofile.write("$end\n")


def get_linkcorrlist(linkatoms, qmatomlist, m1list, m2list, connlist):
    linkcorrlist = []
    m3list = []
    m4list = []
    q1list = []
    q2list = []
    q3list = []
    # get q1 and m2
    for element in m1list:
        q1line = []
        for entry in connlist:
            if int(element) in np.array(entry).astype(int):
                if int(element) != int(entry[0]):
                    if int(entry[0]) in np.array(qmatomlist).astype(int):
                        q1line.append(int(entry[0]))
                else:
                    for i in range(1, len(entry)):
                        if int(entry[i]) in np.array(qmatomlist).astype(int):
                            q1line.append(int(entry[i]))
        q1list.append(q1line)
    # get q2
    q1list = list(_flatten(q1list))
    for element in q1list:
        q2line = []
        for conn in connlist:
            if int(element) in np.array(conn).astype(int):
                if (
                    int(element) != int(conn[0])
                    and (int(conn[0]) in np.array(qmatomlist).astype(int))
                    and (int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)]))
                ):
                    q2line.append(int(conn[0]))
                elif int(element) == int(conn[0]):
                    for i in range(1, len(conn)):
                        if (int(conn[i]) in np.array(qmatomlist).astype(int)) and (
                            int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                        ):
                            q2line.append(int(conn[i]))
        q2list.append(q2line)
    # get q3
    for element in q2list:
        q3lineline = []
        for entry in element:
            q3line = []
            for conn in connlist:
                if int(entry) in np.array(conn).astype(int):
                    if (
                        int(entry) != int(conn[0])
                        and (int(conn[0]) in np.array(qmatomlist).astype(int))
                        and int(conn[0]) not in np.array(q1list).astype(int)
                        and int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)])
                    ):
                        q3line.append(int(conn[0]))
                    elif int(entry) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if (
                                int(conn[i]) in np.array(qmatomlist).astype(int)
                                and int(conn[i]) not in np.array(q1list).astype(int)
                                and int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                            ):
                                q3line.append(int(conn[i]))
            q3lineline.append(q3line)
        q3list.append(q3lineline)
    # get m3
    for element in m2list:
        m3lineline = []
        for entry in element:
            m3line = []
            for conn in connlist:
                if int(entry) in np.array(conn).astype(int):
                    if (
                        int(entry) != int(conn[0])
                        and (int(conn[0]) not in np.array(qmatomlist).astype(int))
                        and int(conn[0]) not in np.array(m1list).astype(int)
                    ):
                        m3line.append(int(conn[0]))
                    elif int(entry) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if int(conn[i]) not in np.array(qmatomlist).astype(
                                int
                            ) and int(conn[i]) not in np.array(m1list).astype(int):
                                m3line.append(int(conn[i]))
            m3lineline.append(m3line)
        m3list.append(m3lineline)
    # get m4
    for element in m3list:
        m4linelineline = []
        for entry in element:
            m4lineline = []
            for stuff in entry:
                m4line = []
                for conn in connlist:
                    if int(stuff) in np.array(conn).astype(int):
                        if (
                            int(stuff) != int(conn[0])
                            and (int(conn[0]) not in np.array(qmatomlist).astype(int))
                            and int(conn[0]) not in np.array(m1list).astype(int)
                        ):
                            found = False
                            for morestuff in m2list:
                                if int(conn[0]) in np.array(morestuff).astype(int):
                                    found = True
                                    break
                            if not found:
                                m4line.append(int(conn[0]))
                        elif int(stuff) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if int(conn[i]) not in np.array(qmatomlist).astype(
                                    int
                                ) and int(conn[i]) not in np.array(m1list).astype(int):
                                    found = False
                                    for morestuff in m2list:
                                        if int(conn[i]) in np.array(morestuff).astype(
                                            int
                                        ):
                                            found = True
                                            break
                                    if not found:
                                        m4line.append(int(conn[i]))
                m4lineline.append(m4line)
            m4linelineline.append(m4lineline)
        m4list.append(m4linelineline)
    # set up link atom charge corr pairs: q1-m2, q1-m3, q2-m2, l-m2, l-m3, l-m4 - l are represented as their m1 counterparts!
    count = 0
    for element in m1list:
        linkpairline = []
        for entry in m2list[count]:
            if int(element) < int(entry):
                linkpairline.append([element, entry])
            else:
                linkpairline.append([entry, element])
        for entry in _flatten(m3list[count]):
            if int(element) < int(entry):
                linkpairline.append([element, entry])
            else:
                linkpairline.append([entry, element])
        for entry in _flatten(m4list[count]):
            if int(element) < int(entry):
                linkpairline.append([element, entry])
            else:
                linkpairline.append([entry, element])
        linkcorrlist.append(linkpairline)
        count += 1
    count = 0
    for element in q1list:
        linkpairline = []
        for stuff in m2list[count]:
            if int(element) < int(stuff):
                linkpairline.append([element, stuff])
            else:
                linkpairline.append([stuff, element])
        for stuff in list(_flatten(m3list[count])):
            if int(element) < int(stuff):
                linkpairline.append([element, stuff])
            else:
                linkpairline.append([stuff, element])
        linkcorrlist.append(linkpairline)
        count += 1
    count = 0
    for element in q2list:
        for entry in element:
            linkpairline = []
            for stuff in m2list[count]:
                if int(entry) < int(stuff):
                    linkpairline.append([entry, stuff])
                else:
                    linkpairline.append([stuff, entry])
            linkcorrlist.append(linkpairline)
        count += 1
    reshaped_linkcorrlist = np.array(list(_flatten(linkcorrlist))).reshape(-1, 2)
    final_linkcorrlist = []
    for element in reshaped_linkcorrlist:
        found = False
        for i in range(0, len(final_linkcorrlist)):
            if (
                final_linkcorrlist[i][0] == element[0]
                and final_linkcorrlist[i][1] == element[1]
            ):
                found = True
                break
        if found:
            continue
        final_linkcorrlist.append(element)
    return (
        sorted(final_linkcorrlist, key=lambda l: l[0]),
        q1list,
        q2list,
        q3list,
        m3list,
    )


def write_highprec(gro, jobname, logfile):
    filename = str(jobname + ".g96")
    with open(filename, "w") as ofile:
        with open(gro) as ifile:
            count = 0
            for line in ifile:
                count += 1
                if count == 2:
                    break
            ofile.write("TITLE\nProtein\nEND\nPOSITION\n")
            counter = 0
            finalline = ""
            for line in ifile:
                match = re.search(
                    r"^([\s\d]{5})(.{5})(.{5})([\s\d]{5})\s*([-]*\d+\.\d*)\s*([-]*\d+\.\d*)\s*([-]*\d+\.\d*)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    logger(
                        logfile,
                        str(
                            "Successfully wrote "
                            + str(int(counter))
                            + " atoms to internal precision format file.\n"
                        ),
                    )
                    finalline = line
                    break
                else:
                    counter += 1
                    ofile.write(
                        str(match.group(1))
                        + " "
                        + str(match.group(2))
                        + " "
                        + str(match.group(3))
                        + " {:>6d} ".format(int(counter))
                        + "{:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                            float(match.group(5)),
                            float(match.group(6)),
                            float(match.group(7)),
                        )
                    )

            ofile.write("END\nBOX\n")
            match = re.search(
                r"^\s*(\d+\.*\d*)\s*(\d+\.*\d*)\s*(\d+\.*\d*)",
                finalline,
                flags=re.MULTILINE,
            )
            if not match:
                logger(
                    logfile,
                    str(
                        "Unexpected line instead of box vectors. Exiting. Last line:\n"
                    ),
                )
                logger(logfile, line)
                exit(1)
            else:
                ofile.write(
                    str(
                        " {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                        )
                    )
                )
            ofile.write("END")
    return filename


def write_m1_index(jobname, m1list):
    m1listfile = str(jobname + ".m1.ndx")
    with open(m1listfile, "w") as ofile:
        ofile.write("[ M1 ]\n")
        for element in m1list:
            ofile.write(str(str(int(element)) + " "))


def display_help():
    print(
        "gmx2qmmm (GNU licenced) by Jan Philipp Goetze\n"
        "gmx2qmmm is a free QM/MM interface for Gromacs. Compatible with most versions, starting from 5.X.\n"
        "gmx2qmmm uses an additive QM/MM scheme with a numerical charge shift procedure and correction of cut bonds via model potentials.\n"
        "\nSYNOPSIS\n"
        "\ngmx2qmmm <options>\n"
        "\nOPTIONS\n"
        "\nOptions to specify input files:\n\n"
        " -c      [<.g96>]           (conf.g96)\n"
        "           Structure file: g96\n"
        " -p      [<.top>]           (topol.top)\n"
        "           Topology file\n"
        " -n      [<.ndx>]           (qmatoms.ndx)\n"
        "           Index file (containing only the QM atom indices)\n"
        " -qm     [<.dat>]           (qm.dat)\n"
        "           QM input parameters\n"
        " -mm     [<.dat>]           (mm.dat)\n"
        "           MM input parameters\n"
        " -qmmm   [<.dat>]           (qmmm.dat)\n"
        "           QM/MM input parameters\n"
        " -path   [<.dat>]           (path.dat)\n"
        "           QM & MM software path input parameters\n"
        " -act    [<.ndx>]           (act.ndx)        (Opt.)\n"
        "           Index file (containing only indices of non-frozen atoms)\n"
        "\nOptions to specify output files:\n\n"
        " -g      [<.log>]           (logfile)   (Opt.)\n"
        "           Log file\n"
        "\ngmx2qmmm reminds you: This is not an official Gromacs program, so no reminder!\n\n"
    )


def find_optentry(command, options):
    workwith = options
    defaults = [
        "-c",
        "conf.g96",
        "-p",
        "topol.top",
        "-n",
        "qmatoms.ndx",
        "-qm",
        "qm.dat",
        "-mm",
        "mm.dat",
        "-qmmm",
        "qmmm.dat",
        "-act",
        "act.ndx",
        "-g",
        "logfile",
        "-path",
        "path.dat",
    ]
    if command not in options:
        workwith = defaults
    else:
        for i in range(0, len(options)):
            if options[i] == command:
                if i + 1 >= len(options):
                    workwith = defaults
                    break
                if options[i + 1][0] == "-":
                    workwith = defaults
                    break
    for i in range(0, len(workwith) - 1):
        if workwith[i] == command:
            return workwith[i + 1]


def read_options(cmd_options):
    new_cmd_options = []
    for i in range(0, len(cmd_options)):
        new_cmd_options.append(cmd_options[i])
    if "-h" in new_cmd_options:
        display_help()
        exit(0)
    gro = find_optentry("-c", new_cmd_options)
    top = find_optentry("-p", new_cmd_options)
    qmatoms = find_optentry("-n", new_cmd_options)
    qmparams = find_optentry("-qm", new_cmd_options)
    mmparams = find_optentry("-mm", new_cmd_options)
    qmmmparams = find_optentry("-qmmm", new_cmd_options)
    activefile = find_optentry("-act", new_cmd_options)
    logfile = find_optentry("-g", new_cmd_options)
    path = find_optentry("-path", new_cmd_options)
    return gro, top, qmatoms, qmparams, mmparams, qmmmparams, activefile, logfile, path


def add_gmxpath(basedir, path):
    readfile = open(basedir + "/pointcharges/generate_pcf_from_top.py", "r")
    gmxpath = 'gmxpath="' + path + '"'
    oline = readfile.readlines()
    insert_index = oline.index("#GROMACS path\n") + 1
    oline[insert_index] = gmxpath + "\n"
    readfile.close()

    src = open(basedir + "/pointcharges/generate_pcf_from_top.py", "w")
    src.writelines(oline)
    src.close()


def gmx2qmmm(cmd_options):
    (
        gro,
        top,
        qmatoms,
        qmparams,
        mmparams,
        qmmmparams,
        activefile,
        logfile,
        path,
    ) = read_options(cmd_options)
    basedir = os.path.dirname(os.path.abspath(__file__))
    if os.path.isfile(logfile):
        subprocess.call(["rm", logfile])
    logger(logfile, "Reading QM/MM parameters...")
    qmmminfo = read_qmmmparams(qmmmparams)
    logger(logfile, "done.\n")
    logger(logfile, "Reading QM/MM software path...\n")
    pathinfo = read_pathparams(path)

    if pathinfo[0] != "":
        logger(logfile, "The path of Gaussain is: %s \n" % pathinfo[0])
    if pathinfo[1] != "":
        logger(logfile, "The path of Turbomole is: %s \n" % pathinfo[1])
    if pathinfo[2] != "":
        logger(logfile, "The path of ORCA is: %s \n" % pathinfo[2])
    if pathinfo[3] != "":
        logger(logfile, "The path of GROMACS is: %s \n" % pathinfo[3])
        gmxpath = pathinfo[3]
    if pathinfo[4] != "":
        logger(logfile, "The executed command of Gaussain is: %s \n" % pathinfo[4])
    if pathinfo[5] != "":
        logger(logfile, "The executed command of Turbomole is: %s \n" % pathinfo[5])
    if pathinfo[6] != "":
        logger(logfile, "The executed command of ORCA is: %s \n" % pathinfo[6])
    if pathinfo[7] != "":
        logger(logfile, "The executed command of GROMACS is: %s \n" % pathinfo[7])
    if pathinfo[8] != "":
        logger(logfile, "The path of force field in GROMACS is: %s \n" % pathinfo[8])

    logger(logfile, "done.\n")

    jobname = qmmminfo[0]
    jobtype = qmmminfo[1]
    step = qmmminfo[7]
    gro = stepper(gro, step)

    logger(logfile, "Initializing dependencies...")

    logger(logfile, "complete.\n")
    chargevec = []
    logger(logfile, "Trying to understand your MM files.\n")
    logger(logfile, "List of molecules...")
    mollist = make_pcf.readmols(top)
    logger(logfile, "done.\n")
    logger(logfile, "Reading charges...")
    for element in mollist:
        chargevec.extend(make_pcf.readcharges(element, top, pathinfo))
    logger(logfile, "done.\n")

    if gro[-4:] == ".gro":
        logger(logfile, "Reading geometry (.gro)...")
        geo = make_pcf.readgeo(gro)
    elif gro[-4:] == ".g96":
        logger(logfile, "Reading geometry (.g96)...")
        geo = make_pcf.readg96(gro)

    logger(logfile, "%s\n" % geo)
    logger(logfile, "done.\n")
    logger(logfile, "Reading connectivity matrix...")
    connlist = read_conn_list_from_top(top, mollist, basedir)
    logger(logfile, "done.\n")
    logger(logfile, "Trying to understand your QM, MM and QM/MM parameters.\n")
    logger(logfile, "Reading QM atom list...")
    qmatomlist = prep_pcf.read_qmatom_list(qmatoms)
    logger(logfile, "done.\n")
    logger(logfile, "Reading QM parameters...")
    qminfo = read_qmparams(qmparams)
    logger(logfile, "done.\n")
    logger(logfile, "Reading MM parameters...")
    mminfo, flaglist = read_mmparams(mmparams)
    logger(logfile, "done.\n")

    if gro[-4:] == ".gro":
        logger(logfile, "Writing high-precision coordinate file...")
        grohigh = write_highprec(gro, jobname, logfile)
        gro = jobname + ".g96"
        logger(logfile, "done.\n")

    logger(
        logfile,
        'Setting up a calculation named "'
        + str(jobname)
        + '" of type "'
        + str(jobtype)
        + '".\n',
    )
    qmmmtop = str(jobname + ".qmmm.top")
    pcffile = str(jobname + ".pointcharges")
    pcffile = stepper(pcffile, step)
    logger(
        logfile,
        "Starting the preparation of the point charge field used in the calculation.\n",
    )
    logger(logfile, "Creating the full xyzq Matrix...")
    xyzq = make_xyzq(geo, chargevec)
    logger(logfile, "done.\n")
    logger(
        logfile,
        "Preparing the point charge field for a numerically optimized charge shift...",
    )
    (
        qmcoordlist,
        m1list,
        m2list,
        updated_chargelist,
    ) = prep_pcf.prepare_pcf_for_shift_fieldsonly(xyzq, qmatomlist, qminfo[3], connlist)
    logger(logfile, "done.\n")
    logger(logfile, "Setting up link atoms...")
    linkatoms = get_linkatoms_ang(xyzq, qmatomlist, m1list, connlist, [])
    linkcorrlist, q1list, q2list, q3list, m3list = get_linkcorrlist(
        linkatoms, qmatomlist, m1list, m2list, connlist
    )
    logger(logfile, "done.\n")
    if qmmminfo[7] > 0:
        jobname = str(jobname + "." + str(qmmminfo[7]))
    if not os.path.isfile(str(jobname + ".pointcharges")):
        logger(logfile, "Shifting...")
        final_pcf.generate_charge_shift_fieldsonly(
            updated_chargelist, m1list, qmcoordlist, m2list, jobname, basedir
        )
    else:
        logger(
            logfile,
            "NOTE: Shifting omitted due to "
            + str(jobname + ".pointcharges")
            + " being an existing file!\n",
        )
    logger(logfile, "done.\n")
    logger(logfile, "Preparing the QM/MM top file...")
    topprep.generate_top_listsonly(
        top,
        qmatoms,
        qmmmtop,
        flaglist,
        q1list,
        q2list,
        q3list,
        m1list,
        m2list,
        m3list,
        basedir,
        logfile,
        pathinfo,
    )
    logger(logfile, "done.\n")
    active = []
    if jobtype != "SINGLEPOINT":
        logger(logfile, "Reading indices of active atoms...")
        active = prep_pcf.read_qmatom_list(activefile)
        logger(logfile, "done.\n")
    logger(logfile, "Starting the job.\n")
    qmmm.perform_job(
        gro,
        top,
        xyzq,
        connlist,
        qmatomlist,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        jobtype,
        linkcorrlist,
        flaglist,
        pcffile,
        linkatoms,
        active,
        logfile,
        basedir,
        step,
        pathinfo,
    )


if __name__ == "__main__":
    gmx2qmmm(sys.argv)
