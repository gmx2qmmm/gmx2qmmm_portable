"""
This will take a point charge field (turbomole style, without a leading
"$point_charges" line) and shift the charges such that the M1 atom
charge is replaced by a distribution of charge on the neighbouring
atoms. Will also generate a corrected dipole field at all Q atoms by
setting pairs of displaced dipoles next to the M2 atoms and numerically
optimizing them. Requires input point_charge field, M1 atom list, file
containing the XYZ coordinates of all QM atoms (one per row), M2 atoms
list(format see below), output file format for M2 file: each row
corresponds to M1 file rows, each row contains all M2 atom indices
separated by a space/tab whatever

FOR A QM/MM PCF, YOU NEED TO FIRST CUT THE QM ATOMS AND MOVE THEIR
SUMMED UP CHARGE TO THE M1 ATOM, UNLESS THE CHARGE IS REPRESENTED IN THE
QM PART
"""

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$"  # during a rain storm

import math
import os
import re
import sys

import numpy as np

from gmx2qmmm._helper import _flatten
from gmx2qmmm.pointcharges import sum_pcf_tm as make_p_sum


def uvec(vec):
    v = np.array(vec) / np.linalg.norm(np.array(vec))
    return v


def create_corr_charges(m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist):

    corr_charge_list = []
    count = 0
    for i in range(0, len(m2_nolist)):
        short_disp_vec = []
        for j in range(0, len(m2_nolist[i])):
            short_disp_vec.extend([disp_vec[count + j]])
        charge_list = write_disp_charges(
            m1coordsq[i], m2coordsqlist[i], short_disp_vec, disp_charge_vec[i]
        )
        corr_charge_list.extend(charge_list)
        count += len(m2_nolist[i])
    corr_charge = np.array(corr_charge_list).reshape(-1, 4)
    return corr_charge


def write_disp_charges(m1, m2coordsqlist, disp_vec, dispcharge):
    charge_list = []
    for i in range(0, len(m2coordsqlist)):
        ab_vec = np.array(
            [
                float(m2coordsqlist[i][0]),
                float(m2coordsqlist[i][1]),
                float(m2coordsqlist[i][2])
            ]
        ) - np.array([float(m1[0]), float(m1[1]), float(m1[2])])

        length_ab = np.linalg.norm(ab_vec)
        normvec = uvec(ab_vec)
        ab_long = np.array(normvec) * (length_ab + disp_vec[i])
        ab_short = np.array(normvec) * (length_ab - disp_vec[i])
        new_dispcharge = 0.0
        if disp_vec[i] >= 0.000001:
            new_dispcharge = 0.5 - (0.5 / (2.0 * disp_vec[i] + 1.0))
        if dispcharge > 0.0:
            chargeminus = [
                float(m1[0]) + float(ab_long[0]),
                float(m1[1]) + float(ab_long[1]),
                float(m1[2]) + float(ab_long[2]),
                float(-1.0 * new_dispcharge),
            ]
            chargeplus = [
                float(m1[0]) + float(ab_short[0]),
                float(m1[1]) + float(ab_short[1]),
                float(m1[2]) + float(ab_short[2]),
                float(new_dispcharge),
            ]
        else:
            chargeplus = [
                float(m1[0]) + float(ab_long[0]),
                float(m1[1]) + float(ab_long[1]),
                float(m1[2]) + float(ab_long[2]),
                float(new_dispcharge),
            ]
            chargeminus = [
                float(m1[0]) + float(ab_short[0]),
                float(m1[1]) + float(ab_short[1]),
                float(m1[2]) + float(ab_short[2]),
                float(-1.0 * new_dispcharge),
            ]
        charge_list.append(chargeminus)
        charge_list.append(chargeplus)
    return charge_list


def make_new_field(m2coordsqlist, corr_charge_list):
    new_field = []
    for element in m2coordsqlist:
        for entry in element:
            new_field.append(entry)
    for element in corr_charge_list:
        new_field.append(element)
    return new_field


def write_new_pcf(
    inp, out, m1line, m1, m2list, m2atoms, disp_vec, dispcharge, distrib_charge
):
    with open(inp) as ifile:
        ofile = open(out, "w")
        count = 0
        for line in ifile:
            count += 1
            if int(count) == int(m1line):
                ofile.write("QM\n")
                continue
            if int(count) in m2list:
                for i in range(0, len(m2list)):
                    if int(count) == int(m2list[i]):
                        ofile.write(
                            str(m2atoms[i][0])
                            + " "
                            + str(m2atoms[i][1])
                            + " "
                            + str(m2atoms[i][2])
                            + " "
                            + str(m2atoms[i][3])
                            + "\n"
                        )
                        continue
            match = re.search("\$end", line)
            if match:
                for i in range(0, len(m2list)):
                    # print m2atoms[i]
                    ab_vecq = np.array(m2atoms[i]) - np.array(m1)
                    ab_vec = np.array([ab_vecq[0], ab_vecq[1], ab_vecq[2]])
                    length_ab = np.linalg.norm(ab_vec)
                    normvec = uvec(ab_vec)
                    ab_long = np.array(normvec) * (length_ab + disp_vec[i])
                    ab_short = np.array(normvec) * (length_ab - disp_vec[i])
                    if distrib_charge > 0.0:
                        ofile.write(
                            str(float(m1[0]) + float(ab_long[0]))
                            + " "
                            + str(float(m1[1]) + float(ab_long[1]))
                            + " "
                            + str(float(m1[2]) + float(ab_long[2]))
                            + " "
                            + str(dispcharge)
                            + "\n"
                        )
                        ofile.write(
                            str(float(m1[0]) + float(ab_short[0]))
                            + " "
                            + str(float(m1[1]) + float(ab_short[1]))
                            + " "
                            + str(float(m1[2]) + float(ab_short[2]))
                            + " "
                            + str(-1.0 * dispcharge)
                            + "\n"
                        )

                    else:
                        ofile.write(
                            str(float(m1[0]) + float(ab_long[0]))
                            + " "
                            + str(float(m1[1]) + float(ab_long[1]))
                            + " "
                            + str(float(m1[2]) + float(ab_long[2]))
                            + " "
                            + str(-1 * dispcharge)
                            + "\n"
                        )
                        ofile.write(
                            str(float(m1[0]) + float(ab_short[0]))
                            + " "
                            + str(float(m1[1]) + float(ab_short[1]))
                            + " "
                            + str(float(m1[2]) + float(ab_short[2]))
                            + " "
                            + str(dispcharge)
                            + "\n"
                        )
                ofile.write("$end\n")
            else:
                ofile.write(line)
        ofile.close()


def write_new_field_to_disk_listsonly(inp, ofilename, new_field, getlist, m2_nolist):
    ofile = open(ofilename, "w")
    count = 0
    m2list = np.array(m2_nolist).reshape(-1)
    m2count = 0
    for element in inp:
        count += 1
        if int(count) in np.array(list(_flatten(m2list))).astype(int):
            if len(element) != 4:
                print("Line " + str(count) + " does not contain data. Exiting.")
                exit(1)
            ofile.write(
                "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                    float(element[0]),
                    float(element[1]),
                    float(element[2]),
                    float(new_field[m2count][3]) + float(element[3]),
                )
            )
            m2count += 1
            continue
        if int(count) in np.array(getlist).astype(int):
            ofile.write("QM\n")
            continue
        else:
            for i in range(0, len(element)):
                ofile.write(str(element[i]))
                if i != len(element) - 1:
                    ofile.write(" ")
            ofile.write("\n")
    for i in range(len(list(_flatten(m2list))), len(new_field)):
        ofile.write(
            "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                float(new_field[i][0]),
                float(new_field[i][1]),
                float(new_field[i][2]),
                float(new_field[i][3]),
            )
        )
    ofile.write("$end\n")
    ofile.close()


def write_new_field_to_disk(inp, ofilename, new_field, getlist, m2_nolist):
    ofile = open(ofilename, "w")
    count = 0
    m2list = np.array(m2_nolist).reshape(-1)

    with open(inp) as ifile:
        for line in ifile:
            count += 1
            if count in m2list.astype(int):
                match = re.search(
                    r"^\s*([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    print("Line " + str(count) + " does not contain data. Exiting.")
                    exit(1)
                m2pos = 0
                for i in range(0, len(m2list)):
                    if count == int(m2list[i]):
                        m2pos = i
                        break
                ofile.write(
                    "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                        float(new_field[i][3]) + float(match.group(4)),
                    )
                )
                continue
            if int(count) in np.array(getlist).astype(int):
                ofile.write("QM\n")
                continue
            match = re.search(r"^\$end", line, flags=re.MULTILINE)
            if match:

                for i in range(len(m2list), len(new_field)):
                    ofile.write(
                        "{:<.10f} {:<.10f} {:<.10f} {:<.10f}\n".format(
                            float(new_field[i][0]),
                            float(new_field[i][1]),
                            float(new_field[i][2]),
                            float(new_field[i][3]),
                        )
                    )
                ofile.write("$end\n")
                continue
            else:
                ofile.write(line)
    ofile.close()


def read_xyzq(inp, qline):
    xyzq = []
    with open(inp) as ifile:
        count = 0
        for line in ifile:
            count += 1
            if int(count) == int(qline):
                match = re.search(
                    r"^\s*([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)\s+([-]*\d+\.\d+)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    print("Line " + str(qline) + " does not contain data. Exiting.")
                    exit(1)
                xyzq.extend(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                        float(match.group(4)),
                    ]
                )
                break
    return xyzq


def get_qmlist(inp):
    xyz = []
    with open(inp) as ifile:
        for line in ifile:
            match = re.search(
                r"^\s*([-]*\d+\.\d+[e]*[-]*[+]*\d*)\s+([-]*\d+\.\d+[e]*[-]*[+]*\d*)\s+([-]*\d+\.\d+[e]*[-]*[+]*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                xyz.append(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                    ]
                )
    return xyz


def get_m2vec(line, inp):
    m2list = re.findall("\d+", line)
    m2coordsq = []
    for i in range(0, len(m2list)):
        m2coordsq.append(read_xyzq(inp, m2list[i]))
    return m2list, m2coordsq


def get_m2vec_fieldsonly(m2entry, pcf):
    m2coordsq = []
    for i in range(0, len(m2entry)):
        m2coordsq.append(pcf[int(m2entry[i] - 1)])
    return m2coordsq


def make_droplist(inp, m1line, m2list):
    droplist = []
    count = 0
    with open(inp) as ifile:
        for line in ifile:
            count += 1
            if count == int(m1line) or count in m2list:
                continue
            else:
                droplist.append(count)
    return droplist


def generate_charge_shift_fieldsonly(pcf, m1list, qmcoords, m2list, jobname, basedir):

    getlist = m1list
    orgfield = []
    qmlist = qmcoords
    target_sum = np.array([0.0, 0.0, 0.0])
    for i in range(0, len(getlist)):
        orgcoordsq = pcf[int(getlist[i]) - 1]

        orgfield.append(orgcoordsq)

    for element in qmlist:
        target_sum += np.array(
            make_p_sum.sum_pcf_tm_nofile(orgfield, element[0], element[1], element[2])
        )

    m1coordsq = []
    m2_nolist = m2list
    m2coordsqlist = []
    dispcharge = 0.1944  # for 0.023333333333 shifted charge on three M2 atoms; should work in any case unless distances become HUGE
    disp_charge_vec = []
    count = 0
    for element in m2list:
        m1 = pcf[int(getlist[count]) - 1]
        m1coordsq.append(m1)

        m2coordsqthing = get_m2vec_fieldsonly(element, pcf)
        m2coordsq = []
        for i in range(0, len(m2coordsqthing)):
            m2coordsq.append(list(m2coordsqthing[i]))

        if m1 == [] or m2coordsq == []:
            print("Ding")
            continue
        for entry in m2coordsq:

            entry[3] = m1[3] / float(len(m2coordsq))

        effect_dispcharge = float(dispcharge) * (
            m1[3] / float(len(element)) / 0.023333333333
        )
        if math.sqrt(float(effect_dispcharge) * float(effect_dispcharge)) > 0.5:
            if float(effect_dispcharge) < 0.0:
                effect_dispcharge = -0.5
            else:
                effect_dispcharge = 0.5
        disp_charge_vec.append(effect_dispcharge)

        count += 1
        m2coordsqlist.append(m2coordsq)

    disp = 0.1
    disp_vec = []
    for i in range(0, len(m2_nolist)):
        for j in range(0, len(m2_nolist[i])):
            disp_vec.append(disp)

    corr_charge_list = create_corr_charges(
        m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist
    )

    new_field = np.array(make_new_field(m2coordsqlist, corr_charge_list))
    curr_sum = np.array([0.0, 0.0, 0.0])
    for element in qmlist:
        curr_sum += np.array(
            make_p_sum.sum_pcf_tm_nofile(new_field, element[0], element[1], element[2])
        )

    curr_delta = math.sqrt(
        (curr_sum[0] - target_sum[0]) * (curr_sum[0] - target_sum[0])
        + (curr_sum[1] - target_sum[1]) * (curr_sum[1] - target_sum[1])
        + (curr_sum[2] - target_sum[2]) * (curr_sum[2] - target_sum[2])
    )
    max_change_disp = disp
    change_disp = []
    change_disp = np.array(disp_vec)
    startfac = 10.0
    v_opt_tracker = []
    for element in disp_vec:
        changestat = False
        point_vals = [0.0, 0.0, changestat]
        v_opt_tracker.append(point_vals)
    iter_count = 0
    strike = 0

    while max_change_disp > 0.0000001 and iter_count < 200:

        max_change_disp = 0.0
        iter_count += 1
        count = 0
        mod_disp_vec = np.array(disp_vec)
        for i in range(0, len(m2_nolist)):
            for j in range(0, len(m2_nolist[i])):
                if not v_opt_tracker[count][2]:
                    change_disp[count] /= startfac
                    v_opt_tracker[count][2] = True
                if change_disp[count] <= 0.00001:
                    count += 1
                    continue
                new_disp_vec = np.array(disp_vec)
                new_disp_vec[count] += change_disp[count]
                curr_corr_charge_list = create_corr_charges(
                    m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, m2_nolist
                )
                curr_field = np.array(
                    make_new_field(m2coordsqlist, curr_corr_charge_list)
                )
                new_sum = np.array([0.0, 0.0, 0.0])
                for element in qmlist:
                    new_sum += np.array(
                        make_p_sum.sum_pcf_tm_nofile(
                            curr_field, element[0], element[1], element[2]
                        )
                    )
                v_opt_tracker[count][0] = math.sqrt(
                    (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                    + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                    + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
                )
                new_disp_vec = np.array(disp_vec)
                new_disp_vec[count] -= change_disp[count]
                curr_corr_charge_list = create_corr_charges(
                    m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, m2_nolist
                )
                curr_field = np.array(
                    make_new_field(m2coordsqlist, curr_corr_charge_list)
                )
                new_sum = np.array([0.0, 0.0, 0.0])
                for element in qmlist:
                    new_sum += np.array(
                        make_p_sum.sum_pcf_tm_nofile(
                            curr_field, element[0], element[1], element[2]
                        )
                    )
                v_opt_tracker[count][1] = math.sqrt(
                    (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                    + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                    + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
                )
                if (float(v_opt_tracker[count][1]) < float(curr_delta)) and (
                    float(v_opt_tracker[count][1]) < float(v_opt_tracker[count][0])
                ):
                    mod_disp_vec[count] -= change_disp[count]
                    if (
                        math.sqrt(mod_disp_vec[count] * mod_disp_vec[count])
                    ) < 0.000001:
                        mod_disp_vec[count] = 0.000001
                    v_opt_tracker[count][2] = True
                elif (float(v_opt_tracker[count][0]) < float(curr_delta)) and (
                    float(v_opt_tracker[count][0]) < float(v_opt_tracker[count][1])
                ):
                    mod_disp_vec[count] += change_disp[count]
                    if (
                        math.sqrt(mod_disp_vec[count] * mod_disp_vec[count])
                    ) < 0.000001:
                        mod_disp_vec[count] = 0.000001
                    v_opt_tracker[count][2] = True
                else:
                    v_opt_tracker[count][2] = False
                count += 1
        for element in mod_disp_vec:
            if element > max_change_disp:
                max_change_disp = element
        new_disp_vec = np.array(mod_disp_vec)
        curr_corr_charge_list = create_corr_charges(
            m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, m2_nolist
        )
        curr_field = np.array(make_new_field(m2coordsqlist, curr_corr_charge_list))
        new_sum = np.array([0.0, 0.0, 0.0])
        for element in qmlist:
            new_sum += np.array(
                make_p_sum.sum_pcf_tm_nofile(
                    curr_field, element[0], element[1], element[2]
                )
            )
        new_delta = math.sqrt(
            (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
            + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
            + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
        )
        if new_delta < curr_delta:
            disp_vec = new_disp_vec
            curr_delta = new_delta
            curr_sum = new_sum
            strike = 0
        elif new_delta == curr_delta:
            disp_vec = new_disp_vec
            curr_delta = new_delta
            curr_sum = new_sum
            strike += 1
            if strike > 3:
                max_change_disp = 0.0000001
                break
        else:
            count = 0
            for i in range(0, len(m2_nolist)):
                for j in range(0, len(m2_nolist[i])):
                    change_disp[count] /= startfac
                    count += 1
        if max_change_disp < 0.0000001:
            break
    corr_charge_list = create_corr_charges(
        m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist
    )
    new_field = np.array(make_new_field(m2coordsqlist, corr_charge_list))
    new_sum = np.array([0.0, 0.0, 0.0])
    for element in qmlist:
        new_sum += np.array(
            make_p_sum.sum_pcf_tm_nofile(new_field, element[0], element[1], element[2])
        )
    new_delta = math.sqrt(
        (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
        + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
        + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
    )
    outname = jobname + ".pointcharges"

    write_new_field_to_disk_listsonly(pcf, outname, new_field, getlist, m2_nolist)


def generate_charge_shift(syscmds):
    basedir = os.path.dirname(os.path.abspath(__file__))
    inp = syscmds[1]
    m1list = syscmds[2]
    qmcoords = syscmds[3]
    m2listlist = syscmds[4]
    ofilename = syscmds[5]
    getlist = []
    with open(m1list) as ifile:
        for line in ifile:
            getlist.extend([int(line)])
    orgfield = []

    for i in range(0, len(getlist)):
        orgcoordsq = read_xyzq(inp, getlist[i])

        orgfield.append(orgcoordsq)
    qmlist = get_qmlist(qmcoords)
    target_sum = np.array([0.0, 0.0, 0.0])

    for element in qmlist:
        target_sum += np.array(
            make_p_sum.sum_pcf_tm_nofile(orgfield, element[0], element[1], element[2])
        )

    m1coordsq = []
    m2_nolist = []
    m2coordsqlist = []
    dispcharge = 0.1944  # for 0.023333333333 shifted charge on three M2 atoms; should work in any case unless distances become HUGE

    disp_charge_vec = []
    count = 0
    with open(m2listlist) as ifile:
        for line in ifile:

            m1 = read_xyzq(inp, getlist[count])
            m1coordsq.append(m1)
            m2_no, m2coordsq = get_m2vec(line, inp)
            if m1 == [] or m2coordsq == []:
                continue
            for element in m2coordsq:
                element[3] = m1[3] / float(len(m2_no))
            disp_charge_vec.append(
                float(dispcharge) * (m1[3] / float(len(m2_no)) / 0.023333333333)
            )
            count += 1
            m2_nolist.append(m2_no)
            m2coordsqlist.append(m2coordsq)

    # initialize displacement
    disp = 0.1
    disp_vec = []
    for i in range(0, len(m2_nolist)):
        for j in range(0, len(m2_nolist[i])):
            disp_vec.append(disp)

    corr_charge_list = create_corr_charges(
        m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist
    )

    new_field = np.array(make_new_field(m2coordsqlist, corr_charge_list))

    curr_sum = np.array([0.0, 0.0, 0.0])
    for element in qmlist:
        curr_sum += np.array(
            make_p_sum.sum_pcf_tm_nofile(new_field, element[0], element[1], element[2])
        )

    curr_delta = math.sqrt(
        (curr_sum[0] - target_sum[0]) * (curr_sum[0] - target_sum[0])
        + (curr_sum[1] - target_sum[1]) * (curr_sum[1] - target_sum[1])
        + (curr_sum[2] - target_sum[2]) * (curr_sum[2] - target_sum[2])
    )
    # opt distances
    max_change_disp = disp
    change_disp = []
    change_disp = np.array(disp_vec)
    startfac = 10.0
    v_opt_tracker = []
    for element in disp_vec:
        changestat = False
        point_vals = [0.0, 0.0, changestat]
        v_opt_tracker.append(point_vals)
    iter_count = 0
    strike = 0
    while max_change_disp > 0.0000001 and iter_count < 200:

        max_change_disp = 0.0
        iter_count += 1
        count = 0
        mod_disp_vec = np.array(disp_vec)
        for i in range(0, len(m2_nolist)):
            for j in range(0, len(m2_nolist[i])):
                if not v_opt_tracker[count][2]:
                    change_disp[count] /= startfac
                    v_opt_tracker[count][2] = True
                if change_disp[count] <= 0.00001:
                    count += 1
                    continue
                new_disp_vec = np.array(disp_vec)
                new_disp_vec[count] += change_disp[count]
                curr_corr_charge_list = create_corr_charges(
                    m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, m2_nolist
                )
                curr_field = np.array(
                    make_new_field(m2coordsqlist, curr_corr_charge_list)
                )
                new_sum = np.array([0.0, 0.0, 0.0])
                for element in qmlist:
                    new_sum += np.array(
                        make_p_sum.sum_pcf_tm_nofile(
                            curr_field, element[0], element[1], element[2]
                        )
                    )
                v_opt_tracker[count][0] = math.sqrt(
                    (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                    + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                    + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
                )
                new_disp_vec = np.array(disp_vec)
                new_disp_vec[count] -= change_disp[count]
                curr_corr_charge_list = create_corr_charges(
                    m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, m2_nolist
                )
                curr_field = np.array(
                    make_new_field(m2coordsqlist, curr_corr_charge_list)
                )
                new_sum = np.array([0.0, 0.0, 0.0])
                for element in qmlist:
                    new_sum += np.array(
                        make_p_sum.sum_pcf_tm_nofile(
                            curr_field, element[0], element[1], element[2]
                        )
                    )
                v_opt_tracker[count][1] = math.sqrt(
                    (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
                    + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
                    + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
                )
                if (float(v_opt_tracker[count][1]) < float(curr_delta)) and (
                    float(v_opt_tracker[count][1]) < float(v_opt_tracker[count][0])
                ):
                    mod_disp_vec[count] -= change_disp[count]
                    v_opt_tracker[count][2] = True
                elif (float(v_opt_tracker[count][0]) < float(curr_delta)) and (
                    float(v_opt_tracker[count][0]) < float(v_opt_tracker[count][1])
                ):
                    mod_disp_vec[count] += change_disp[count]
                    v_opt_tracker[count][2] = True
                else:
                    v_opt_tracker[count][2] = False
                count += 1
        for element in mod_disp_vec:
            if element > max_change_disp:
                max_change_disp = element
        new_disp_vec = np.array(mod_disp_vec)
        curr_corr_charge_list = create_corr_charges(
            m1coordsq, m2coordsqlist, new_disp_vec, disp_charge_vec, m2_nolist
        )
        curr_field = np.array(make_new_field(m2coordsqlist, curr_corr_charge_list))
        new_sum = np.array([0.0, 0.0, 0.0])
        for element in qmlist:
            new_sum += np.array(
                make_p_sum.sum_pcf_tm_nofile(
                    curr_field, element[0], element[1], element[2]
                )
            )
        new_delta = math.sqrt(
            (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
            + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
            + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
        )

        if new_delta < curr_delta:
            disp_vec = new_disp_vec
            curr_delta = new_delta
            curr_sum = new_sum
            strike = 0
        elif new_delta == curr_delta:
            disp_vec = new_disp_vec
            curr_delta = new_delta
            curr_sum = new_sum
            strike += 1
            if strike > 3:
                max_change_disp = 0.0000001
                break
        else:
            count = 0
            for i in range(0, len(m2_nolist)):
                for j in range(0, len(m2_nolist[i])):
                    change_disp[count] /= startfac
                    count += 1
        if max_change_disp < 0.0000001:
            break
    corr_charge_list = create_corr_charges(
        m1coordsq, m2coordsqlist, disp_vec, disp_charge_vec, m2_nolist
    )
    new_field = np.array(make_new_field(m2coordsqlist, corr_charge_list))
    new_sum = np.array([0.0, 0.0, 0.0])
    for element in qmlist:
        new_sum += np.array(
            make_p_sum.sum_pcf_tm_nofile(new_field, element[0], element[1], element[2])
        )
    new_delta = math.sqrt(
        (new_sum[0] - target_sum[0]) * (new_sum[0] - target_sum[0])
        + (new_sum[1] - target_sum[1]) * (new_sum[1] - target_sum[1])
        + (new_sum[2] - target_sum[2]) * (new_sum[2] - target_sum[2])
    )
    write_new_field_to_disk(inp, ofilename, new_field, getlist, m2_nolist)


if __name__ == "__main__":
    generate_charge_shift(sys.argv)
