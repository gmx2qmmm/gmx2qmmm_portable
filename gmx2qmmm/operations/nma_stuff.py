# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "jangoetze"
__date__ = "$06-Feb-2018 12:45:17$"

import os
import subprocess

import numpy as np

from gmx2qmmm._helper import _flatten
from gmx2qmmm.operations import hes_xyz_g09RevD_01_fchk as write_hess


def make_disp_g96_name(gro, curr_coord, curr_atom, disp):
    new_g96_name = gro + "." + str(curr_atom) + "at"
    if disp < 0.0:
        new_g96_name += "Minus"
    else:
        new_g96_name += "Plus"
    new_g96_name += curr_coord + "ofAtom" + str(curr_atom) + ".g96"
    return new_g96_name


def make_disp_g96(gro, coord_index, curr_atom, disp):
    coord_list = ["X", "Y", "Z"]
    new_g96_name = make_disp_g96_name(gro, coord_list[coord_index], curr_atom, disp)
    nm_disp = float(disp) * 0.052917721
    with open(new_g96_name, "w") as ofile:
        with open(gro) as ifile:
            count = 0
            for line in ifile:
                ofile.write(line)
                count += 1
                if count == 4:
                    break
            count = 0
            for line in ifile:
                ofile.write(line)
                count += 1
                if count == int(curr_atom) - 1:
                    break
            for line in ifile:
                match = re.search(
                    r"^(.{24})\s*([-]*\d+\.\d+)\s*([-]*\d+\.\d+)\s*([-]*\d+\.\d+)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    print("Error in automatically generated file " + str(
                        gro
                    ) + " for atom " + str(
                        curr_atom
                    ) + ". This should never happen. Exiting. Last line:")
                    print(line)
                    exit(1)
                ofile.write(match.group(1))
                for i in range(0, 3):
                    if i != coord_index:
                        ofile.write("{:>16.9f}".format(float(match.group(i + 2))))
                    else:
                        ofile.write(
                            "{:>16.9f}".format(float(match.group(i + 2)) + nm_disp)
                        )
                ofile.write("\n")
                break
            for line in ifile:
                ofile.write(line)
    return new_g96_name


def eval_gro(
    gro,
    top,
    active,
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
    linkcorrlist,
    flaglist,
    pcffile,
    linkatoms,
    logfile,
    basedir,
):

    jobname = gro[:-4]
    new_pcffile = str(jobname + ".pointcharges")
    archivename = str(jobname) + ".tar.gz"
    if os.path.isfile(archivename):
        subprocess.call(["tar", "-xf", archivename])
        subprocess.call(["rm", archivename])
    new_xyzq, m1list, m2list, new_links = qmmm_prep(
        gro, top, jobname, 0, qminfo, qmatomlist, connlist, linkatoms, basedir, logfile
    )
    # we may skip a new pcf generation if all act atoms are in qm region - for this case the pcf stays the same
    qmenergy, mmenergy, qm_corrdata = get_energy(
        gro,
        jobname,
        new_xyzq,
        connlist,
        qmatomlist,
        m1list,
        m2list,
        q1list,
        qmmmtop,
        qminfo,
        mminfo,
        qmmminfo,
        linkcorrlist,
        flaglist,
        str(jobname + ".pointcharges"),
        new_links,
        0,
        logfile,
        basedir,
    )
    sp_energy = float(qmenergy) + float(mmenergy)
    curr_forces = read_forces(
        qmatomlist,
        m1list,
        qmmmtop,
        qminfo,
        jobname,
        0,
        logfile,
        linkcorrlist,
        new_xyzq,
        new_pcffile,
        qm_corrdata,
        m2list,
        q1list,
        new_links,
        active,
        basedir,
    )
    archive = ["tar", "-cf", str(jobname) + ".tar"]
    files = [
        new_pcffile,
        str(jobname) + ".edr",
        str(jobname) + ".edr.xvg",
        str(jobname) + ".trr",
        str(jobname) + ".xvg",
        str(jobname) + ".gmx.log",
        str(jobname) + ".g96",
        str(jobname) + ".gjf",
        str(jobname) + ".mdp",
        str(jobname) + ".boxlarge.g96",
        str(jobname) + ".tpr",
        str(jobname) + ".gjf.log",
        str(jobname) + ".fort.7",
    ]
    archive.extend(files)
    subprocess.call(archive)
    subprocess.call(["gzip", str(jobname) + ".tar"])
    delete = ["rm"]
    delete.extend(files)
    subprocess.call(delete)
    return sp_energy, curr_forces


def prepare_hess(hessian_full, active):
    red_hess = []
    for element in hessian_full:
        red_hess_line = []
        for i in range(0, len(element)):
            if i + 1 in np.array(active).astype(int):
                red_hess_line.extend(element[i])
        red_hess.append(red_hess_line)
    return red_hess


def log_nma(qmmminfo, logfile, evals, nm_matrix, active, qmmmtop, xyzq, hess):
    logger(logfile, "------Results of normal mode analysis------\n")
    jobname = str(qmmminfo[0] + "_final")
    write_pseudofchk_file(
        jobname, evals, nm_matrix, hess, active, qmmmtop, logfile, xyzq
    )
    logger(
        logfile,
        "Wrote new Gaussian03 style pseudo-fchk file (" + jobname + ".pseudofchk).\n",
    )
    write_hess.hes_xyz_fchk(str(jobname + ".pseudofchk"), str(jobname + ".hess"))
    logger(logfile, "Wrote orca format " + jobname + ".hess file.\n")


def write_pseudofchk_file(
    jobname, evals, nm_matrix, hess, active, qmmmtop, logfile, xyzq
):
    flat_nm = list(_flatten(nm_matrix))
    outname = jobname + ".pseudofchk"
    n_a = len(active)
    dof = len(active) * 3
    atom_map = {
        "1": "H",
        "2": "He",
        "3": "Li",
        "4": "Be",
        "5": "B",
        "6": "C",
        "7": "N",
        "8": "O",
        "9": "F",
        "10": "Ne",
        "11": "Na",
        "12": "Mg",
        "13": "Al",
        "14": "Si",
        "15": "P",
        "16": "S",
        "17": "Cl",
        "18": "Ar",
        "19": "K",
        "20": "Ca",
        "21": "Sc",
        "22": "Ti",
        "23": "V",
        "24": "Cr",
        "25": "Mn",
        "26": "Fe",
        "27": "Co",
        "28": "Ni",
        "29": "Cu",
        "30": "Zn",
        "31": "Ga",
        "32": "Ge",
        "33": "As",
        "34": "Se",
        "35": "Br",
        "36": "Kr",
        "37": "Rb",
        "38": "Sr",
        "39": "Y",
        "40": "Zr",
        "41": "Nb",
        "42": "Mo",
        "43": "Tc",
        "44": "Ru",
        "45": "Rh",
        "46": "Pd",
        "47": "Ag",
        "48": "Cd",
        "49": "In",
        "50": "Sn",
        "51": "Sb",
        "52": "Te",
        "53": "I",
        "54": "Xe",
        "55": "Cs",
        "56": "Ba",
        "57": "La",
        "58": "Ce",
        "59": "Pr",
        "60": "Nd",
        "61": "Pm",
        "62": "Sm",
        "63": "Eu",
        "64": "Gd",
        "65": "Tb",
        "66": "Dy",
        "67": "Ho",
        "68": "Er",
        "69": "Tm",
        "70": "Yb",
        "71": "Lu",
        "72": "Hf",
        "73": "Ta",
        "74": "W",
        "75": "Re",
        "76": "Os",
        "77": "Ir",
        "78": "Pt",
        "79": "Au",
        "80": "Hg",
        "81": "Tl",
        "82": "Pb",
        "83": "Bi",
        "84": "Po",
        "85": "At",
        "86": "Rn",
        "87": "Fr",
        "88": "Ra",
        "89": "Ac",
        "90": "Th",
        "91": "Pa",
        "92": "U",
        "93": "Np",
        "94": "Pu",
        "95": "Am",
        "96": "Cm",
        "97": "Bk",
        "98": "Cf",
        "99": "Es",
        "100": "Fm",
        "101": "Md",
        "102": "No",
        "103": "Lr",
        "104": "Rf",
        "105": "Db",
        "106": "Sg",
        "107": "Bh",
        "108": "Hs",
        "109": "Mt",
        "110": "Ds",
        "111": "Rg",
        "112": "Cn",
        "113": "Nh",
        "114": "Fl",
        "115": "Mc",
        "116": "Lv",
        "117": "Ts",
        "118": "Og",
    }
    number_map = {value: key for key, value in atom_map.items()}
    mass_map = {
        "H": "1.008",
        "He": "4.0026",
        "Li": "6.94",
        "Be": "9.0122",
        "B": "10.81",
        "C": "12.011",
        "N": "14.007",
        "O": "15.999",
        "F": "18.998",
        "Ne": "20.180",
        "Na": "22.990",
        "Mg": "24.305",
        "Al": "26.982",
        "Si": "28.085",
        "P": "30.974",
        "S": "32.06",
        "Cl": "35.45",
        "Ar": "39.948",
        "K": "39.098",
        "Ca": "40.0784",
        "Sc": "44.956",
        "Ti": "47.867",
        "V": "50.942",
        "Cr": "51.996",
        "Mn": "54.938",
        "Fe": "55.8452",
        "Co": "58.933",
        "Ni": "58.693",
        "Cu": "63.5463",
        "Zn": "65.382",
        "Ga": "69.723",
        "Ge": "72.6308",
        "As": "74.922",
        "Se": "78.9718",
        "Br": "79.904",
        "Kr": "83.7982",
        "Rb": "85.468",
        "Sr": "87.62",
        "Y": "88.906",
        "Zr": "91.2242",
        "Nb": "92.906",
        "Mo": "95.95",
        "Tc": "98.906254721",
        "Ru": "101.072",
        "Rh": "102.91",
        "Pd": "106.42",
        "Ag": "107.87",
        "Cd": "112.41",
        "In": "114.82",
        "Sn": "118.71",
        "Sb": "121.76",
        "Te": "127.603",
        "I": "126.90",
        "Xe": "131.29",
        "Cs": "132.91",
        "Ba": "137.33",
        "La": "138.91",
        "Ce": "140.12",
        "Pr": "140.91",
        "Nd": "144.24",
        "Pm": "144.9127493",
        "Sm": "150.362",
        "Eu": "151.96",
        "Gd": "157.253",
        "Tb": "158.93",
        "Dy": "162.50",
        "Ho": "164.93",
        "Er": "167.26",
        "Tm": "168.93",
        "Yb": "173.05",
        "Lu": "174.97",
        "Hf": "178.492",
        "Ta": "180.95",
        "W": "183.84",
        "Re": "186.21",
        "Os": "190.233",
        "Ir": "192.22",
        "Pt": "195.08",
        "Au": "196.97",
        "Hg": "200.59",
        "Tl": "204.38",
        "Pb": "207.2",
        "Bi": "208.98",
        "Po": "208.982430420",
        "At": "209.9871488",
        "Rn": "222.017577725",
        "Fr": "223.019735926",
        "Ra": "226.025409825",
        "Ac": "227.027752126",
        "Th": "232.04",
        "Pa": "231.04",
        "U": "238.03",
        "Np": "237.04817342",
        "Pu": "244.0642045",
        "Am": "243.061381125",
        "Cm": "247.0703545",
        "Bk": "247.0703076",
        "Cf": "251.0795875",
        "Es": "252.082985",
        "Fm": "257.0951067",
        "Md": "258.0984315",
        "No": "259.1010311",
        "Lr": "266.1198356",
        "Rf": "267.1217962",
        "Db": "268.1256757",
        "Sg": "269.1286339",
        "Bh": "270.1333631",
        "Hs": "277.1519058",
        "Mt": "278.1563168",
        "Ds": "281.1645159",
        "Rg": "282.1691272",
        "Cn": "285.177126",
        "Nh": "286.1822172",
        "Fl": "289.190426",
        "Mc": "289.1936389",
        "Lv": "293.204496",
        "Ts": "294.2104674",
        "Og": "295.2162469",
    }
    atoms = get_atoms(qmmmtop, logfile)
    with open(outname, "w") as ofile:
        ofile.write("Atomic numbers N= " + str(n_a) + "\n")
        linecount = 0
        for i in range(0, n_a):
            linecount += 1
            if linecount > 6:
                ofile.write("\n")
                linecount = 1
            ofile.write("{:>12d}".format(int(number_map[atoms[active[i]]])))
        ofile.write("\nCurrent cartesian coordinates N= " + str(dof) + "\n")
        linecount = 0
        for i in range(0, n_a):
            for j in range(0, 3):
                linecount += 1
                if linecount > 5:
                    ofile.write("\n")
                    linecount = 1
                ofile.write(" {:>16.8E}".format(float(xyzq[int(active[i]) - 1][j])))
        ofile.write(
            "\nCartesian Force Constants N= " + str((dof * dof + dof) / 2) + "\n"
        )
        linecount = 0
        for i in range(0, dof):
            for j in range(0, i + 1):
                linecount += 1
                if linecount > 5:
                    ofile.write("\n")
                    linecount = 1
                curr_element = (float(hess[i][j]) + float(hess[j][i])) / 2.0
                ofile.write(" {:>16.8E}".format(float(curr_element)))
        ofile.write("\nVib-E2 N= " + str((dof - 6) * 14) + "\n")
        linecount = 0
        for i in range(0, (dof - 6) * 14):
            linecount += 1
            if linecount > 5:
                ofile.write("\n")
                linecount = 1
            if i >= len(evals):
                ofile.write(" {:>16.8E}".format(float(1000.0)))
            else:
                if float(evals[i]) < 0.0:
                    ofile.write(
                        " {:>16.8E}".format(
                            float(-1.0 * sqr(abs(float(evals[i])) * float(26424608)))
                        )
                    )
                else:
                    ofile.write(
                        " {:>16.8E}".format(
                            float(sqr(float(evals[i]) * float(26424608)))
                        )
                    )
        ofile.write("\nVib-Modes N= " + str((dof - 6) * dof) + "\n")
        linecount = 0
        for i in range(0, (dof - 6) * dof):
            linecount += 1
            if linecount > 5:
                ofile.write("\n")
                linecount = 1
            if i >= len(flat_nm):
                ofile.write(" {:>16.8E}".format(float(1000.0)))
            else:
                ofile.write(" {:>16.8E}".format(float(flat_nm[i])))
        ofile.write("\nVib-AtMass N= " + str(n_a) + "\n")
        linecount = 0
        for i in range(0, n_a):
            linecount += 1
            if linecount > 5:
                ofile.write("\n")
                linecount = 1
            ofile.write(" {:>16.8E}".format(float(mass_map[atoms[active[i]]])))
        ofile.write("\nDipole Derivatives N= " + str(9 * n_a) + "\n")
        linecount = 0
        for i in range(0, 9 * n_a):
            linecount += 1
            if linecount > 5:
                ofile.write("\n")
                linecount = 1
            ofile.write(" {:>16.8E}".format(float(1000.0)))
        ofile.write("\n")


def get_xyz_2nd_deriv(
    curr_atom,
    start_energy,
    start_forces,
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
    linkcorrlist,
    flaglist,
    pcffile,
    linkatoms,
    active,
    logfile,
    basedir,
):
    # DO FOR 3 COORDS AND PLUS AND MINUS
    disp = float(qmmminfo[6])
    grad_vectors = []
    logger(
        logfile,
        "------Getting second derivatives for atom " + str(curr_atom) + "------\n",
    )
    for i in range(0, 3):
        new_g96 = make_disp_g96(gro, i, curr_atom, disp)
        plus_energy, plus_forces = eval_gro(
            new_g96,
            top,
            active,
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
            linkcorrlist,
            flaglist,
            pcffile,
            linkatoms,
            logfile,
            basedir,
        )
        plus_grad = np.array(plus_forces) * -1.0
        if plus_energy < start_energy:
            logger(
                logfile,
                "Note: Atom "
                + str(curr_atom)
                + " displaced by "
                + str(disp)
                + " a.u. along coordinate "
                + str(i)
                + " yields lower energy than initial point:\n",
            )
            logger(
                logfile,
                str(float(plus_energy))
                + " a.u. vs. "
                + str(float(start_energy))
                + " a.u.\n",
            )
        new_g96 = make_disp_g96(gro, i, curr_atom, -1.0 * disp)
        minus_energy, minus_forces = eval_gro(
            new_g96,
            top,
            active,
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
            linkcorrlist,
            flaglist,
            pcffile,
            linkatoms,
            logfile,
            basedir,
        )
        minus_grad = np.array(minus_forces) * -1.0
        if minus_energy < start_energy:
            logger(
                logfile,
                "Note: Atom "
                + str(curr_atom)
                + " displaced by -"
                + str(disp)
                + " a.u. along coordinate "
                + str(i)
                + " yields lower energy than initial point:\n",
            )
            logger(
                logfile,
                str(float(minus_energy))
                + " a.u. vs. "
                + str(float(start_energy))
                + " a.u.\n",
            )
        average_grad = list(_flatten(
            (np.array(list(_flatten(plus_grad))) - np.array(list(_flatten(minus_grad)))) / (float(disp) * 2.0)
        ))
        grad_vectors.append(average_grad)
    return grad_vectors


if __name__ == "__main__":
    print("This is a library for normal mode analysis functions. Do not execute directly.")
